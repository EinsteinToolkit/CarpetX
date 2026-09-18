#include "driver.hxx"
#include "fillpatch.hxx"
#include "io.hxx"
#include "loop.hxx"
#include "schedule.hxx"
#include "task_manager.hxx"
#include "timer.hxx"
#include "valid.hxx"

#include <cctk.h>
#include <cctk_Parameters.h>
#include <cctk_Schedule.h>
#include <cctki_GHExtensions.h>
#include <cctki_ScheduleBindings.h>
#include <cctki_WarnLevel.h>

#include <AMReX_MultiFabUtil.H>

#if defined _OPENMP
#include <omp.h>
#elif defined __HIPCC__
#define omp_get_max_threads() 1
#define omp_get_num_threads() 1
#define omp_get_thread_num() 0
#define omp_in_parallel() 0
#else
static inline int omp_get_max_threads() { return 1; }
static inline int omp_get_num_threads() { return 1; }
static inline int omp_get_thread_num() { return 0; }
static inline int omp_in_parallel() { return 0; }
#endif

#include <sys/time.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace CarpetX {

#ifndef CCTK_HAVE_CGH_LEVEL
#error                                                                         \
    "The Cactus flesh does not support cctk_level in the cGH structure. Update the flesh."
#endif
#ifndef CCTK_HAVE_CGH_PATCH
#error                                                                         \
    "The Cactus flesh does not support cctk_patch in the cGH structure. Update the flesh."
#endif
#ifndef CCTK_HAVE_CGH_COMPONENT
#error                                                                         \
    "The Cactus flesh does not support cctk_component in the cGH structure. Update the flesh."
#endif
#ifndef CCTK_HAVE_CGH_TILE
#error                                                                         \
    "The Cactus flesh does not support cctk_tile_min etc. in the cGH structure. Update the flesh."
#endif

#if defined _OPENMP
#if !defined AMREX_USE_OMP
#error                                                                         \
    "Cactus is configured with OpenMP, but AMReX is configured without OpenMP. This does not work."
#endif
#endif

namespace {
double gettime() {
  timeval tv;
  gettimeofday(&tv, nullptr);
  return tv.tv_sec + tv.tv_usec / 1.0e+6;
}
} // namespace

// Used to pass active levels from AMReX's regridding functions
std::optional<active_levels_t> active_levels;

////////////////////////////////////////////////////////////////////////////////
//
// AMR-C7.  THE SCHEDULE-ORDER INSTRUMENT (`CARPETX_LOG_SCHED`).
//
// WHAT IT IS FOR.  Two sentences this project has been making from a code
// reading and wants to make from a measurement instead: that `Restrict` runs
// BEFORE the interpatch fill inside a sync, and that `Reflux` moves nothing.
// `[P186]`: a deadness claim backed by one grep is one typo away from being
// wrong, and an ordering claim read off the source is one refactor away from
// it.
//
// WHY A SEQUENCE NUMBER AND NOT THE ITERATION.  `[P422]`: an instrument that
// delimits its stream by `cctk_iteration` MERGES the two traversals that both
// live inside `ScheduleTraverseGH iteration 0`, which is exactly where the
// initial restriction and the initial interpatch fills are.  So the clock here
// is a single monotonic counter shared by every site below, and the iteration
// is printed as CONTEXT beside it, never as a delimiter.  Two events are
// ordered iff their `seq` values are.
//
// WHY STDERR.  `[P417]`: the flesh `freopen`s stdout to the null device on
// every non-root process (`flesh/src/main/CommandLine.c:783`) unless the run is
// given `-r`, so no `CCTK_VINFO` is a per-process channel and an ordering claim
// read off INFO lines is a claim about rank 0 and nothing else.  stderr the
// flesh leaves alone.  This is the same reasoning, and the same mechanism, as
// `CapyrX_MultiPatch`'s `CAPYRX_LOG_LEVELS` instrument, whose `MPLEVELS` line
// marks the interpatch fill from the OTHER repo and is the independent second
// reading of the same ordering.
//
// COST WHEN OFF.  One `std::getenv` per process and one predictable branch per
// site.  No line, no string, no allocation.  It is ALWAYS COMPILED rather than
// `CCTK_DEBUG`-gated for AMR-A4's reason: the numbers have to be readable in
// the optimized build, which is the only build in which the production geometry
// runs in minutes.
//
// CONCURRENCY.  Every site below is reached in global mode -- `SyncGroupsByDirI`
// asserts `in_global_mode`, and `Restrict`/`Reflux` are called only from the
// driver's own serial points -- so the counter is not shared between threads.
// Per B10's rule the static is written only on the path that also reads it,
// i.e. only when the instrument is on, so a run without it is untouched.
//
// THE LINE IS ASSEMBLED AND EMITTED WITH ONE `<<`, so that a rank's own stream
// carries whole lines in program order.  That order IS the measurement.
//
////////////////////////////////////////////////////////////////////////////////

namespace {

bool log_sched_on() {
  static const bool on = std::getenv("CARPETX_LOG_SCHED") != nullptr;
  return on;
}

void log_sched(const cGH *const cctkGH, const std::string &fields) {
  // `proc` is on every line because at more than one rank the launcher MERGES
  // the per-process stderr streams into one file, and an ordering claim read
  // off an interleaved stream is not a claim about any process.  With it, the
  // reader splits exactly; without it, it would have to guess from the `seq`
  // values, and a single dropped line would make the guess wrong silently.
  static const int myproc = CCTK_MyProc(cctkGH);
  static long seq = 0;
  std::ostringstream line;
  line << "SCHED seq=" << ++seq
       << " it=" << (cctkGH ? cctkGH->cctk_iteration : -1)
       << " proc=" << myproc << " " << fields << "\n";
  std::cerr << line.str();
}

// THE FLUX CENSUS, and it is the driver's own reader rather than a grep.
//
// A group gets a flux register -- and is therefore the only kind of group
// `Reflux` can move -- iff its TAGS carry `fluxes="<gx> <gy> <gz>"`
// (`driver.cxx`, `GroupData::GroupData`, `level > 0`).  This walks every group
// the run actually declared and asks `get_group_fluxes` about each one, so the
// answer covers thorns that are active, thorns that are compiled in, tags
// spelled with the wrong case and tags this file has never heard of.  `[P186]`.
//
// It is ONE SHOT and it prints the EXAMINED count as well as the hit count,
// because `[P135]`/`[P184]`: a zero is only a zero if the instrument spoke.
//
// NOTE that `get_group_fluxes` asserts on a malformed `fluxes` tag.  On a run
// with a refined level the driver calls it anyway and the assert is already
// reachable; on a single-level run this instrument would be the first caller,
// which is one more reason it is off unless asked for.
void log_sched_flux_census(const cGH *const cctkGH) {
  static bool done = false;
  if (done)
    return;
  done = true;

  const int ngroups = CCTK_NumGroups();
  int gf_groups = 0, with_tag = 0;
  std::ostringstream tagged;
  bool first = true;
  for (int gi = 0; gi < ngroups; ++gi) {
    if (CCTK_GroupTypeI(gi) != CCTK_GF)
      continue;
    ++gf_groups;
    const std::array<int, dim> fluxes = get_group_fluxes(gi);
    if (fluxes[0] < 0)
      continue;
    ++with_tag;
    if (!first)
      tagged << ",";
    first = false;
    const char *const name = CCTK_FullGroupName(gi);
    tagged << (name ? name : "?");
  }

  std::ostringstream fields;
  fields << "site=fluxcensus ngroups=" << ngroups << " gf_groups=" << gf_groups
         << " with_fluxes_tag=" << with_tag
         << " tagged=" << (with_tag == 0 ? std::string("-") : tagged.str());
  log_sched(cctkGH, fields.str());
}

} // namespace

void Reflux(const cGH *cctkGH, int level);
void Restrict(const cGH *cctkGH, int level, const std::vector<int> &groups,
              const char *site);
void Restrict(const cGH *cctkGH, int level, const char *site);

namespace {
// Convert a (direction, face) pair to an AMReX Orientation
amrex::Orientation orient(int d, int f) {
  return amrex::Orientation(d, amrex::Orientation::Side(f));
}
int GroupStorageCrease(const cGH *cctkGH, int n_groups, const int *groups,
                       const int *requested_tls, int *status, const bool inc);

#ifdef CCTK_DEBUG
// mp_slave_4.md §6 instrumentation: at every call site that invokes
// MultiPatch_Interpolate, log the current AMR epoch (CarpetX_GetEpoch(),
// the same value CapyrX::MultiPatch1_Interpolate's Step 1 compares against
// its cache's stored epoch to decide whether to rebuild) and whether
// CoordinatesX::vcoordx -- and hence this call's own write-back -- is among
// the variables this particular call will (re)settle. This answers the
// scheduling-order question mp_slave_4.md left open: does the interpolation
// cache's one-time rebuild read CoordinatesX's own ghost-zone coordinates
// before or after some call has corrected them, and if the correction only
// ever happens inside the very call that also rebuilds the cache, the fresh
// (post-write) coordinate can never be the one the (never-rebuilt-again)
// cache's routing decisions were based on.
void log_mp_interpolate_call(const char *site,
                             const std::vector<CCTK_INT> &cactusvarinds) {
  static const bool log_donors = std::getenv("CAPYRX_LOG_DONORS") != nullptr;
  if (!log_donors)
    return;
  static const int vcoordx_varind = CCTK_VarIndex("CoordinatesX::vcoordx");
  const bool has_coordinatesx =
      vcoordx_varind >= 0 &&
      std::find(cactusvarinds.begin(), cactusvarinds.end(), vcoordx_varind) !=
          cactusvarinds.end();
#pragma omp critical
  CCTK_VINFO("MPINTERP_CALL site=\"%s\" epoch=%d nvars=%zu "
             "has_coordinatesx=%d",
             site, int(CarpetX_GetEpoch()), cactusvarinds.size(),
             int(has_coordinatesx));
}
#endif

// AMR-B4a (`[P123]`, `[P127]`).  THE INTER-BOX CONSISTENCY CHECK.
//
// THE INVARIANT.  A global grid point that several AMReX boxes of the same
// patch and level hold is stored several times, and the copies must agree.
// Nothing in this driver checked that, and `[P127]` is the finding that no
// shipped checker could: every one of them judges an output row on its own and
// never compares two rows for the same cell, so `[P123]` -- a slaved interior
// cell whose inter-box ghost copies hold the patch's own PRE-slave value --
// sat inside a green battery from before Phase B until C1 ran `rowdiff.py` by
// hand.
//
// "THE BOX THAT OWNS A NODE" IS NOT A THING, AND ASSUMING IT IS OVERCOUNTS.
// A vertex-centred group is stored on a NODAL `BoxArray` (`driver.cxx:1008`,
// `convert(..., IndexType::NODE)`), and adjacent nodal boxes SHARE the plane
// of nodes between them: both hold it, and both hold it as VALID data, not as
// a ghost.  Measured on `color_ghost.par` at one rank: 3103 of 23491 nodes lie
// in more than one box's valid region.  A first version of this check assumed
// a unique owner, compared each ghost copy once per covering box and folded
// the flags onto "the" owner; it reported 1936 disagreeing cells on that rig
// where `rowdiff.py`, on the same run's own output, reported 1624.  So:
//
//   * a copy is compared ONCE per box that holds it, not once per
//     (holder, coverer) pair -- component 1 of `flag` is the visited marker
//     that enforces it;
//   * a node's valid copies are compared against each other too, because on a
//     shared plane there is no ghost involved and the old formulation was
//     blind to it;
//   * and a node is COUNTED at the lowest-indexed box whose valid region holds
//     it, so `bad_cells` is a count of distinct global nodes.
//
// WHAT IT COUNTS.  Three numbers, because they are three different questions:
//
//   checked      copies (cell x component) that at least one OTHER box also
//                holds in its valid region -- the redundant copies, the only
//                ones this invariant is about.  It is the instrument's own
//                "did it speak" column: at ONE AMReX box per patch nothing is
//                redundant and every count below is vacuously zero (`[P35]`,
//                `[N13]`, `[P177]`).  A gate that reads `bad_values = 0`
//                without also reading `checked > 0` has measured nothing.
//   bad_values   of those, how many differ from the reference copy.  IT IS
//                REFERENCE-DEPENDENT and deliberately not pinned: for a node
//                holding {a, a, b} it is 1 or 2 depending on which copy
//                `ParallelCopy` happened to deliver.  Zero is
//                reference-INdependent, which is the only thing gated.
//   bad_cells    distinct global NODES not all of whose copies agree.  That is
//                exactly the quantity `rowdiff.py` calls a multi-valued cell,
//                so it is the one comparable with C1 gate 5b's pinned
//                400 / 1624 / 5040.
//
// HOW THE REFERENCE COPY IS OBTAINED, AND WHY IT IS NOT `FillBoundary`.
// `FillBoundary` is exactly the operation AMR-B4b adds as the FIX, and a
// detector built out of the fix's own call would be checking its own
// arithmetic.  `ParallelCopy` from the source's VALID region into the
// destination's valid+ghost region states the invariant through AMReX's other
// communication path (a `CPC` plan instead of an `FB` plan).  Cells that no
// box's valid region covers -- the interpatch ghosts `I` and the outer ghosts
// `O`, which lie outside the patch domain -- are not reached by it, keep the
// value `MultiFab::Copy` put there, and are neither checked nor counted.  That
// is deliberate: they have exactly one writer each and this invariant does not
// apply to them.
//
// THE DEDUPLICATION IS AN `ADD` PARALLELCOPY, NOT A GATHER.  A disagreeing
// copy writes a 1 into `flag`; `flag` is then `ParallelCopy`-ed into VALID
// regions with `FabArrayBase::ADD`, which lands every copy of a node on every
// box whose valid region holds it, and the count is taken at the lowest-indexed
// of those.  No index gather and no hash table.
//
// PERIODICITY.  The comparison follows `geom.periodicity()`, so a ghost filled
// across a periodic face is checked.  The DEDUPLICATION does not: a periodic
// image carries a different index, and `bad_cells` counts distinct indices.
// On a periodic geometry `bad_values` can therefore exceed what `bad_cells`
// accounts for, which is why both are gated at zero and neither alone.
//
// COMPARISON IS EXACT, NOT TOLERANT.  Two copies of the same number produced
// by the same write are bit-identical, so any difference at all is the defect
// and a tolerance would only hide the small ones.  Two NaNs are treated as
// agreeing (`NaN == NaN` is false, and a poisoned cell copied to a ghost is
// consistent, not broken); a NaN against a number is a disagreement that
// cannot contribute to `max_absdiff`, so it is counted separately.
struct interbox_report_t {
  long checked = 0;
  long bad_values = 0;
  long bad_cells = 0;
  long nan_mismatch = 0;
  long boxes = 0;
  double max_absdiff = 0.0;

  void operator+=(const interbox_report_t &o) {
    checked += o.checked;
    bad_values += o.bad_values;
    bad_cells += o.bad_cells;
    nan_mismatch += o.nan_mismatch;
    boxes += o.boxes;
    using std::max;
    max_absdiff = max(max_absdiff, o.max_absdiff);
  }
};

interbox_report_t interbox_check_one(const amrex::MultiFab &mfab,
                                     const amrex::Geometry &geom) {
  interbox_report_t rep;
  const amrex::IntVect ng = mfab.nGrowVect();
  if (ng.max() <= 0)
    return rep;

  const int ncomp = mfab.nComp();
  const amrex::BoxArray &ba = mfab.boxArray();
  const amrex::DistributionMapping &dm = mfab.DistributionMap();
  const amrex::Periodicity period = geom.periodicity();

  // The reference copy, delivered into every cell some valid region covers.
  // The `Copy` first, so that cells NO valid region covers compare equal and
  // cannot be mistaken for a disagreement; the `ParallelCopy` then overwrites
  // exactly the covered ones.
  amrex::MultiFab ref(ba, dm, ncomp, ng, amrex::MFInfo(), mfab.Factory());
  amrex::MultiFab::Copy(ref, mfab, 0, 0, ncomp, ng);
  ref.ParallelCopy(mfab, 0, 0, ncomp, amrex::IntVect(0), ng, period);

  // component 0: this copy disagrees with the reference.
  // component 1: this copy has already been visited.  A cell can be covered by
  // several other boxes at once (nodal shared planes, and ordinary corners),
  // and must be compared once per HOLDER, not once per (holder, coverer) pair.
  amrex::MultiFab flag(ba, dm, 2, ng, amrex::MFInfo(), mfab.Factory());
  flag.setVal(0.0, 0, 2, ng);

  std::vector<std::pair<int, amrex::Box> > isects;
  for (amrex::MFIter mfi(mfab); mfi.isValid(); ++mfi) {
    const int me = mfi.index();
    const amrex::Box &vbx = mfi.validbox();
    const amrex::Box gbx = amrex::grow(vbx, ng);
    const amrex::Array4<const CCTK_REAL> have = mfab.const_array(mfi);
    const amrex::Array4<const CCTK_REAL> want = ref.const_array(mfi);
    const amrex::Array4<CCTK_REAL> flg = flag.array(mfi);
    ++rep.boxes;

    // `shiftIntVect()` is `{(0,0,0)}` on a non-periodic geometry, so this loop
    // degenerates to one pass there.  It is written out because the ghost fill
    // this check is a postcondition for uses `geom.periodicity()` too, and a
    // check that ignored a periodic direction would report a clean patch that
    // the fill had left stale.
    for (const amrex::IntVect &shift : period.shiftIntVect()) {
      amrex::Box sbx(gbx);
      sbx.shift(-shift);
      ba.intersections(sbx, isects);
      for (const auto &is : isects) {
        if (shift == amrex::IntVect(0) && is.first == me)
          continue; // a box does not make its own copy redundant
        amrex::Box ibx(is.second);
        ibx.shift(shift);
        ibx &= gbx;
        const auto lo = amrex::lbound(ibx);
        const auto hi = amrex::ubound(ibx);
        for (int k = lo.z; k <= hi.z; ++k) {
          for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
              if (flg(i, j, k, 1) > 0.5)
                continue; // already compared, through another coverer
              flg(i, j, k, 1) = 1.0;
              bool any = false;
              for (int c = 0; c < ncomp; ++c) {
                ++rep.checked;
                const CCTK_REAL a = have(i, j, k, c);
                const CCTK_REAL b = want(i, j, k, c);
                if (a == b)
                  continue;
                using std::isnan;
                const bool na = isnan(a), nb = isnan(b);
                if (na && nb)
                  continue;
                ++rep.bad_values;
                any = true;
                if (na || nb) {
                  ++rep.nan_mismatch;
                } else {
                  using std::fabs;
                  using std::max;
                  rep.max_absdiff =
                      max(rep.max_absdiff, double(fabs(double(a) - double(b))));
                }
              }
              if (any)
                flg(i, j, k, 0) = 1.0;
            }
          }
        }
      }
    }
  }

  // Fold every disagreeing copy of a node onto every box whose valid region
  // holds that node, then count it at the lowest-indexed of them.  The `ADD`
  // is what carries a ghost copy's verdict back to the node; the lowest-index
  // rule is what makes the count one per NODE rather than one per valid box
  // that happens to share the plane it sits on.
  amrex::MultiFab own(ba, dm, 1, 0, amrex::MFInfo(), mfab.Factory());
  own.setVal(0.0, 0, 1, 0);
  own.ParallelCopy(flag, 0, 0, 1, ng, amrex::IntVect(0),
                   amrex::Periodicity::NonPeriodic(),
                   amrex::FabArrayBase::ADD);
  std::vector<std::pair<int, amrex::Box> > isects2;
  std::vector<amrex::Box> claimed;
  for (amrex::MFIter mfi(own); mfi.isValid(); ++mfi) {
    const int me = mfi.index();
    const amrex::Box &vbx = mfi.validbox();
    const amrex::Array4<const CCTK_REAL> m = own.const_array(mfi);
    claimed.clear();
    ba.intersections(vbx, isects2);
    for (const auto &is : isects2)
      if (is.first < me)
        claimed.push_back(is.second);
    const auto lo = amrex::lbound(vbx);
    const auto hi = amrex::ubound(vbx);
    for (int k = lo.z; k <= hi.z; ++k) {
      for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
          if (m(i, j, k, 0) <= 0.5)
            continue;
          const amrex::IntVect iv(i, j, k);
          bool taken = false;
          for (const amrex::Box &b : claimed)
            if (b.contains(iv)) {
              taken = true;
              break;
            }
          if (!taken)
            ++rep.bad_cells;
        }
      }
    }
  }

  return rep;
}

// One line per (site, group), on the I/O process, AFTER a global reduction --
// which is why `CCTK_VINFO` is safe here and was not in A5's census: the
// numbers printed are already the whole grid's (`[P223]`, `[P242]`).
void interbox_check(const char *const site, const std::vector<int> &groups) {
  DECLARE_CCTK_PARAMETERS;
  if (!check_interbox_consistency)
    return;
  if (!active_levels)
    return;

  static Timer timer("CarpetX::interbox_check");
  Interval interval(timer);

  for (const int gi : groups) {
    interbox_report_t tot;
    active_levels->loop_serially([&](auto &restrict leveldata) {
      const auto &restrict groupdata = *leveldata.groupdata.at(gi);
      const int ntls = groupdata.mfab.size();
      const int sync_tl = ntls > 1 ? ntls - 1 : ntls;
      const amrex::Geometry &geom =
          ghext->patchdata.at(leveldata.patch).amrcore->Geom(leveldata.level);
      for (int tl = 0; tl < sync_tl; ++tl)
        tot += interbox_check_one(*groupdata.mfab.at(tl), geom);
    });

    amrex::ParallelDescriptor::ReduceLongSum(tot.checked);
    amrex::ParallelDescriptor::ReduceLongSum(tot.bad_values);
    amrex::ParallelDescriptor::ReduceLongSum(tot.bad_cells);
    amrex::ParallelDescriptor::ReduceLongSum(tot.nan_mismatch);
    amrex::ParallelDescriptor::ReduceLongSum(tot.boxes);
    amrex::ParallelDescriptor::ReduceRealMax(tot.max_absdiff);

    if (amrex::ParallelDescriptor::IOProcessor())
      CCTK_VINFO("INTERBOX site=\"%s\" group=%s boxes=%ld checked=%ld "
                 "bad_values=%ld bad_cells=%ld nan_mismatch=%ld "
                 "max_absdiff=%.17g%s",
                 site,
                 ghext->patchdata.at(0)
                     .leveldata.at(0)
                     .groupdata.at(gi)
                     ->groupname.c_str(),
                 tot.boxes, tot.checked,
                 tot.bad_values, tot.bad_cells, tot.nan_mismatch,
                 tot.max_absdiff,
                 tot.checked == 0 ? " VACUOUS(nothing-is-redundant)" : "");
  }
}

// AMR-B4b (`[P123]`).  RE-SHARE THE INTERPOLATOR'S INTERIOR WRITES.
//
// THE DEFECT.  Within one sync the order is: AMReX fills the inter-box ghosts
// (`tasks1/2/3`), BC pass 1, `MultiPatch_Interpolate`, BC pass 2.  With
// `CapyrX_MultiPatch::slave_overlap = yes` the interpolator does not only fill
// interpatch ghosts, it OVERWRITES interior cells -- the overlap-band cells a
// patch holds but does not own.  Nothing re-fills the inter-box ghosts after
// that write, so for the rest of the sync every ghost copy of a slaved cell
// holds the patch's own pre-slave value and any stencil in the neighbouring
// box reads it.  Measured cell by cell in `[P123]`: 400 / 1624 / 5040
// disagreeing cells on `color` / `color_ghost` / `color_ghost_overlap`, zero
// with slaving off, at ONE rank as well as two -- it is a multi-BOX defect,
// not a multi-rank one -- and `[P157]` found the same signature on the BBH
// (10.8 % of shared nodes, 97 % of them in the overlap band).
//
// THE FIX IS THE FILL THAT IS MISSING, WITH THE ARGUMENTS THE SYNC'S OWN FILL
// USED.  `FillBoundary` writes exactly the ghost cells another box's valid
// region covers within this patch and level.  Interpatch ghosts and outer
// ghosts lie outside the patch domain, are covered by no box, and are
// therefore untouched -- so the pass the interpolator has just made is not
// clobbered, and neither is BC pass 1's outer-ghost write.  The strongest form
// of that argument is not about box coverage at all: this is the same call
// with the same arguments as `FillPatch_Sync`'s own `FillBoundary_nowait`, so
// it writes the same cell set, no more and no less.
//
// WHY BEFORE THE CORNERS-ONLY SECOND BC PASS AND NOT AFTER.  The corner cells
// that pass writes are computed from cells inside the box; running the
// re-share first means they are computed from settled data.  It also has to be
// before the validity marks, which follow that pass.
//
// WHY IT IS GATED AND NOT UNCONDITIONAL.  With no interior write between the
// two fills the second one is a bit-for-bit re-copy, so it is inert -- but it
// is a collective, and paying it on every synced group of every multipatch run
// that does not use `slave_overlap` is a real cost for nothing.  The gate is
// the multipatch thorn's own answer, so a `slave_overlap = no` run is
// bit-identical to what it was, which is also what keeps `[P142]`'s column
// alive.
//
// IF THE QUERY IS NOT ALIASED, RE-SHARE ANYWAY.  A multipatch thorn that does
// not answer has not said "no".  The failure mode of re-sharing when it was
// unnecessary is one redundant copy; the failure mode of not re-sharing when
// it was necessary is `[P123]`, silently.
bool multipatch_interpolate_writes_interior() {
  static const bool answer = []() -> bool {
    if (!CCTK_IsFunctionAliased("MultiPatch_Interpolate"))
      return false;
    if (!CCTK_IsFunctionAliased("MultiPatch_InterpolateWritesInterior")) {
      CCTK_VWARN(CCTK_WARN_ALERT,
                 "The multipatch thorn provides MultiPatch_Interpolate but not "
                 "MultiPatch_InterpolateWritesInterior, so this driver cannot "
                 "tell whether the interpatch fill also writes interior cells. "
                 "Assuming it does, and re-sharing the inter-box ghost zones "
                 "after every interpatch fill. This is correct but not free; "
                 "provide the function to switch it off.");
      return true;
    }
    return MultiPatch_InterpolateWritesInterior() != 0;
  }();
  return answer;
}

void reshare_interior_writes(const std::vector<int> &groups) {
  if (!multipatch_interpolate_writes_interior())
    return;
  assert(active_levels);

  static Timer timer("CarpetX::reshare_interior_writes");
  Interval interval(timer);

  // Serial over groups and then over (level, patch), for the reason stated at
  // the sync's own fill loop: every one of these is collective, and the ranks
  // have to enter them in the same order.
  //
  // `tl = 0` only, because that is the time level `MultiPatch_Interpolate`
  // writes. A re-share at `tl >= 1` would copy cells nothing had changed.
  for (const int gi : groups) {
    active_levels->loop_serially([&](auto &restrict leveldata) {
      auto &restrict groupdata = *leveldata.groupdata.at(gi);
      amrex::MultiFab &mfab = *groupdata.mfab.at(0);
      const amrex::Geometry &geom =
          ghext->patchdata.at(leveldata.patch).amrcore->Geom(leveldata.level);
      mfab.FillBoundary(0, mfab.nComp(), mfab.nGrowVect(), geom.periodicity());
    });
  }
}

// BUGFIX_TODO.md step D3 (C10).  The regrid path's interpatch repair, in ONE
// place instead of two byte-identical copies -- one in `Initialise`, one in
// `Evolve`.
//
// WHAT THIS PASS IS FOR.  `FillPatch_NewLevel` / `FillPatch_RemakeLevel` apply
// boundary conditions to a freshly created or remade level BEFORE
// `MultiPatch_Interpolate` has run on it, so a ghost cell that is
// simultaneously in an interpatch direction and on an outer-BC face -- an
// "interpatch corner" -- was written from an interpatch ghost zone nobody had
// filled.  This runs the interpolator and then rewrites exactly those corners,
// mirroring what step B2 did on the sync path.
//
// TWO THINGS IT USED TO DO THAT IT HAD NO BUSINESS DOING, both measured on
// `evidence/fix/d3/pars/d3_cart_edge_L2_Pyes.par`:
//
//  1. IT RAN `bc_pass_t::all`.  That is every face, edge and corner of every
//     box -- a second full outer-boundary write over a level FillPatch had just
//     written -- when the only cells with a stale source are the corners.  It is
//     now `bc_pass_t::interpatch_corners_only`, which is what the comment above
//     it has always claimed and what `bc_pass_t` was added for.  On a
//     single-patch grid no cell is an interpatch corner, so the pass now writes
//     nothing instead of rewriting the whole boundary.
//
//  2. IT RAN OVER EVERY `CCTK_GF` GROUP, including groups whose interior does
//     not hold a value yet.  On a level created by this very regrid,
//     `CapyrX_MultiPatch`'s `vertex_Jacobians` and `vertex_dJacobians` and
//     `CoordinatesX`'s three coordinate groups are unwritten: they are written
//     `(everywhere)` at `CCTK_BASEGRID`, which is traversed FOUR LINES BELOW
//     this call.  Both the interpolation and the boundary pass therefore read
//     them.  With `poison_undefined_values = yes` that read is a NaN --
//     measured, `MULTIPATCH::VERTEX_JACOBIANS` patch 0 level 1
//     `bc_pass = all`, `Assertion !isnan(val) failed` at
//     `boundaries_impl.hxx:738`, backtrace in
//     `evidence/fix/d3/report/gdb_poison_group.txt`.  With poisoning off the
//     same read returns whatever the allocator left there.  Nothing observable
//     came of it, because `CCTK_BASEGRID` overwrites the result immediately --
//     but a grid function whose interior is not valid must not be a SOURCE, and
//     `vertex_dJacobians` has 18 components, which is also how this pass reached
//     the boundary kernel's component-count guard (C10) on a single-patch rig.
//
// THE PREDICATE IS VALIDITY, NOT A LIST OF GROUP NAMES.  `poison_invalid_gf`
// poisons exactly the regions the validity flags call invalid
// (`valid.cxx:186-207`), so "it held poison" and "its interior is not valid" are
// the same statement, and asking the flags asks the driver's own record instead
// of compiling another thorn's schedule into it.  A group is skipped when ANY
// active level, time level or variable reports an invalid interior.  Excluding
// is the safe direction: everything this pass would have done to such a group is
// redone by the `CCTK_BASEGRID` and `CCTK_POSTREGRID` traverses that follow it.
//
// WHAT IS NOT DISJOINT HERE, AND WHY C-AMR IS WHY IT DOES NOT MATTER.
//
// The paragraph this replaces said that the interpatch corners of a freshly
// created level "would still be written twice", and excused the absence of a
// failing-before test by citing a PARAMCHECK refusal of multipatch together
// with `max_num_levels > 1`.  The excuse is gone -- that refusal was deleted
// and replaced by the C-AMR contract in `fillpatch.cxx`.  The claim itself is
// TRUE, and it is measured: on `patch_system = "Thornburg06"`, whose every
// patch carries interpatch angular faces AND a physical radial boundary, a
// level-1 fill of the whole of patch 0 writes 2304 interpatch-corner cells at
// `FillPatch_NewLevel` before `MultiPatch_Interpolate` has run on that level,
// and the corners-only pass below then writes them again.
//
// `FillPatch_NewLevel` and `FillPatch_RemakeLevel` nevertheless still run
// `bc_pass_t::all` over the real `mfab`, and making them
// `skip_interpatch_corners` -- the exact mirror of what this pass does -- is
// deliberately NOT done.  Three reasons, in the order of how much they carry:
//
//  1. C-AMR ALREADY EMPTIES THE SET, AND IT IS ENFORCED AT THOSE VERY SITES.
//     A region is an interpatch corner only if the fine FAB box -- the valid
//     box grown by `nghostzones` -- extends OUTSIDE the patch domain in an
//     interpatch direction.  C-AMR forbids exactly that: it requires the
//     coarse temporary, which is `CoarseBox(coarsen(grow(box, ng) & domain))`,
//     to stay inside the coarse domain across every interpatch face, and
//     since `CoarseBox` grows by the prolongation stencil, a non-negative
//     clearance means the grow was never clipped in the first place.  So
//     under the contract there is no such region to write, on ANY patch
//     system.  Measured on the one patch system where the corner set is
//     non-empty when the contract is broken, one knob apart: contract
//     violated (whole patch refined, `multipatch_amr_contract = "warn"`) ->
//     40 boxes reached, 2304 corner cells; contract held (the refined region
//     moved off the angular faces and left spanning both radial ones,
//     clearance +7) -> 2 boxes reached, 4356 cells written on the physical
//     radial faces, and ZERO interpatch corners.  `has_interpatch = 1` and
//     `has_outerbc = 1` in both columns, so nothing about the geometry
//     changed except the contract.
//
//  2. AND ON THE CUBED SPHERE THE SET IS EMPTY A SECOND WAY, independently of
//     the contract.  `BoxInBox_Setup` returns early for `cctk_patch != 0`
//     (`BoxInBox/src/boxinbox.cxx`), so the only patch that can carry a
//     refined level is patch 0, and patch 0 of a cubed sphere has all six
//     faces `symmetry_t::interpatch` -- which, since step B7 stopped storing
//     the configured outer BC there, means `boundary_t::none` on all six.
//     Measured: `has_outerbc = 0` at both sites on every cubed-sphere rig,
//     including one that violates C-AMR and refines the whole of patch 0.
//     This leg is narrower than reason 1 and is second for that reason.
//
//  3. AND THE CHANGE WOULD NOT BE FREE.  `MakeNewLevelFromCoarse` marks
//     `valid_int | valid_ghosts | outer_valid` and calls `check_valid_gf`
//     with `forbid_nans` ten lines after `FillPatch_NewLevel` returns
//     (`driver.cxx`), while this corners-only pass does not run until after
//     `CCTK_Traverse("CCTK_BASEGRID")`.  That ordering is measured, not read:
//     a poisoned Thornburg06 leg stops at `valid.cxx`, "MakeNewLevelFromCoarse
//     after prolongation", before this function is ever reached.  The sync
//     path is disjoint only because step B2 also DEFERRED its validity marks
//     past its second pass; the regrid path never got that half.  So on a
//     patch that did refine while carrying an outer BC, skipping the corner
//     would leave it poisoned under a validity claim -- and for
//     `boundary_t::dirichlet`, which reads no ghost at all and would
//     otherwise have written the corner correctly, that is an abort where
//     today there is none.
//
// WHAT WOULD HAVE TO CHANGE THIS.  Admitting a level > 0 fill whose fine box
// leaves the patch domain across an interpatch face -- i.e. relaxing C-AMR,
// or the `warn` hatch used as a mode of operation rather than for diagnosis.
// On that day this pass and the validity marking in `driver.cxx` have to move
// together, and moving only this one would be worse than leaving both alone.
// AMR-B3.  THE REPAIR IS IN TWO HALVES, AND THE SPLIT *IS* THE FIX.
//
// THE DEFECT (AMR-D5).  This pass ran BEFORE `CCTK_Traverse("CCTK_BASEGRID")`
// four lines below, and `MultiPatch_Interpolate` interpolates AT vertex
// coordinates that `CoordinatesX_Setup` writes IN that traverse.  On a level
// this regrid just created or remade those coordinates hold poison
// (`poison_undefined_values = yes`) or allocator residue, and CapyrX read them
// as physical coordinates: `[P212]` exit 1 at `cubed_sphere.cxx:177`, `[P246]`
// 9126 interpatch queries at (0,0,0) answered from the cube's centre with no
// message at all, `[P257]` exit 134 at `interpolate.cxx:293` on the production
// geometry in 13 of 24 attempts with poisoning OFF.
//
// WHY NOT SIMPLY MOVE THE CALL AFTER THE TRAVERSE.  Because the variable list
// and the coordinates become available at DIFFERENT points of the regrid step,
// and the list is the half that must be chosen early:
//
//   * `CoordinatesX_Setup` declares `WRITES: cell_coords(everywhere)` and
//     `cell_volume(everywhere)`, and `CallFunction` marks every written
//     variable valid on every active level (`:2363-2387`).  A call placed
//     after the traverse therefore passes the `interior_is_valid` filter below
//     for two `CENTERING={ccc}` groups, and CapyrX's own B6 refusal
//     `CCTK_VERROR`s the moment `npoints > 0`
//     (`CapyrX_MultiPatch/src/interpolate.cxx:880-897`).  EVERY two-level
//     multipatch run would abort at startup.
//   * It would also admit `CoordinatesX::vertex_coords`,
//     `CapyrX::vertex_Jacobians` and `vertex_dJacobians` -- all `{vvv}`, all
//     written `(everywhere)` at BASEGRID -- and replace their analytically
//     exact ghost values with order-4 interpolated ones.
//
// The list as chosen HERE is the prolongated, checkpointed groups: the ones
// that hold a value on a level this regrid made, and the only ones that need
// this repair at all.  The groups BASEGRID writes are written EVERYWHERE,
// ghost zones included, and never needed it.  So the two halves:
// `_select` runs where the old call did and answers "which variables?";
// `_apply` runs after the traverse and answers "with which coordinates?".
//
// WHAT IS STILL NOT TESTED, said here rather than found upstream.  The
// `Evolve` twin at `:1931` is not exercised by any rig in this tree:
// `[P211]`/`[P270]` measured `RemakeLevel` called ZERO times across A3's, A6's
// and B1's whole matrices, because the AMR rigs that move boxes declare no
// checkpointed group and so never reach a fill site.  Half of this change
// ships on a code reading.
std::vector<CCTK_INT> regrid_interpatch_repair_select(const char *const site) {
  assert(active_levels);
  const int ngroups = CCTK_NumGroups();

  // Which groups hold a value everywhere this pass would read one?
  std::vector<char> interior_is_valid(ngroups, 1);
  active_levels->loop_serially([&](auto &restrict leveldata) {
    for (int gi = 0; gi < ngroups; ++gi) {
      if (!interior_is_valid.at(gi))
        continue;
      if (CCTK_GroupTypeI(gi) != CCTK_GF)
        continue;
      const auto &restrict groupdata = *leveldata.groupdata.at(gi);
      const int ntls = groupdata.mfab.size();
      const int sync_tl = ntls > 1 ? ntls - 1 : ntls;
      for (int tl = 0; tl < sync_tl; ++tl)
        for (int vi = 0; vi < groupdata.numvars; ++vi)
          if (!groupdata.valid.at(tl).at(vi).get().valid_int)
            interior_is_valid.at(gi) = 0;
    }
  });

  std::vector<CCTK_INT> cactusvarinds;
  for (int gi = 0; gi < ngroups; ++gi) {
    if (CCTK_GroupTypeI(gi) != CCTK_GF)
      continue;
    if (!interior_is_valid.at(gi))
      continue;
    const auto &groupdata =
        *ghext->patchdata.at(0).leveldata.at(0).groupdata.at(gi);
    for (int var = 0; var < groupdata.numvars; ++var)
      cactusvarinds.push_back(groupdata.firstvarindex + var);
  }

#ifdef CCTK_DEBUG
  log_mp_interpolate_call(site, cactusvarinds);
#else
  (void)site;
#endif

  return cactusvarinds;
}

// AMR-B3, the second half: run what `regrid_interpatch_repair_select` chose,
// after `CCTK_Traverse("CCTK_BASEGRID")` has written the coordinates
// `MultiPatch_Interpolate` interpolates at.  `active_levels` is still the
// regrid's own range -- the caller sets it once and clears it after
// `CCTK_POSTREGRID` -- so this half covers exactly the levels `_select` asked
// about.
void regrid_interpatch_repair_apply(
    cGH *const cctkGH, const std::vector<CCTK_INT> &cactusvarinds) {
  assert(active_levels);

  if (cactusvarinds.empty())
    return;

  const int ngroups = CCTK_NumGroups();

  // The group set is RECOVERED FROM THE VARIABLE LIST rather than recomputed.
  // Recomputing `interior_is_valid` here would re-ask the validity question at
  // the new position and get a DIFFERENT, larger answer -- which is the whole
  // hazard this split exists to avoid. `_select` pushes every variable of
  // every admitted group, so this reconstruction is exact rather than
  // conservative.
  std::vector<char> in_list(ngroups, 0);
  for (const CCTK_INT varind : cactusvarinds) {
    const int gi = CCTK_GroupIndexFromVarI(varind);
    assert(gi >= 0 && gi < ngroups);
    in_list.at(gi) = 1;
  }

  // Standalone call: no later pass in this regrid step reads its output.
  //
  // Two of this comment's clauses had gone stale and are replaced rather than
  // patched.  "B8 refuses the configuration that reaches this block at all"
  // is no longer true: that PARAMCHECK was deleted and replaced by the C-AMR
  // contract at the top of `fillpatch.cxx` and the C-AMR2 pre-pass in
  // `CapyrX_MultiPatch`, which ADMIT the configuration where the contracts
  // hold and refuse it by name where they do not.  And this call no longer
  // covers every level: `active_levels` is the regrid's own
  // `(first_modified, last_modified + 1)` range, and on a contract-holding
  // geometry it writes nothing at all.  What survives unchanged is A8's
  // measurement that it is dead at `max_num_levels = 1` ([P28]).
  if (log_sched_on()) {
    std::ostringstream fields;
    fields << "site=regrid/mpinterp-pre nvars=" << cactusvarinds.size();
    log_sched(cctkGH, fields.str());
  }
  MultiPatch_Interpolate(cctkGH, cactusvarinds.size(), cactusvarinds.data());
  if (log_sched_on()) {
    std::ostringstream fields;
    fields << "site=regrid/mpinterp-post nvars=" << cactusvarinds.size();
    log_sched(cctkGH, fields.str());
  }

  // AMR-B4b, the second site.  `[P276]` measured this call writing no slaved
  // cell on today's geometry (`nslaved = 0`, `[P217]`), so this is predicted
  // inert here -- and it ships anyway, because that zero is a property of a
  // geometry and of correct coordinates, not of this code: a refinement box in
  // the overlap band would give level 1 slaved cells and this site would then
  // leave them stale exactly as the sync did.  AMR-B4a's own check runs at the
  // end of this function and measures whether it was inert.
  {
    std::vector<int> repaired;
    for (int gi = 0; gi < ngroups; ++gi)
      if (in_list.at(gi))
        repaired.push_back(gi);
    reshare_interior_writes(repaired);
  }

  active_levels->loop_serially([&](auto &restrict leveldata) {
    for (int gi = 0; gi < ngroups; ++gi) {
      if (CCTK_GroupTypeI(gi) != CCTK_GF)
        continue;
      if (!in_list.at(gi))
        continue;
      auto &restrict groupdata = *leveldata.groupdata.at(gi);
      const int ntls = groupdata.mfab.size();
      const int sync_tl = ntls > 1 ? ntls - 1 : ntls;
      for (int tl = 0; tl < sync_tl; ++tl)
        groupdata.apply_boundary_conditions(*groupdata.mfab.at(tl),
                                            bc_pass_t::interpatch_corners_only);
    }
  });

  // AMR-B4a: the same postcondition the sync path checks, at the OTHER site
  // that runs `MultiPatch_Interpolate`.  `[P276]` measured this call writing
  // no slaved cell at all on today's geometry (`nslaved = 0`, because
  // `[P217]`'s level-1 candidate set is empty), so this is expected to read
  // zero -- and a zero that is measured is worth more than a zero that is
  // inferred from a level count.
  {
    std::vector<int> repaired;
    for (int gi = 0; gi < ngroups; ++gi)
      if (in_list.at(gi))
        repaired.push_back(gi);
    interbox_check("regrid_interpatch_repair", repaired);
  }
}
} // namespace

////////////////////////////////////////////////////////////////////////////////

GridDesc::GridDesc(const GHExt::PatchData::LevelData &leveldata,
                   const MFPointer &mfp) {
  DECLARE_CCTK_PARAMETERS;

  // The number of ghostzones in each direction
  // for (int d = 0; d < dim; ++d)
  //   nghostzones[d] = mfp.nGrowVect()[d];
  nghostzones = {ghost_size >= 0 ? ghost_size : ghost_size_x,
                 ghost_size >= 0 ? ghost_size : ghost_size_y,
                 ghost_size >= 0 ? ghost_size : ghost_size_z};

  const auto &patchdata = ghext->patchdata.at(leveldata.patch);
  const amrex::IntVect ng(nghostzones[0], nghostzones[1], nghostzones[2]);
  const amrex::Box &domain = patchdata.amrcore->Geom(leveldata.level).Domain();
  const amrex::Box &vbx = mfp.validbox(); // interior region (without ghosts)
  const amrex::Box &fbx = mfp.fabbox(ng); // allocated array
  const amrex::Box &gbx = mfp.growntilebox(ng); // current region (with ghosts)

  for (int d = 0; d < dim; ++d)
    assert(domain.type(d) == amrex::IndexType::CELL);

  // Level, patch, and component
  level = leveldata.level;
  patch = leveldata.patch;
  component = mfp.index();

  // Global shape
  for (int d = 0; d < dim; ++d)
    gsh[d] = domain[orient(d, 1)] + 1 - domain[orient(d, 0)] + 1 +
             2 * nghostzones[d];

  // Local shape
  for (int d = 0; d < dim; ++d)
    lsh[d] = fbx[orient(d, 1)] - fbx[orient(d, 0)] + 1 + 1;

  // Allocated shape
  for (int d = 0; d < dim; ++d)
    ash[d] = lsh[d];

  // Local extent
  for (int d = 0; d < dim; ++d) {
    lbnd[d] = fbx[orient(d, 0)] + nghostzones[d];
    ubnd[d] = fbx[orient(d, 1)] + 1 + nghostzones[d];
  }

  // Boundaries
  for (int d = 0; d < dim; ++d)
    for (int f = 0; f < 2; ++f)
      bbox[f][d] = vbx[orient(d, f)] == domain[orient(d, f)];

  // Thread tile box
  for (int d = 0; d < dim; ++d) {
    tmin[d] = gbx[orient(d, 0)] - fbx[orient(d, 0)];
    // For vertex centred grids, the allocated box is 1 vertex larger
    // than the number of cells, and AMReX assigns this extra vertex
    // to the final tile
    assert(gbx[orient(d, 1)] <= fbx[orient(d, 1)]);
    tmax[d] = gbx[orient(d, 1)] + 1 - fbx[orient(d, 0)] +
              (gbx[orient(d, 1)] == fbx[orient(d, 1)]);
  }

  const amrex::Geometry &geom = patchdata.amrcore->Geom(0);
  const CCTK_REAL *restrict const global_x0 = geom.ProbLo();
  const CCTK_REAL *restrict const global_dx = geom.CellSize();
  for (int d = 0; d < dim; ++d) {
    const int levfac = 1 << leveldata.level;
    // Offset between this level's and the coarsest level's origin as
    // multiple of the grid spacing
    const int levoff = (1 - levfac) * (1 - 2 * nghostzones[d]);
    const int levoffdenom = 2;
    // Vertex-centred coordinates on coarse level
    const CCTK_REAL origin_space =
        global_x0[d] + (1 - 2 * nghostzones[d]) * global_dx[d] / 2;
    const CCTK_REAL delta_space = global_dx[d];
    // Vertex-centred coordinates on current level
    dx[d] = delta_space / levfac;
    x0[d] = origin_space + dx[d] * levoff / levoffdenom;
  }

  // Check constraints
  for (int d = 0; d < dim; ++d) {
    // Domain size
    assert(gsh[d] >= 0);

    // Local size
    assert(lbnd[d] >= 0);
    assert(lsh[d] >= 0);
    assert(lbnd[d] + lsh[d] <= gsh[d]);
    assert(ubnd[d] == lbnd[d] + lsh[d] - 1);

    // Internal representation
    assert(ash[d] >= 0);
    assert(ash[d] >= lsh[d]);

    // Ghost zones
    assert(nghostzones[d] >= 0);
    assert(2 * nghostzones[d] <= lsh[d]);

    // Tiles
    assert(tmin[d] >= 0);
    assert(tmin[d] <= tmax[d]);
    assert(tmax[d] <= lsh[d]);
  }
}

GridDesc::GridDesc(const GHExt::PatchData::LevelData &leveldata,
                   const int component) {
  // `global_component` is the global component index.
  // There is no tiling.

  const auto &patchdata = ghext->patchdata.at(leveldata.patch);

  const amrex::FabArrayBase &fab = *leveldata.fab;

  const amrex::Box &fbx = fab.fabbox(component); // allocated array
  const amrex::Box &vbx =
      fab.box(component);      // interior region (without ghosts)
  const amrex::Box &gbx = fbx; // current region (with ghosts)
  const amrex::Box &domain = patchdata.amrcore->Geom(leveldata.level).Domain();

  for (int d = 0; d < dim; ++d)
    assert(domain.type(d) == amrex::IndexType::CELL);

  // Level, patch, and component
  level = leveldata.level;
  patch = leveldata.patch;
  this->component = component;

  // The number of ghostzones in each direction
  for (int d = 0; d < dim; ++d)
    nghostzones[d] = fab.nGrowVect()[d];

  // Global shape
  for (int d = 0; d < dim; ++d)
    gsh[d] = domain[orient(d, 1)] + 1 - domain[orient(d, 0)] + 1 +
             2 * nghostzones[d];

  // Local shape
  for (int d = 0; d < dim; ++d)
    lsh[d] = fbx[orient(d, 1)] - fbx[orient(d, 0)] + 1 + 1;

  // Allocated shape
  for (int d = 0; d < dim; ++d)
    ash[d] = lsh[d];

  // Local extent
  for (int d = 0; d < dim; ++d) {
    lbnd[d] = fbx[orient(d, 0)] + nghostzones[d];
    ubnd[d] = fbx[orient(d, 1)] + 1 + nghostzones[d];
  }

  // Boundaries
  const auto &symmetries = ghext->patchdata.at(leveldata.patch).symmetries;
  for (int d = 0; d < dim; ++d)
    for (int f = 0; f < 2; ++f)
      bbox[f][d] = vbx[orient(d, f)] == domain[orient(d, f)] &&
                   symmetries[f][d] != symmetry_t::none;

  // Thread tile box
  for (int d = 0; d < dim; ++d) {
    tmin[d] = gbx[orient(d, 0)] - fbx[orient(d, 0)];
    // For vertex centred grids, the allocated box is 1 vertex larger
    // than the number of cells, and AMReX assigns this extra vertex
    // to the final tile
    assert(gbx[orient(d, 1)] <= fbx[orient(d, 1)]);
    tmax[d] = gbx[orient(d, 1)] + 1 - fbx[orient(d, 0)] +
              (gbx[orient(d, 1)] == fbx[orient(d, 1)]);
  }

  const amrex::Geometry &geom = patchdata.amrcore->Geom(0);
  const CCTK_REAL *restrict const global_x0 = geom.ProbLo();
  const CCTK_REAL *restrict const global_dx = geom.CellSize();
  for (int d = 0; d < dim; ++d) {
    const int levfac = 1 << leveldata.level;
    // Offset between this level's and the coarsest level's origin as
    // multiple of the grid spacing
    const int levoff = (1 - levfac) * (1 - 2 * nghostzones[d]);
    const int levoffdenom = 2;
    // Vertex-centred coordinates on coarse level
    const CCTK_REAL origin_space =
        global_x0[d] + (1 - 2 * nghostzones[d]) * global_dx[d] / 2;
    const CCTK_REAL delta_space = global_dx[d];
    // Vertex-centred coordinates on current level
    dx[d] = delta_space / levfac;
    x0[d] = origin_space + dx[d] * levoff / levoffdenom;
  }

  // Check constraints
  for (int d = 0; d < dim; ++d) {
    // Domain size
    assert(gsh[d] >= 0);

    // Local size
    assert(lbnd[d] >= 0);
    assert(lsh[d] >= 0);
    assert(lbnd[d] + lsh[d] <= gsh[d]);
    assert(ubnd[d] == lbnd[d] + lsh[d] - 1);

    // Internal representation
    assert(ash[d] >= 0);
    assert(ash[d] >= lsh[d]);

    // Ghost zones
    assert(nghostzones[d] >= 0);
    assert(2 * nghostzones[d] <= lsh[d]);

    // Tiles
    assert(tmin[d] >= 0);
    assert(tmin[d] <= tmax[d]);
    assert(tmax[d] <= lsh[d]);
  }
}

GridPtrDesc::GridPtrDesc(const GHExt::PatchData::LevelData &leveldata,
                         const MFPointer &mfp)
    : GridDesc(leveldata, mfp) {
  const amrex::IntVect ng(nghostzones[0], nghostzones[1], nghostzones[2]);
  const amrex::Box &fbx = mfp.fabbox(ng); // allocated array
  cactus_offset = lbound(fbx);
}

GridPtrDesc1::GridPtrDesc1(
    const GHExt::PatchData::LevelData &leveldata,
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    const MFPointer &mfp)
    : GridDesc(leveldata, mfp) {
  DECLARE_CCTK_PARAMETERS;
  const amrex::IntVect ng(nghostzones[0], nghostzones[1], nghostzones[2]);
  const amrex::Box &fbx = mfp.fabbox(ng); // allocated array
  cactus_offset = lbound(fbx);
  for (int d = 0; d < dim; ++d) {
    assert(groupdata.nghostzones.at(d) >= 0);
    assert(groupdata.nghostzones.at(d) <= nghostzones[d]);
  }
  for (int d = 0; d < dim; ++d)
    gimin[d] = nghostzones[d] - groupdata.nghostzones.at(d);
  for (int d = 0; d < dim; ++d)
    gimax[d] = lsh[d] - groupdata.indextype.at(d) -
               (nghostzones[d] - groupdata.nghostzones.at(d));
  for (int d = 0; d < dim; ++d)
    gash[d] = ash[d] - groupdata.indextype.at(d) -
              2 * (nghostzones[d] - groupdata.nghostzones.at(d));
}

////////////////////////////////////////////////////////////////////////////////

cGH *copy_cctkGH(const cGH *restrict const sourceGH) {
  cGH *restrict const cctkGH = new cGH;

  // Copy all fields by default
  *cctkGH = *sourceGH;

  // Allocate most pointers anew
  const auto copy_array = [](const auto *restrict const srcptr, const int sz) {
    using T = std::decay_t<decltype(*srcptr)>;
    T *restrict const ptr = new T[sz];
    std::copy(srcptr, srcptr + sz, ptr);
    return ptr;
  };
  cctkGH->cctk_gsh = copy_array(sourceGH->cctk_gsh, dim);
  cctkGH->cctk_lsh = copy_array(sourceGH->cctk_lsh, dim);
  cctkGH->cctk_lbnd = copy_array(sourceGH->cctk_lbnd, dim);
  cctkGH->cctk_ubnd = copy_array(sourceGH->cctk_ubnd, dim);
  cctkGH->cctk_tile_min = copy_array(sourceGH->cctk_tile_min, dim);
  cctkGH->cctk_tile_max = copy_array(sourceGH->cctk_tile_max, dim);
  cctkGH->cctk_ash = copy_array(sourceGH->cctk_ash, dim);
  cctkGH->cctk_to = copy_array(sourceGH->cctk_to, dim);
  cctkGH->cctk_from = copy_array(sourceGH->cctk_from, dim);
  cctkGH->cctk_delta_space = copy_array(sourceGH->cctk_delta_space, dim);
  cctkGH->cctk_origin_space = copy_array(sourceGH->cctk_origin_space, dim);
  cctkGH->cctk_bbox = copy_array(sourceGH->cctk_bbox, 2 * dim);
  cctkGH->cctk_levfac = copy_array(sourceGH->cctk_levfac, dim);
  cctkGH->cctk_levoff = copy_array(sourceGH->cctk_levoff, dim);
  cctkGH->cctk_levoffdenom = copy_array(sourceGH->cctk_levoffdenom, dim);
  cctkGH->cctk_nghostzones = copy_array(sourceGH->cctk_nghostzones, dim);

  const int numvars = CCTK_NumVars();
  cctkGH->data = new void **[numvars];
  for (int vi = 0; vi < numvars; ++vi)
    cctkGH->data[vi] =
        copy_array(sourceGH->data[vi], CCTK_DeclaredTimeLevelsVI(vi));

  return cctkGH;
}

void delete_cctkGH(cGH *cctkGH) {
  delete[] cctkGH->cctk_gsh;
  delete[] cctkGH->cctk_lsh;
  delete[] cctkGH->cctk_lbnd;
  delete[] cctkGH->cctk_ubnd;
  delete[] cctkGH->cctk_tile_min;
  delete[] cctkGH->cctk_tile_max;
  delete[] cctkGH->cctk_ash;
  delete[] cctkGH->cctk_to;
  delete[] cctkGH->cctk_from;
  delete[] cctkGH->cctk_delta_space;
  delete[] cctkGH->cctk_origin_space;
  delete[] cctkGH->cctk_bbox;
  delete[] cctkGH->cctk_levfac;
  delete[] cctkGH->cctk_levoff;
  delete[] cctkGH->cctk_levoffdenom;
  delete[] cctkGH->cctk_nghostzones;
  const int numvars = CCTK_NumVars();
  for (int vi = 0; vi < numvars; ++vi)
    delete[] cctkGH->data[vi];
  delete[] cctkGH->data;
#ifdef CCTK_DEBUG
  memset(cctkGH, 0, sizeof *cctkGH);
#endif
  delete cctkGH;
}

enum class mode_t { unknown, local, patch, level, global, meta };

mode_t current_mode(const cGH *restrict cctkGH) {
  const bool have_local = cctkGH->cctk_component != undefined;
  const bool have_patch = cctkGH->cctk_patch != undefined;
  const bool have_level = cctkGH->cctk_level != undefined;
  const bool have_global = cctkGH->cctk_nghostzones[0] != undefined;
  if (have_local && have_patch && have_level && have_global)
    return mode_t::local;
  else if (!have_local && have_patch && have_level && have_global)
    return mode_t::patch;
  else if (!have_local && !have_patch && have_level && have_global)
    return mode_t::level;
  else if (!have_local && !have_patch && !have_level && have_global)
    return mode_t::global;
  else if (!have_local && !have_patch && !have_level && !have_global)
    return mode_t::meta;
  else
    assert(0);
}

bool in_local_mode(const cGH *restrict cctkGH) {
  return current_mode(cctkGH) == mode_t::local;
}

bool in_patch_mode(const cGH *restrict cctkGH) {
  return current_mode(cctkGH) == mode_t::patch;
}

bool in_level_mode(const cGH *restrict cctkGH) {
  return current_mode(cctkGH) == mode_t::level;
}

bool in_global_mode(const cGH *restrict cctkGH) {
  return current_mode(cctkGH) == mode_t::global;
}

bool in_meta_mode(const cGH *restrict cctkGH) {
  return current_mode(cctkGH) == mode_t::meta;
}

// Initialize cctkGH entries
void setup_cctkGH(cGH *restrict cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  // Dimensions
  cctkGH->cctk_dim = 3;

  // Grid function alignment
  // TODO: Check whether AMReX guarantees a particular alignment
  cctkGH->cctk_alignment = 1;
  cctkGH->cctk_alignment_offset = 0;

  // The refinement factor in time over the top level (coarsest) grid
  cctkGH->cctk_timefac = 1; // no subcycling

  // The total number of patches
  cctkGH->cctk_npatches = ghext->num_patches();

  // The convergence level (numbered from zero upwards)
  cctkGH->cctk_convlevel = 0; // no convergence tests

  // Initialize grid spacing
  for (int d = 0; d < dim; ++d) {
    cctkGH->cctk_origin_space[d] = NAN;
    cctkGH->cctk_delta_space[d] = NAN;
  }

  // Initialize time stepping
  cctkGH->cctk_time = 0;
  cctkGH->cctk_delta_time = NAN;

  // init into meta mode
  cctkGH->cctk_component = undefined;
  cctkGH->cctk_level = undefined;
  cctkGH->cctk_patch = undefined;
  cctkGH->cctk_nghostzones[0] = undefined;
  assert(in_meta_mode(cctkGH));
}

// Update fields that carry state and change over time
void update_cctkGH(cGH *const cctkGH, const cGH *const sourceGH) {
  if (cctkGH == sourceGH)
    return;
  cctkGH->cctk_iteration = sourceGH->cctk_iteration;
  cctkGH->cctk_time = sourceGH->cctk_time;
  cctkGH->cctk_delta_time = sourceGH->cctk_delta_time;
  // for (int d = 0; d < dim; ++d)
  //   cctkGH->cctk_origin_space[d] = sourceGH->cctk_origin_space[d];
  // for (int d = 0; d < dim; ++d)
  //   cctkGH->cctk_delta_space[d] = sourceGH->cctk_delta_space[d];
}

// Set cctkGH entries for global mode
void enter_global_mode(cGH *restrict cctkGH) {
  DECLARE_CCTK_PARAMETERS;
  assert(in_meta_mode(cctkGH));

  // The number of ghostzones in each direction
  // TODO: Get this from mfab (mfab.fb_ghosts)
  cctkGH->cctk_nghostzones[0] = ghost_size >= 0 ? ghost_size : ghost_size_x;
  cctkGH->cctk_nghostzones[1] = ghost_size >= 0 ? ghost_size : ghost_size_y;
  cctkGH->cctk_nghostzones[2] = ghost_size >= 0 ? ghost_size : ghost_size_z;

  // Grid scalar pointers
  {
    auto &restrict globaldata = ghext->globaldata;
    const int num_groups = CCTK_NumGroups();
    for (int gi = 0; gi < num_groups; ++gi) {
      cGroup group;
      int ierr = CCTK_GroupData(gi, &group);
      assert(!ierr);

      if (group.grouptype != CCTK_SCALAR && group.grouptype != CCTK_ARRAY) {
        continue;
      } else { // CCTK_ARRAY or CCTK_SCALAR
        auto &restrict arraygroupdata = *globaldata.arraygroupdata.at(gi);
        for (int tl = 0; tl < int(arraygroupdata.data.size()); ++tl) {
          const auto &restrict vars = arraygroupdata.data.at(tl);
          for (int vi = 0; vi < arraygroupdata.numvars; ++vi)
            cctkGH->data[arraygroupdata.firstvarindex + vi][tl] =
                const_cast<void *>(
                    vars.data_at(vi * arraygroupdata.array_size));
        }
      }
    }
  }

  assert(in_global_mode(cctkGH));
}
void leave_global_mode(cGH *restrict cctkGH) {
  assert(in_global_mode(cctkGH));

  // Grid scalar pointers
  {
    auto &restrict globaldata = ghext->globaldata;
    const int num_groups = CCTK_NumGroups();
    for (int gi = 0; gi < num_groups; ++gi) {
      cGroup group;
      int ierr = CCTK_GroupData(gi, &group);
      assert(!ierr);

      if (group.grouptype != CCTK_SCALAR && group.grouptype != CCTK_ARRAY) {
        continue;
      } else { // CCTK_ARRAY or CCTK_SCALAR
        auto &restrict arraygroupdata = *globaldata.arraygroupdata.at(gi);
        for (int tl = 0; tl < int(arraygroupdata.data.size()); ++tl)
          for (int vi = 0; vi < arraygroupdata.numvars; ++vi)
            cctkGH->data[arraygroupdata.firstvarindex + vi][tl] = nullptr;
      }
    }
  }

  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_nghostzones[d] = undefined;

  assert(in_meta_mode(cctkGH));
}

// Set cctkGH entries for level mode
void enter_level_mode(cGH *restrict cctkGH, const int level) {
  DECLARE_CCTK_PARAMETERS;
  assert(in_global_mode(cctkGH));

  cctkGH->cctk_level = level;
  for (int d = 0; d < dim; ++d) {
    // The refinement factor over the top level (coarsest) grid
    const int levfac = 1 << level;
    cctkGH->cctk_levfac[d] = levfac;
    // Offset between this level's and the coarsest level's origin as multiple
    // of the grid spacing
    const int levoff = (1 - levfac) * (1 - 2 * cctkGH->cctk_nghostzones[d]);
    const int levoffdenom = 2;
    cctkGH->cctk_levoff[d] = levoff;
    cctkGH->cctk_levoffdenom[d] = levoffdenom;
  }

  assert(in_level_mode(cctkGH));
}
void leave_level_mode(cGH *restrict cctkGH, const int level) {
  assert(in_level_mode(cctkGH));
  cctkGH->cctk_level = undefined;
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_levfac[d] = undefined;
  for (int d = 0; d < dim; ++d) {
    cctkGH->cctk_levoff[d] = undefined;
    cctkGH->cctk_levoffdenom[d] = 0;
  }
  assert(in_global_mode(cctkGH));
}

// Set cctkGH entries for patch mode
void enter_patch_mode(cGH *restrict cctkGH,
                      const GHExt::PatchData::LevelData &restrict leveldata) {
  DECLARE_CCTK_PARAMETERS;
  assert(in_level_mode(cctkGH));

  const auto &patchdata = ghext->patchdata.at(leveldata.patch);

  cctkGH->cctk_patch = leveldata.patch;
  const amrex::Box &domain = patchdata.amrcore->Geom(leveldata.level).Domain();
  const amrex::Geometry &geom = patchdata.amrcore->Geom(0);
  const CCTK_REAL *restrict const global_x0 = geom.ProbLo();
  const CCTK_REAL *restrict const global_dx = geom.CellSize();
  for (int d = 0; d < dim; ++d) {
    // Global shape
    assert(cctkGH->cctk_nghostzones[d] != undefined);
    assert(domain.type(d) == amrex::IndexType::CELL);
    cctkGH->cctk_gsh[d] = domain[orient(d, 1)] + 1 - domain[orient(d, 0)] + 1 +
                          2 * cctkGH->cctk_nghostzones[d];
    // Vertex-centred coarse level coordinates
    const CCTK_REAL origin_space =
        global_x0[d] + (1 - 2 * cctkGH->cctk_nghostzones[d]) * global_dx[d] / 2;
    const CCTK_REAL delta_space = global_dx[d];
    cctkGH->cctk_delta_space[d] = delta_space;
    cctkGH->cctk_origin_space[d] = origin_space;
  }

  assert(in_patch_mode(cctkGH));
}
void leave_patch_mode(cGH *restrict cctkGH,
                      const GHExt::PatchData::LevelData &restrict leveldata) {
  assert(in_patch_mode(cctkGH));
  cctkGH->cctk_patch = undefined;
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_gsh[d] = undefined;
  for (int d = 0; d < dim; ++d) {
    cctkGH->cctk_origin_space[d] = NAN;
    cctkGH->cctk_delta_space[d] = NAN;
  }
  assert(in_level_mode(cctkGH));
}

// Set cctkGH entries for local mode
// TODO: Have separate cctkGH for each patch, level, and local box
void enter_local_mode(cGH *restrict cctkGH,
                      const GHExt::PatchData::LevelData &restrict leveldata,
                      const MFPointer &mfp) {
  assert(in_patch_mode(cctkGH));
  const GridPtrDesc grid(leveldata, mfp);

  cctkGH->cctk_component = mfp.index();
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_lsh[d] = grid.lsh[d];
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_ash[d] = grid.ash[d];
  for (int d = 0; d < dim; ++d) {
    cctkGH->cctk_lbnd[d] = grid.lbnd[d];
    cctkGH->cctk_ubnd[d] = grid.ubnd[d];
  }
  for (int d = 0; d < dim; ++d) {
    cctkGH->cctk_tile_min[d] = grid.tmin[d];
    cctkGH->cctk_tile_max[d] = grid.tmax[d];
  }
  for (int d = 0; d < dim; ++d)
    for (int f = 0; f < 2; ++f)
      cctkGH->cctk_bbox[2 * d + f] = grid.bbox[f][d];

  // Grid function pointers
  const int num_groups = CCTK_NumGroups();
  for (int gi = 0; gi < num_groups; ++gi) {
    cGroup group;
    int ierr = CCTK_GroupData(gi, &group);
    assert(!ierr);

    if (group.grouptype != CCTK_GF)
      continue;

    auto &restrict groupdata = *leveldata.groupdata.at(gi);
    const GridPtrDesc1 grid1(leveldata, groupdata, mfp);
    for (int tl = 0; tl < int(groupdata.mfab.size()); ++tl) {
      const amrex::Array4<CCTK_REAL> vars =
          groupdata.mfab.at(tl)->array(mfp.index());
      for (int vi = 0; vi < groupdata.numvars; ++vi)
        cctkGH->data[groupdata.firstvarindex + vi][tl] = grid1.ptr(vars, vi);
    }
  }

  // Check constraints
  for (int d = 0; d < dim; ++d) {
    // Domain size
    assert(cctkGH->cctk_gsh[d] >= 0);

    // Local size
    assert(cctkGH->cctk_lbnd[d] >= 0);
    assert(cctkGH->cctk_lsh[d] >= 0);
    assert(cctkGH->cctk_lbnd[d] + cctkGH->cctk_lsh[d] <= cctkGH->cctk_gsh[d]);
    assert(cctkGH->cctk_ubnd[d] ==
           cctkGH->cctk_lbnd[d] + cctkGH->cctk_lsh[d] - 1);

    // Tile box
    assert(cctkGH->cctk_tile_min[d] >= 0);
    assert(cctkGH->cctk_tile_min[d] <= cctkGH->cctk_tile_max[d]);
    assert(cctkGH->cctk_tile_max[d] <= cctkGH->cctk_lsh[d]);

    // Internal representation
    assert(cctkGH->cctk_ash[d] >= 0);
    assert(cctkGH->cctk_ash[d] >= cctkGH->cctk_lsh[d]);

    // Ghost zones
    assert(cctkGH->cctk_nghostzones[d] >= 0);
    assert(2 * cctkGH->cctk_nghostzones[d] <= cctkGH->cctk_lsh[d]);
  }

  assert(in_local_mode(cctkGH));
}
void leave_local_mode(cGH *restrict cctkGH,
                      const GHExt::PatchData::LevelData &restrict leveldata,
                      const MFPointer &mfp) {
  assert(in_local_mode(cctkGH));
  cctkGH->cctk_component = undefined;
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_lsh[d] = undefined;
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_ash[d] = undefined;
  for (int d = 0; d < dim; ++d) {
    cctkGH->cctk_lbnd[d] = undefined;
    cctkGH->cctk_ubnd[d] = undefined;
  }
  for (int d = 0; d < dim; ++d) {
    cctkGH->cctk_tile_min[d] = undefined;
    cctkGH->cctk_tile_max[d] = undefined;
  }
  for (int d = 0; d < dim; ++d)
    for (int f = 0; f < 2; ++f)
      cctkGH->cctk_bbox[2 * d + f] = undefined;
  const int num_groups = CCTK_NumGroups();
  for (int gi = 0; gi < num_groups; ++gi) {
    cGroup group;
    int ierr = CCTK_GroupData(gi, &group);
    assert(!ierr);

    if (group.grouptype != CCTK_GF)
      continue;

    auto &restrict groupdata = *leveldata.groupdata.at(gi);
    for (int tl = 0; tl < int(groupdata.mfab.size()); ++tl)
      for (int vi = 0; vi < groupdata.numvars; ++vi)
        cctkGH->data[groupdata.firstvarindex + vi][tl] = nullptr;
  }
  assert(in_patch_mode(cctkGH));
}

// Should this be passed in `cGH`?
int CallFunction_count = -1;
extern "C" CCTK_INT CarpetX_GetCallFunctionCount() {
#ifndef AMREX_USE_GPU
  // CPU: use hardware thread index
  return omp_get_thread_num();
#else
  // GPU: count tiles
  assert(CallFunction_count >= 0);
  return CallFunction_count;
#endif
}

active_levels_t::active_levels_t(const int min_level, const int max_level,
                                 const int min_patch, const int max_patch)
    : min_level(min_level), max_level(max_level), min_patch(min_patch),
      max_patch(max_patch) {
  assert(min_level >= 0);
  assert(max_level <= ghext->num_levels());
  assert(min_patch >= 0);
  assert(max_patch <= ghext->num_patches());
}
active_levels_t::active_levels_t(const int min_level, const int max_level)
    : active_levels_t(min_level, max_level, 0, ghext->num_patches()) {}
active_levels_t::active_levels_t() : active_levels_t(0, ghext->num_levels()) {}

void active_levels_t::assert_consistent_iterations() const {
  rat64 good_iteration = -1;
  for (int level = min_level; level < max_level; ++level) {
    for (int patch = min_patch; patch < max_patch; ++patch) {
      const auto &patchdata = ghext->patchdata.at(patch);
      if (level < int(patchdata.leveldata.size())) {
        const auto &leveldata = patchdata.leveldata.at(level);
        const auto &iteration = leveldata.iteration;
        if (good_iteration == -1)
          good_iteration = iteration;
        assert(iteration == good_iteration);
      }
    }
  }
}

// Loop over all active patches of all active levels from coarsest to
// finest
void active_levels_t::loop_coarse_to_fine(
    const std::function<void(GHExt::PatchData::LevelData &level)> &kernel)
    const {
  assert(omp_get_num_threads() == 1);
  assert_consistent_iterations();
  for (int level = min_level; level < max_level; ++level) {
    for (int patch = min_patch; patch < max_patch; ++patch) {
      auto &patchdata = ghext->patchdata.at(patch);
      if (level < int(patchdata.leveldata.size()))
        kernel(patchdata.leveldata.at(level));
    }
  }
}

// Loop over all active patches of all active levels from finest to
// coarsest
void active_levels_t::loop_fine_to_coarse(
    const std::function<void(GHExt::PatchData::LevelData &level)> &kernel)
    const {
  assert(omp_get_num_threads() == 1);
  assert_consistent_iterations();
  for (int level = max_level - 1; level >= min_level; --level) {
    for (int patch = min_patch; patch < max_patch; ++patch) {
      auto &patchdata = ghext->patchdata.at(patch);
      if (level < int(patchdata.leveldata.size()))
        kernel(patchdata.leveldata.at(level));
    }
  }
}

// Loop over all components of all active patches of all active levels
// in parallel
void active_levels_t::loop_parallel(
    const std::function<void(int patch, int level, int index, int component,
                             const cGH *cctkGH)> &kernel) const {
  assert(omp_get_num_threads() == 1);
  task_manager tasks;

  loop_coarse_to_fine([&](const auto &restrict leveldata) {
    int component = 0;
    const auto mfitinfo = amrex::MFItInfo().DisableDeviceSync().EnableTiling();
    const auto &fab0 = *leveldata.fab;
    for (amrex::MFIter mfi(fab0, mfitinfo); mfi.isValid(); ++mfi, ++component) {
      const int index = mfi.index();
      tasks.submit_serially([&kernel, &leveldata, index, component]() {
        const int patch = leveldata.patch;
        const int level = leveldata.level;
        cGH *restrict const localGH = leveldata.get_local_cctkGH(component);
        kernel(patch, level, index, component, localGH);
      });
    }
  });

  // Run all tasks
  tasks.run_tasks();
  // There is an implicit OpenMP barrier here.
}

// Loop over all components of all active patches of all active levels
// serially
void active_levels_t::loop_serially(
    const std::function<void(int patch, int level, int index, int component,
                             const cGH *cctkGH)> &kernel) const {
  loop_coarse_to_fine([&](const auto &restrict leveldata) {
    int component = 0;
    const auto mfitinfo = amrex::MFItInfo().DisableDeviceSync().EnableTiling();
    const auto &fab0 = *leveldata.fab;
    for (amrex::MFIter mfi(fab0, mfitinfo); mfi.isValid(); ++mfi, ++component) {
      const int index = mfi.index();
      const int patch = leveldata.patch;
      const int level = leveldata.level;
      cGH *restrict const localGH = leveldata.get_local_cctkGH(component);
      kernel(patch, level, index, component, localGH);
    }
  });
}

void synchronize() {
#ifdef AMREX_USE_GPU
  // TODO: Synchronize only if GPU kernels were actually launched
  amrex::Gpu::streamSynchronizeAll();
  AMREX_GPU_ERROR_CHECK();
#endif
}

void update_cctkGHs(cGH *restrict const cctkGH) {
  update_cctkGH(ghext->global_cctkGH.get(), cctkGH);
  for (auto &restrict level_cctkGH : ghext->level_cctkGHs) {
    update_cctkGH(level_cctkGH.get(), cctkGH);
  }
  for (auto &patch : ghext->patchdata) {
    for (auto &restrict level : patch.leveldata) {
      update_cctkGH(level.patch_cctkGH.get(), cctkGH);
    }
  }
  for (auto &patch : ghext->patchdata) {
    for (auto &restrict level : patch.leveldata) {
      for (auto &restrict local_cctkGH : level.local_cctkGHs) {
        update_cctkGH(local_cctkGH.get(), cctkGH);
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void CarpetX_GetLoopBoxAll(const void *restrict const cctkGH_,
                                      const CCTK_INT size,
                                      CCTK_INT *restrict const loop_min,
                                      CCTK_INT *restrict const loop_max) {
  const cGH *restrict const cctkGH = static_cast<const cGH *>(cctkGH_);
  assert(size == dim);
  const GridDescBase grid(cctkGH);
  vect<int, dim> imin, imax;
  // Assumd cell centred box
  grid.box_all<1, 1, 1>(grid.nghostzones, imin, imax);
  for (int d = 0; d < dim; ++d)
    loop_min[d] = imin[d];
  for (int d = 0; d < dim; ++d)
    loop_max[d] = imax[d];
}

extern "C" void CarpetX_GetLoopBoxInt(const void *restrict const cctkGH_,
                                      const CCTK_INT size,
                                      CCTK_INT *restrict const loop_min,
                                      CCTK_INT *restrict const loop_max) {
  const cGH *restrict const cctkGH = static_cast<const cGH *>(cctkGH_);
  assert(size == dim);
  const GridDescBase grid(cctkGH);
  vect<int, dim> imin, imax;
  // Assumd cell centred box
  grid.box_int<1, 1, 1>(grid.nghostzones, imin, imax);
  for (int d = 0; d < dim; ++d)
    loop_min[d] = imin[d];
  for (int d = 0; d < dim; ++d)
    loop_max[d] = imax[d];
}

////////////////////////////////////////////////////////////////////////////////

mode_t decode_mode(const cFunctionData *restrict attribute) {
  bool local_mode = attribute->local;
  bool level_mode = attribute->level;
  bool global_mode = attribute->global;
  bool meta_mode = attribute->meta;
  assert(int(local_mode) + int(level_mode) + int(global_mode) +
             int(meta_mode) <=
         1);
  if (attribute->local)
    return mode_t::local;
  if (attribute->level)
    return mode_t::level;
  if (attribute->global)
    return mode_t::global;
  if (attribute->meta)
    return mode_t::meta;
  return mode_t::local; // default
}

enum class rdwr_t { read, write, invalid };
std::ostream &operator<<(std::ostream &os, const rdwr_t rdwr) {
  switch (rdwr) {
  case rdwr_t::read:
    return os << "read";
  case rdwr_t::write:
    return os << "write";
  case rdwr_t::invalid:
    return os << "invalid";
  default:
    assert(0);
  }
}

struct clause_t {
  int gi, vi, tl;
  valid_t valid;

  friend bool operator==(const clause_t &x, const clause_t &y) {
    return std::make_tuple(x.gi, x.vi, x.tl, x.valid) ==
           std::make_tuple(y.gi, y.vi, y.tl, y.valid);
  }
  friend bool operator<(const clause_t &x, const clause_t &y) {
    return std::make_tuple(x.gi, x.vi, x.tl, x.valid) <
           std::make_tuple(y.gi, y.vi, y.tl, y.valid);
  }

  friend std::ostream &operator<<(std::ostream &os, const clause_t &cl) {
    return os << "clause_t{gi:" << cl.gi << ",vi:" << cl.vi << ",tl:" << cl.tl
              << ",valid:" << cl.valid << "}";
  }
};

std::vector<clause_t> decode_clauses(const cFunctionData *restrict attribute,
                                     const rdwr_t rdwr) {
  std::vector<clause_t> result;
  result.reserve(attribute->n_RDWR);
  for (int n = 0; n < attribute->n_RDWR; ++n) {
    const RDWR_entry &restrict RDWR = attribute->RDWR[n];
    int gi = CCTK_GroupIndexFromVarI(RDWR.varindex);
    assert(gi >= 0);
    int vi = RDWR.varindex - CCTK_FirstVarIndexI(gi);
    assert(vi >= 0 && vi < CCTK_NumVarsInGroupI(gi));
    int tl = RDWR.timelevel;
    int where;
    switch (rdwr) {
    case rdwr_t::read:
      where = RDWR.where_rd;
      break;
    case rdwr_t::write:
      where = RDWR.where_wr;
      break;
    case rdwr_t::invalid:
      where = RDWR.where_inv;
      break;
    default:
      assert(0);
    }
    // Ignore clauses that have no effect
    if (where == 0)
      continue;
    valid_t valid;
    valid.valid_int = where & CCTK_VALID_INTERIOR;
    valid.valid_outer = where & CCTK_VALID_BOUNDARY;
    valid.valid_ghosts = where & CCTK_VALID_GHOSTS;
    result.push_back({gi, vi, tl, valid});
  }
  return result;
}

// Schedule initialisation
int Initialise(tFleshConfig *config) {
  DECLARE_CCTK_PARAMETERS;

  static Timer timer("Initialise");
  Interval interval(timer);

  cGH *restrict const cctkGH = CCTK_SetupGH(config, 0);
  CCTKi_AddGH(config, 0, cctkGH);

  // Check presync mode
  if (!CCTK_EQUALS(presync_mode, "mixed-error") &&
      !CCTK_EQUALS(presync_mode, "presync-only"))
    CCTK_ERROR("CarpetX currently requires Cactus::presync_mode = "
               "\"mixed-error\" or \"presync-only\"");

  // Initialise iteration and time
  cctkGH->cctk_iteration = 0;
  cctkGH->cctk_time = *static_cast<const CCTK_REAL *>(
      CCTK_ParameterGet("cctk_initial_time", "Cactus", nullptr));

  // Initialise schedule
  CCTKi_ScheduleGHInit(cctkGH);

  // Initialise all grid extensions
  CCTKi_InitGHExtensions(cctkGH);

  // Set up cctkGH
  setup_cctkGH(cctkGH);
  enter_global_mode(cctkGH);
  ghext->global_cctkGH = GHExt::cctkGHptr(copy_cctkGH(cctkGH));

  for (const auto &patchdata : ghext->patchdata)
    assert(patchdata.leveldata.empty());
  assert(!active_levels);
  active_levels = std::make_optional<active_levels_t>(0, 0);

  CCTK_Traverse(cctkGH, "CCTK_WRAGH");
  CCTK_Traverse(cctkGH, "CCTK_PARAMCHECK");
  CCTKi_FinaliseParamWarn();

  active_levels = std::optional<active_levels_t>();

  // Set the initial value of max_grid_size for all levels
  CCTK_VINFO("Setting initial values for max_grid_size values for all levels");
  for (const auto &patchdata : ghext->patchdata) {
    amrex::Vector<amrex::IntVect> max_grid_sizes_vec(20);
    for (int level = 0; level < 20; ++level) {
      int size_x = max_grid_sizes_x[level];
      if (size_x == -1)
        size_x = max_grid_size_x;
      int size_y = max_grid_sizes_y[level];
      if (size_y == -1)
        size_y = max_grid_size_y;
      int size_z = max_grid_sizes_z[level];
      if (size_z == -1)
        size_z = max_grid_size_z;
      const amrex::IntVect max_grid_size_vec{size_x, size_y, size_z};
      max_grid_sizes_vec.at(level) = max_grid_size_vec;
    }
    patchdata.amrcore->SetMaxGridSize(max_grid_sizes_vec);
  }

  if (config->recovered) {
    // Recover
#pragma omp critical
    CCTK_VINFO("Recovering from checkpoint...");

    RecoverGridStructure(cctkGH);

    assert(!active_levels);
    active_levels = std::make_optional<active_levels_t>();

    CCTK_Traverse(cctkGH, "CCTK_BASEGRID");

    const char *recovery_mode = *static_cast<const char *const *>(
        CCTK_ParameterGet("recovery_mode", "Cactus", nullptr));
    if (!CCTK_Equals(recovery_mode, "strict")) {
      // Set up initial conditions
      CCTK_Traverse(cctkGH, "CCTK_INITIAL");
      CCTK_Traverse(cctkGH, "CCTK_POSTINITIAL");
      CCTK_Traverse(cctkGH, "CCTK_POSTPOSTINITIAL");
    }

    // Recover
    RecoverGH(cctkGH);
    CCTK_Traverse(cctkGH, "CCTK_RECOVER_VARIABLES");
    CCTK_Traverse(cctkGH, "CCTK_POST_RECOVER_VARIABLES");

    active_levels = std::optional<active_levels_t>();

    // Enable regridding
    for (auto &patchdata : ghext->patchdata)
      patchdata.amrcore->cactus_is_initialized = true;

    // Determine time step size
    if (CCTK_EQUALS(timestep_choice, "timestep")) {
      cctkGH->cctk_delta_time = timestep;
    } else if (CCTK_EQUALS(timestep_choice, "dtfac")) {
      CCTK_REAL mindx = 1.0 / 0.0;
      for (const auto &patchdata : ghext->patchdata) {
        const amrex::Geometry &geom = patchdata.amrcore->Geom(0);
        const CCTK_REAL *restrict const dx = geom.CellSize();
        CCTK_REAL mindx1 = 1.0 / 0.0;
        for (int d = 0; d < dim; ++d)
          mindx1 = fmin(mindx1, dx[d]);
        mindx1 = ldexp(mindx1, -(int(patchdata.leveldata.size()) - 1));
        mindx = fmin(mindx, mindx1);
      }
      cctkGH->cctk_delta_time = dtfac * mindx;
    } else {
      abort();
    }
    using std::isfinite;
    assert(isfinite(cctkGH->cctk_delta_time));
#pragma omp critical
    CCTK_VINFO("Iteration: %d   time: %g   delta_time: %g",
               cctkGH->cctk_iteration, double(cctkGH->cctk_time),
               double(cctkGH->cctk_delta_time));

  } else {
    // Set up initial conditions
#pragma omp critical
    CCTK_VINFO("Setting up initial conditions...");

    // Create coarse grid
    {
      static Timer timer("InitialiseRegrid [coarse]");
      Interval interval(timer);

      const CCTK_REAL time = 0; // dummy time
      for (const auto &patchdata : ghext->patchdata)
        patchdata.amrcore->MakeNewGrids(time);

      // Determine time step size
      if (CCTK_EQUALS(timestep_choice, "timestep")) {
        cctkGH->cctk_delta_time = timestep;
      } else if (CCTK_EQUALS(timestep_choice, "dtfac")) {
        CCTK_REAL mindx = 1.0 / 0.0;
        for (const auto &patchdata : ghext->patchdata) {
          const amrex::Geometry &geom = patchdata.amrcore->Geom(0);
          const CCTK_REAL *restrict const dx = geom.CellSize();
          CCTK_REAL mindx1 = 1.0 / 0.0;
          for (int d = 0; d < dim; ++d)
            mindx1 = fmin(mindx1, dx[d]);
          mindx1 = ldexp(mindx1, -(int(patchdata.leveldata.size()) - 1));
          mindx = fmin(mindx, mindx1);
        }
        cctkGH->cctk_delta_time = dtfac * mindx;
      } else {
        abort();
      }
      using std::isfinite;
      assert(isfinite(cctkGH->cctk_delta_time));
#pragma omp critical
      CCTK_VINFO("Iteration: %d   time: %g   delta_time: %g",
                 cctkGH->cctk_iteration, double(cctkGH->cctk_time),
                 double(cctkGH->cctk_delta_time));

      assert(!active_levels);
      active_levels = std::make_optional<active_levels_t>(0, 1);
      CCTK_Traverse(cctkGH, "CCTK_BASEGRID");
      // CCTK_Traverse(cctkGH, "CCTK_POSTREGRID");
      active_levels = std::optional<active_levels_t>();
    }

    // Output domain information
    if (CCTK_MyProc(nullptr) == 0) {
      for (const auto &patchdata : ghext->patchdata) {
        const int level = 0;
        const cGH *const patchGH =
            ghext->get_patch_cctkGH(level, patchdata.patch);
        const int *restrict const gsh = patchGH->cctk_gsh;
        const int *restrict const nghostzones = patchGH->cctk_nghostzones;
        CCTK_REAL x0[dim], x1[dim], dx[dim];
        for (int d = 0; d < dim; ++d) {
          dx[d] = patchGH->cctk_delta_space[d];
          x0[d] = patchGH->cctk_origin_space[d] +
                  (2 * nghostzones[d] - 1) * dx[d] / 2;
          x1[d] = x0[d] + (gsh[d] - 2 * nghostzones[d] - 1) * dx[d];
        }
#pragma omp critical
        {
          CCTK_VINFO("Patch %d:", patchdata.patch);
          CCTK_VINFO("  Grid extent:");
          CCTK_VINFO("    gsh=[%d,%d,%d]", gsh[0], gsh[1], gsh[2]);
          const auto &bf = patchdata.amrcore->blockingFactor(0);
          CCTK_VINFO("    blocking_factor=[%d,%d,%d]", bf[0], bf[1], bf[2]);
          const auto &mgs = patchdata.amrcore->maxGridSize(0);
          CCTK_VINFO("    max_grid_size=[%d,%d,%d]", mgs[0], mgs[1], mgs[2]);
          const auto mfitinfo = amrex::MFItInfo().EnableTiling();
          const auto &ts = mfitinfo.tilesize;
          CCTK_VINFO("    max_tile_size=[%d,%d,%d]", ts[0], ts[1], ts[2]);
          CCTK_VINFO("  Domain extent:");
          CCTK_VINFO("    xmin=[%.17g,%.17g,%.17g]", double(x0[0]),
                     double(x0[1]), double(x0[2]));
          CCTK_VINFO("    xmax=[%.17g,%.17g,%.17g]", double(x1[0]),
                     double(x1[1]), double(x1[2]));
          CCTK_VINFO("    base dx=[%.17g,%.17g,%.17g]", double(dx[0]),
                     double(dx[1]), double(dx[2]));
          // CCTK_VINFO("  Time stepping:");
          // CCTK_VINFO("    t0=%.17g", double(patchGH->cctk_time));
          // CCTK_VINFO("    dt=%.17g", double(patchGH->cctk_delta_time));
        }
      }
    }

    // Enable regridding. We can only enable regridding after the
    // coarse level has been initialized, since otherwise the error
    // estimate is undefined.
    for (auto &patchdata : ghext->patchdata)
      patchdata.amrcore->cactus_is_initialized = true;

    for (;;) {
      const int level = ghext->num_levels() - 1;
#pragma omp critical
      CCTK_VINFO("Initializing level %d...", level);

      // Check whether a patch has too many levels
      for (const auto &patchdata : ghext->patchdata)
        assert(patchdata.amrcore->finestLevel() <= level);

      assert(!active_levels);
      active_levels = std::make_optional<active_levels_t>(0, level + 1);

      InputGH(cctkGH);
      CCTK_Traverse(cctkGH, "CCTK_INITIAL");
      CCTK_Traverse(cctkGH, "CCTK_POSTINITIAL");
      CCTK_Traverse(cctkGH, "CCTK_POSTPOSTINITIAL");

      active_levels = std::optional<active_levels_t>();

      // Regrid
      bool did_modify_any_level;
      {
#pragma omp critical
        CCTK_VINFO("Regridding...");
        static Timer timer("InitialiseRegrid [refined]");
        Interval interval(timer);

        for (const auto &patchdata : ghext->patchdata) {

          const int old_numlevels = patchdata.amrcore->finestLevel() + 1;
          patchdata.amrcore->level_modified.clear();
          patchdata.amrcore->level_modified.resize(old_numlevels, false);
          const CCTK_REAL time = 0; // dummy time
          patchdata.amrcore->regrid(0, time);

          const int new_numlevels = patchdata.amrcore->finestLevel() + 1;
          const int max_numlevels = patchdata.amrcore->maxLevel() + 1;
          assert(new_numlevels >= 0 && new_numlevels <= max_numlevels);
          assert(new_numlevels == old_numlevels ||
                 new_numlevels == old_numlevels + 1);

#pragma omp critical
          {
            const double pts0 =
                patchdata.leveldata.at(0).fab->boxArray().d_numPts();
            for (const auto &leveldata : patchdata.leveldata) {
              const int sz = leveldata.fab->size();
              const double pts = leveldata.fab->boxArray().d_numPts();
              if (leveldata.level == 0) {
                CCTK_VINFO(
                    "  level %d: %d boxes, %.0f cells (%.4g%%)",
                    leveldata.level, sz, pts,
                    100 * pts /
                        (ldexp(CCTK_REAL(1), dim * leveldata.level) * pts0));
              } else {
                const double ptsc = patchdata.leveldata.at(leveldata.level - 1)
                                        .fab->boxArray()
                                        .d_numPts();
                CCTK_VINFO(
                    "  level %d: %d boxes, %.0f cells (%.4g%%, %.0f%%)",
                    leveldata.level, sz, pts,
                    100 * pts /
                        (ldexp(CCTK_REAL(1), dim * leveldata.level) * pts0),
                    100 * pts / (ldexp(CCTK_REAL(1), dim) * ptsc));
              }
            }
          } // omp critical
        } // for patchdata

        int first_modified_level = INT_MAX;
        int last_modified_level = -1;
        for (const auto &patchdata : ghext->patchdata) {
          for (int lev = 0; lev < int(patchdata.amrcore->level_modified.size());
               ++lev) {
            if (patchdata.amrcore->level_modified.at(lev)) {
              using std::max, std::min;
              first_modified_level = min(first_modified_level, lev);
              last_modified_level = max(last_modified_level, lev);
            }
          }
        }
        did_modify_any_level = last_modified_level >= first_modified_level;

        if (did_modify_any_level) {
          // Determine time step size
          if (CCTK_EQUALS(timestep_choice, "timestep")) {
            cctkGH->cctk_delta_time = timestep;
          } else if (CCTK_EQUALS(timestep_choice, "dtfac")) {
            CCTK_REAL mindx = 1.0 / 0.0;
            for (const auto &patchdata : ghext->patchdata) {
              const amrex::Geometry &geom = patchdata.amrcore->Geom(0);
              const CCTK_REAL *restrict const dx = geom.CellSize();
              CCTK_REAL mindx1 = 1.0 / 0.0;
              for (int d = 0; d < dim; ++d)
                mindx1 = fmin(mindx1, dx[d]);
              mindx1 = ldexp(mindx1, -(int(patchdata.leveldata.size()) - 1));
              mindx = fmin(mindx, mindx1);
            }
            cctkGH->cctk_delta_time = dtfac * mindx;
          } else {
            abort();
          }
          using std::isfinite;
          assert(isfinite(cctkGH->cctk_delta_time));
#pragma omp critical
          CCTK_VINFO("Iteration: %d   time: %g   delta_time: %g",
                     cctkGH->cctk_iteration, double(cctkGH->cctk_time),
                     double(cctkGH->cctk_delta_time));

          assert(!active_levels);
          active_levels = std::make_optional<active_levels_t>(
              first_modified_level, last_modified_level + 1);

          // Regrid path: fill interpatch ghosts + 2nd BC pass for
          // newly created/remade levels.
          //
          // FillPatch_NewLevel / FillPatch_RemakeLevel each call
          // apply_boundary_conditions (1st BC pass) before
          // MultiPatch_Interpolate has run. Corner ghost cells at
          // outer+interpatch face intersections are therefore left with
          // stale values sourced from not-yet-filled interpatch ghosts.
          // Correct them below, mirroring the fix in SyncGroupsByDirI -- but
          // in TWO halves, and see `regrid_interpatch_repair_select` for why.
          // The variable list has to be chosen here, before CCTK_BASEGRID
          // makes six more groups valid; the interpolation has to run after
          // it, because CCTK_BASEGRID is what writes the coordinates it
          // interpolates at (AMR-B3 / AMR-D5).
          static const bool have_multipatch_boundaries =
              CCTK_IsFunctionAliased("MultiPatch_Interpolate");
          std::vector<CCTK_INT> repair_varinds;
          if (have_multipatch_boundaries)
            repair_varinds =
                regrid_interpatch_repair_select("regrid (new/remade levels)");

          CCTK_Traverse(cctkGH, "CCTK_BASEGRID");

          if (have_multipatch_boundaries)
            regrid_interpatch_repair_apply(cctkGH, repair_varinds);

          CCTK_Traverse(cctkGH, "CCTK_POSTREGRID");
          active_levels = std::optional<active_levels_t>();
        }
      } // Regrid

      if (!did_modify_any_level)
        break;
    } // for level
  }
#pragma omp critical
  CCTK_VINFO("Initialized %d levels", ghext->num_levels());

  assert(!active_levels);
  active_levels = std::make_optional<active_levels_t>();

  if (!restrict_during_sync) {
    // Restrict
    assert(active_levels);
    active_levels->loop_fine_to_coarse([&](const auto &leveldata) {
      if (leveldata.level != ghext->num_levels() - 1)
        Restrict(cctkGH, leveldata.level, "init");
    });
    CCTK_Traverse(cctkGH, "CCTK_POSTRESTRICT");
  }

  // Checkpoint, analysis, output
  CCTK_Traverse(cctkGH, "CCTK_POSTSTEP");
  CCTK_Traverse(cctkGH, "CCTK_CPINITIAL");
  CCTK_Traverse(cctkGH, "CCTK_ANALYSIS");
  CCTK_OutputGH(cctkGH);

  active_levels = std::optional<active_levels_t>();

  return 0;
} // namespace CarpetX

bool EvolutionIsDone(cGH *restrict const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  const bool max_iteration_reached = cctkGH->cctk_iteration >= cctk_itlast;

  const bool max_simulation_time_reached =
      cctk_initial_time < cctk_final_time
          ? cctkGH->cctk_time >= cctk_final_time
          : cctkGH->cctk_time <= cctk_final_time;

  const int runtime = CCTK_RunTime();
  const bool max_runtime_reached = runtime >= 60 * max_runtime;

  // Note: Some MPI implementations have been built without support
  // for `MPI_CXX_BOOL`, so we use `int` instead
  int we_are_done;
  if (terminate_next || CCTK_TerminationReached(cctkGH))
    we_are_done = true;
  else if (CCTK_Equals(terminate, "never"))
    we_are_done = false;
  else if (CCTK_Equals(terminate, "iteration"))
    we_are_done = max_iteration_reached;
  else if (CCTK_Equals(terminate, "time"))
    we_are_done = max_simulation_time_reached;
  else if (CCTK_Equals(terminate, "runtime"))
    we_are_done = max_runtime_reached;
  else if (CCTK_Equals(terminate, "any"))
    we_are_done = max_iteration_reached || max_simulation_time_reached ||
                  max_runtime_reached;
  else if (CCTK_Equals(terminate, "all"))
    we_are_done = max_iteration_reached && max_simulation_time_reached &&
                  max_runtime_reached;
  else if (CCTK_Equals(terminate, "either"))
    we_are_done = max_iteration_reached || max_simulation_time_reached;
  else if (CCTK_Equals(terminate, "both"))
    we_are_done = max_iteration_reached && max_simulation_time_reached;
  else
    CCTK_ERROR("internal error");

  // Ensure all processes make the same decision
  MPI_Allreduce(MPI_IN_PLACE, &we_are_done, 1, MPI_INT, MPI_LOR,
                MPI_COMM_WORLD);

  return we_are_done;
}

void InvalidateTimelevels(cGH *restrict const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  static Timer timer("InvalidateTimelevels");
  Interval interval(timer);

  // TODO: Parallelize over groups
  const int num_groups = CCTK_NumGroups();
  for (int gi = 0; gi < num_groups; ++gi) {
    cGroup group;
    int ierr = CCTK_GroupData(gi, &group);
    assert(!ierr);

    if (group.grouptype == CCTK_GF) {
      const auto &patchdata0 = ghext->patchdata.at(0);
      const auto &leveldata0 = patchdata0.leveldata.at(0);
      const auto &groupdata0 = *leveldata0.groupdata.at(gi);
      if (!groupdata0.do_checkpoint) {
        const int ntls0 = groupdata0.mfab.size();
        assert(active_levels);
        active_levels->loop_serially([&](const auto &restrict leveldata) {
          auto &restrict groupdata = *leveldata.groupdata.at(gi);
          // Invalidate all time levels
          const int ntls = groupdata.mfab.size();
          assert(ntls == ntls0);
          for (int tl = 0; tl < ntls; ++tl)
            for (int vi = 0; vi < groupdata.numvars; ++vi)
              groupdata.valid.at(tl).at(vi).set_all(valid_t(), []() {
                return "InvalidateTimelevels (invalidate all non-checkpointed "
                       "variables)";
              });
        });
        // TODO: Parallelize over timelevels and variables
        for (int tl = 0; tl < ntls0; ++tl)
          for (int vi = 0; vi < groupdata0.numvars; ++vi)
            poison_invalid_gf(*active_levels, gi, vi, tl);
      }
    } else { // CCTK_ARRAY or CCTK_SCALAR

      auto &restrict globaldata = ghext->globaldata;
      auto &restrict arraygroupdata = *globaldata.arraygroupdata.at(gi);
      if (!arraygroupdata.do_checkpoint) {
        // Invalidate all time levels
        const int ntls = arraygroupdata.data.size();
        for (int tl = 0; tl < ntls; ++tl)
          for (int vi = 0; vi < arraygroupdata.numvars; ++vi)
            // TODO: handle this more nicely
            arraygroupdata.valid.at(tl).at(vi).set_int(false, []() {
              return "InvalidateTimelevels (invalidate all non-checkpointed "
                     "variables)";
            });
        // TODO: Parallelize over timelevels and variables
        for (int tl = 0; tl < ntls; ++tl)
          for (int vi = 0; vi < arraygroupdata.numvars; ++vi)
            poison_invalid_ga(gi, vi, tl);
      }
    }

  } // for gi
}

void CycleTimelevels(cGH *restrict const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  static Timer timer("CycleTimelevels");
  Interval interval(timer);

  cctkGH->cctk_iteration += 1;
  cctkGH->cctk_time += cctkGH->cctk_delta_time;
  update_cctkGHs(cctkGH);

  // TODO: Parallelize over groups
  const int num_groups = CCTK_NumGroups();
  const bool presync_only = CCTK_EQUALS(presync_mode, "presync-only");
  std::vector<int> presync_groups;
  for (int gi = 0; gi < num_groups; ++gi) {
    cGroup group;
    int ierr = CCTK_GroupData(gi, &group);
    assert(!ierr);

    if (group.grouptype == CCTK_GF) {
      const auto &patchdata0 = ghext->patchdata.at(0);
      const auto &leveldata0 = patchdata0.leveldata.at(0);
      const auto &groupdata0 = *leveldata0.groupdata.at(gi);
      const int ntls0 = groupdata0.mfab.size();
      const nan_handling_t nan_handling = groupdata0.do_checkpoint
                                              ? nan_handling_t::forbid_nans
                                              : nan_handling_t::allow_nans;

      assert(active_levels);
      active_levels->loop_serially([&](auto &restrict leveldata) {
        auto &restrict groupdata = *leveldata.groupdata.at(gi);
        const int ntls = groupdata.mfab.size();
        assert(ntls == ntls0);
        // Rotate time levels and invalidate current time level
        if (ntls > 1) {
          rotate(groupdata.mfab.begin(), groupdata.mfab.end() - 1,
                 groupdata.mfab.end());
          rotate(groupdata.valid.begin(), groupdata.valid.end() - 1,
                 groupdata.valid.end());
          for (int vi = 0; vi < groupdata.numvars; ++vi)
            groupdata.valid.at(0).at(vi).set_all(valid_t(), []() {
              return "CycletimeLevels (invalidate current time level)";
            });
        }
        // All time levels (except the current) must be valid everywhere for
        // checkpointed groups
        if (groupdata.do_checkpoint) {
          for (int tl = (ntls == 1 ? 0 : 1); tl < ntls; ++tl) {
            // it is only possible to sync time-level zero
            if (tl == 0 && presync_only) {
              presync_groups.push_back(gi);
            } else {
              for (int vi = 0; vi < groupdata.numvars; ++vi) {
                error_if_invalid(groupdata, vi, tl, make_valid_all(), []() {
                  return "CycleTimelevels for the state vector";
                });
              }
            }
          }
        }
      });
      for (int vi = 0; vi < groupdata0.numvars; ++vi) {
        if (ntls0 > 1)
          poison_invalid_gf(*active_levels, gi, vi, 0);
        for (int tl = 0; tl < ntls0; ++tl)
          check_valid_gf(*active_levels, gi, vi, tl, nan_handling,
                         []() { return "CycleTimelevels"; });
      }
    } else { // CCTK_ARRAY or CCTK_SCALAR

      auto &restrict globaldata = ghext->globaldata;
      auto &restrict arraygroupdata = *globaldata.arraygroupdata.at(gi);
      const nan_handling_t nan_handling = arraygroupdata.do_checkpoint
                                              ? nan_handling_t::forbid_nans
                                              : nan_handling_t::allow_nans;
      const int ntls = arraygroupdata.data.size();
      // Rotate time levels and invalidate current time level
      if (ntls > 1) {
        rotate(arraygroupdata.data.begin(), arraygroupdata.data.end() - 1,
               arraygroupdata.data.end());
        rotate(arraygroupdata.valid.begin(), arraygroupdata.valid.end() - 1,
               arraygroupdata.valid.end());
        for (int vi = 0; vi < arraygroupdata.numvars; ++vi) {
          arraygroupdata.valid.at(0).at(vi).set_int(false, []() {
            return "CycletimeLevels (invalidate current time level)";
          });
          poison_invalid_ga(gi, vi, 0);
        }
      }
      for (int tl = 0; tl < ntls; ++tl)
        for (int vi = 0; vi < arraygroupdata.numvars; ++vi)
          check_valid_ga(gi, vi, tl, nan_handling,
                         []() { return "CycleTimelevels"; });
    }

  } // for gi
  if (!presync_groups.empty()) {
    SyncGroupsByDirI(cctkGH, presync_groups.size(), presync_groups.data(),
                     nullptr);
  }
}

// Schedule evolution
int Evolve(tFleshConfig *config) {
  DECLARE_CCTK_PARAMETERS;

  static Timer timer("Evolve");
  Interval interval(timer);

  assert(config);
  cGH *restrict const cctkGH = config->GH[0];
  assert(cctkGH);

#pragma omp critical
  CCTK_VINFO("Starting evolution...");

  double total_evolution_time = 0;
  double total_evolution_output_time = 0;
  int total_iterations = 0;
  double total_cell_updates = 0;
  int average_iteration_time_iterations = 0;
  double average_iteration_time = 0;

  std::ofstream performance_file;
  if (out_performance && CCTK_MyProc(NULL) == 0) {
    const int every =
        out_performance_every == -1 ? out_every : out_performance_every;
    if (every > 0) {
      std::ostringstream buf;
      buf << out_dir << "/performance.yaml";
      const std::string filename = buf.str();
      performance_file.open(filename);
      performance_file << "performance:\n" << std::flush;
    }
  }

  while (!EvolutionIsDone(cctkGH)) {

    const double start_time = gettime();

    assert(!active_levels);

    // TODO: Move regridding into a function
    if (regrid_every > 0 && cctkGH->cctk_iteration % regrid_every == 0) {
#pragma omp critical
      CCTK_VINFO("Regridding...");
      static Timer timer("EvolveRegrid");
      Interval interval(timer);

      for (const auto &patchdata : ghext->patchdata) {
        const int old_numlevels = patchdata.amrcore->finestLevel() + 1;
        patchdata.amrcore->level_modified.clear();
        patchdata.amrcore->level_modified.resize(old_numlevels, false);
        const CCTK_REAL time = 0; // dummy time

        // Set the value of max_grid_size before regridding
        CCTK_VINFO(
            "Setting max_grid_size values for all levels before regridding");
        amrex::Vector<amrex::IntVect> max_grid_sizes_vec(20);
        for (int level = 0; level < 20; ++level) {
          int size_x = max_grid_sizes_x[level];
          if (size_x == -1)
            size_x = max_grid_size_x;
          int size_y = max_grid_sizes_y[level];
          if (size_y == -1)
            size_y = max_grid_size_y;
          int size_z = max_grid_sizes_z[level];
          if (size_z == -1)
            size_z = max_grid_size_z;
          const amrex::IntVect max_grid_size_vec{size_x, size_y, size_z};
          max_grid_sizes_vec.at(level) = max_grid_size_vec;
        }
        patchdata.amrcore->SetMaxGridSize(max_grid_sizes_vec);

        patchdata.amrcore->regrid(0, time);

        const int new_numlevels = patchdata.amrcore->finestLevel() + 1;
        const int max_numlevels = patchdata.amrcore->maxLevel() + 1;
        assert(new_numlevels >= 0 && new_numlevels <= max_numlevels);

#pragma omp critical
        {
          CCTK_VINFO("  old levels %d, new levels %d", old_numlevels,
                     new_numlevels);
          double pts0 = patchdata.leveldata.at(0).fab->boxArray().d_numPts();
          assert(!active_levels);
          for (const auto &leveldata : patchdata.leveldata) {
            const int sz = leveldata.fab->size();
            const double pts = leveldata.fab->boxArray().d_numPts();
            if (leveldata.level == 0) {
              CCTK_VINFO(
                  "  level %d: %d boxes, %.0f cells (%.4g%%)", leveldata.level,
                  sz, pts,
                  100 * pts /
                      (ldexp(CCTK_REAL(1), dim * leveldata.level) * pts0));
            } else {
              const double ptsc = patchdata.leveldata.at(leveldata.level - 1)
                                      .fab->boxArray()
                                      .d_numPts();
              CCTK_VINFO(
                  "  level %d: %d boxes, %.0f cells (%.4g%%, %.0f%%)",
                  leveldata.level, sz, pts,
                  100 * pts /
                      (ldexp(CCTK_REAL(1), dim * leveldata.level) * pts0),
                  100 * pts / (ldexp(CCTK_REAL(1), dim) * ptsc));
            }
          }
        } // omp critical
      } // for patchdata

      int first_modified_level = INT_MAX;
      int last_modified_level = -1;
      for (const auto &patchdata : ghext->patchdata) {
        for (int lev = 0; lev < int(patchdata.amrcore->level_modified.size());
             ++lev) {
          if (patchdata.amrcore->level_modified.at(lev)) {
            using std::max, std::min;
            first_modified_level = min(first_modified_level, lev);
            last_modified_level = max(last_modified_level, lev);
          }
        }
      }
      const bool did_modify_any_level =
          last_modified_level >= first_modified_level;

      if (did_modify_any_level) {
        // Determine time step size
        if (CCTK_EQUALS(timestep_choice, "timestep")) {
          cctkGH->cctk_delta_time = timestep;
        } else if (CCTK_EQUALS(timestep_choice, "dtfac")) {
          CCTK_REAL mindx = 1.0 / 0.0;
          for (const auto &patchdata : ghext->patchdata) {
            const amrex::Geometry &geom = patchdata.amrcore->Geom(0);
            const CCTK_REAL *restrict const dx = geom.CellSize();
            CCTK_REAL mindx1 = 1.0 / 0.0;
            for (int d = 0; d < dim; ++d)
              mindx1 = fmin(mindx1, dx[d]);
            mindx1 = ldexp(mindx1, -(int(patchdata.leveldata.size()) - 1));
            mindx = fmin(mindx, mindx1);
          }
          cctkGH->cctk_delta_time = dtfac * mindx;
        } else {
          abort();
        }
        using std::isfinite;
        assert(isfinite(cctkGH->cctk_delta_time));
#pragma omp critical
        CCTK_VINFO("Iteration: %d   time: %g   delta_time: %g",
                   cctkGH->cctk_iteration, double(cctkGH->cctk_time),
                   double(cctkGH->cctk_delta_time));

        assert(!active_levels);
        active_levels = std::make_optional<active_levels_t>(
            first_modified_level, last_modified_level + 1);

        // Regrid path: fill interpatch ghosts + 2nd BC pass for
        // newly created/remade levels.
        //
        // FillPatch_NewLevel / FillPatch_RemakeLevel each call
        // apply_boundary_conditions (1st BC pass) before
        // MultiPatch_Interpolate has run. Corner ghost cells at
        // outer+interpatch face intersections are therefore left with
        // stale values sourced from not-yet-filled interpatch ghosts.
        // Correct them below, mirroring the fix in SyncGroupsByDirI -- but in
        // TWO halves, and see `regrid_interpatch_repair_select` for why. The
        // variable list has to be chosen here, before CCTK_BASEGRID makes six
        // more groups valid; the interpolation has to run after it, because
        // CCTK_BASEGRID is what writes the coordinates it interpolates at
        // (AMR-B3 / AMR-D5).
        //
        // THIS TWIN IS NOT EXERCISED BY ANY RIG IN THIS TREE: `[P211]` and
        // `[P270]` measured `RemakeLevel` called zero times, so the code below
        // ships on a code reading and not on a measurement.
        static const bool have_multipatch_boundaries =
            CCTK_IsFunctionAliased("MultiPatch_Interpolate");
        std::vector<CCTK_INT> repair_varinds;
        if (have_multipatch_boundaries)
          repair_varinds =
              regrid_interpatch_repair_select("regrid (level removal)");

        CCTK_Traverse(cctkGH, "CCTK_BASEGRID");

        if (have_multipatch_boundaries)
          regrid_interpatch_repair_apply(cctkGH, repair_varinds);

        CCTK_Traverse(cctkGH, "CCTK_POSTREGRID");
        active_levels = std::optional<active_levels_t>();
      }
    } // Regrid

    // Find smallest iteration number. Levels at this iteration will
    // be evolved.
    rat64 iteration = ghext->patchdata.at(0).leveldata.at(0).iteration;
    using std::min;
    for (const auto &patchdata : ghext->patchdata)
      for (const auto &leveldata : patchdata.leveldata)
        iteration = min(iteration, leveldata.iteration);

    // Loop over all levels, in batches that combine levels that don't
    // subcycle. The level range is [min_level, max_level).
    int min_level = 0;
    while (min_level < ghext->num_levels()) {
      // Find end of batch
      int max_level = min_level + 1;

      while (max_level < ghext->num_levels()) {
        bool level_is_subcycling_level = false;
        for (const auto &patchdata : ghext->patchdata)
          if (max_level < int(patchdata.leveldata.size()))
            level_is_subcycling_level |=
                patchdata.leveldata.at(max_level).is_subcycling_level;
        if (level_is_subcycling_level)
          break;
        ++max_level;
      }

      // Skip this batch of levels if it is not active at the current
      // iteration
      rat64 level_iteration = -1;
      for (const auto &patchdata : ghext->patchdata)
        if (min_level < int(patchdata.leveldata.size()))
          level_iteration = patchdata.leveldata.at(min_level).iteration;
      assert(level_iteration != -1);
      if (level_iteration > iteration)
        break;

      active_levels = std::make_optional<active_levels_t>(min_level, max_level);

      // Advance iteration number on this batch of levels
      active_levels->loop_serially([&](auto &restrict leveldata) {
        leveldata.iteration += leveldata.delta_iteration;
      });

      // We cannot invalidate all non-evolved variables. ODESolvers
      // calculates things in ODESolvers_Poststep, and we want to use
      // them in the next iteration.
      // InvalidateTimelevels(cctkGH);

      CycleTimelevels(cctkGH);

      CCTK_Traverse(cctkGH, "CCTK_PRESTEP");
      CCTK_Traverse(cctkGH, "CCTK_EVOL");

      // Reflux
      // TODO: These loop bounds are wrong for subcycling
      assert(active_levels);
      for (int level = ghext->num_levels() - 2; level >= 0; --level)
        Reflux(cctkGH, level);

      if (!restrict_during_sync) {
        // Restrict
        //
        // AMR-C7: this is the OTHER schedule position, and on a patch system it
        // is not equivalent to the one inside the sync.  It runs after
        // `CCTK_EVOL` -- i.e. after every interpatch fill of the step -- and it
        // restricts a DIFFERENT set of groups (every `do_checkpoint` group,
        // rather than the groups the sync was asked for).  Both differences are
        // measured in `evidence/amr/c7/`.
        // TODO: These loop bounds are wrong for subcycling
        for (int level = ghext->num_levels() - 2; level >= 0; --level)
          Restrict(cctkGH, level, "evol");
        CCTK_Traverse(cctkGH, "CCTK_POSTRESTRICT");
      }

      CCTK_Traverse(cctkGH, "CCTK_POSTSTEP");
      CCTK_Traverse(cctkGH, "CCTK_CHECKPOINT");
      CCTK_Traverse(cctkGH, "CCTK_ANALYSIS");
      const double output_start_time = gettime();
      CCTK_OutputGH(cctkGH);
      const double output_finish_time = gettime();
      total_evolution_output_time += output_finish_time - output_start_time;

      active_levels = std::optional<active_levels_t>();
    } // for min_level

    const double waiting_start_time = gettime();
    MPI_Barrier(MPI_COMM_WORLD);
    const double waiting_finish_time = gettime();

    const double finish_time = gettime();
    double num_cells = 0;
    for (const auto &patch : ghext->patchdata)
      for (const auto &level : patch.leveldata)
        num_cells += level.fab->boxArray().d_numPts();
    total_cell_updates += num_cells;
    ++total_iterations;
    const double iteration_time = finish_time - start_time;
    total_evolution_time += iteration_time;

    // Calculate statistics
    if (average_iteration_time_iterations < 10)
      ++average_iteration_time_iterations;
    // Calculate exponential moving average
    average_iteration_time =
        ((average_iteration_time_iterations - 1) * average_iteration_time +
         iteration_time) /
        average_iteration_time_iterations;

    const double iterations_per_second = 1 / average_iteration_time;
    const double cell_updates_per_second = num_cells * iterations_per_second;
    CCTK_VINFO("Simulation time: %g   "
               "Iterations per second: %g   "
               "Simulation time per second: %g",
               double(cctkGH->cctk_time), iterations_per_second,
               double(cctkGH->cctk_delta_time * iterations_per_second)

    );
    // This is the same as H-AMR's "cell updates per second":
    CCTK_VINFO("Grid cells: %g   "
               "Grid cell updates per second: %g",
               num_cells, cell_updates_per_second);

    const double total_evolution_compute_time =
        total_evolution_time - total_evolution_output_time;
    CCTK_VINFO("Performance:");
    CCTK_VINFO("  total evolution time:            %g sec",
               total_evolution_time);
    CCTK_VINFO("  total evolution compute time:    %g sec",
               total_evolution_compute_time);
    CCTK_VINFO("  total evolution output time:     %g sec",
               total_evolution_output_time);
    CCTK_VINFO("  total iterations:                %d", total_iterations);
    CCTK_VINFO("  total cells updated:             %g", total_cell_updates);
    CCTK_VINFO("  average iterations per second: %g",
               total_iterations / total_evolution_time);
    CCTK_VINFO("  average cell updates per second: %g",
               total_cell_updates / total_evolution_time);
    // TODO: Output this in a proper I/O method
    if (out_performance && CCTK_MyProc(NULL) == 0) {
      const int every =
          out_performance_every == -1 ? out_every : out_performance_every;
      if (every > 0 && cctkGH->cctk_iteration % every == 0)
        performance_file << "  " << total_iterations << ":\n"
                         << "    evolution-seconds: " << total_evolution_time
                         << "\n"
                         << "    evolution-compute-seconds: "
                         << total_evolution_compute_time << "\n"
                         << "    evolution-output-seconds: "
                         << total_evolution_output_time << "\n"
                         << "    evolution-cell-updates: " << total_cell_updates
                         << "\n"
                         << "    evolution-iterations: " << total_iterations
                         << "\n"
                         << std::flush;
    }

  } // main loop

  if (out_performance && CCTK_MyProc(NULL) == 0)
    performance_file.close();

  return 0;
} // namespace CarpetX

// Schedule shutdown
int Shutdown(tFleshConfig *config) {
  assert(config);
  cGH *restrict const cctkGH = config->GH[0];
  assert(cctkGH);

  static Timer timer("Shutdown");
  Interval interval(timer);

#pragma omp critical
  CCTK_VINFO("Shutting down...");

  assert(!active_levels);
  active_levels = std::make_optional<active_levels_t>();

  CCTK_Traverse(cctkGH, "CCTK_TERMINATE");

  active_levels = std::optional<active_levels_t>();
  active_levels = std::make_optional<active_levels_t>(0, 0);

  CCTK_Traverse(cctkGH, "CCTK_SHUTDOWN");

  active_levels = std::optional<active_levels_t>();
  assert(!ghext);

  return 0;
}

// Call a scheduled function
int CallFunction(void *function, cFunctionData *restrict attribute,
                 void *data) {
  DECLARE_CCTK_PARAMETERS;

  assert(function);
  assert(attribute);
  assert(data);

  cGH *restrict const cctkGH = static_cast<cGH *>(data);

  if (verbose)
#pragma omp critical
    CCTK_VINFO("CallFunction iteration %d %s: %s::%s", cctkGH->cctk_iteration,
               attribute->where, attribute->thorn, attribute->routine);

  static std::map<cFunctionData *restrict, Timer> timers;

  std::map<cFunctionData *restrict, Timer>::iterator timer_iter;
#pragma omp critical(CarpetX_CallFunction)
  {
    timer_iter = timers.find(attribute);
    if (timer_iter == timers.end()) {
      std::ostringstream buf;
      buf << "CallFunction " << attribute->where << ": " << attribute->thorn
          << "::" << attribute->routine;
      timer_iter = std::get<0>(timers.emplace(attribute, buf.str()));
    }
  }
  Timer &timer = timer_iter->second;
  Interval interval(timer);

  assert(active_levels);

  if (CCTK_EQUALS(presync_mode, "presync-only")) {
    const std::vector<clause_t> &reads =
        decode_clauses(attribute, rdwr_t::read);
    std::set<int> sync_set;
    for (const auto &rd : reads) {
      if (CCTK_GroupTypeI(rd.gi) == CCTK_GF) {

        active_levels->loop_serially([&](const auto &restrict leveldata) {
          const auto &restrict groupdata = *leveldata.groupdata.at(rd.gi);
          const valid_t &need = rd.valid;
          valid_t have = groupdata.valid.at(rd.tl).at(rd.vi).get();
          if (need.valid_ghosts && !have.valid_ghosts && have.valid_int)
            sync_set.insert(rd.gi);
        });
      }
    }
    if (!sync_set.empty()) {
      std::vector<int> sync_vec(sync_set.begin(), sync_set.end());
      SyncGroupsByDirI(cctkGH, sync_vec.size(), sync_vec.data(), nullptr);
    }
  }

  // Check whether input variables have valid data
  {
    const std::vector<clause_t> &reads =
        decode_clauses(attribute, rdwr_t::read);
    for (const auto &rd : reads) {
      if (CCTK_GroupTypeI(rd.gi) == CCTK_GF) {
        const auto &patchdata0 = ghext->patchdata.at(0);
        const auto &leveldata0 = patchdata0.leveldata.at(0);
        const auto &groupdata0 = *leveldata0.groupdata.at(rd.gi);
        const nan_handling_t nan_handling = groupdata0.do_checkpoint
                                                ? nan_handling_t::forbid_nans
                                                : nan_handling_t::allow_nans;

        active_levels->loop_serially([&](const auto &restrict leveldata) {
          const auto &restrict groupdata = *leveldata.groupdata.at(rd.gi);
          const valid_t &need = rd.valid;
          error_if_invalid(
              groupdata, rd.vi, rd.tl, need, [attribute, cctkGH]() {
                std::ostringstream buf;
                buf << "CallFunction iteration " << cctkGH->cctk_iteration
                    << " " << attribute->where << ": " << attribute->thorn
                    << "::" << attribute->routine << " checking input";
                return buf.str();
              });
        });
        check_valid_gf(*active_levels, rd.gi, rd.vi, rd.tl, nan_handling,
                       [attribute, cctkGH]() {
                         std::ostringstream buf;
                         buf << "CallFunction iteration "
                             << cctkGH->cctk_iteration << " "
                             << attribute->where << ": " << attribute->thorn
                             << "::" << attribute->routine << " checking input";
                         return buf.str();
                       });
      } else { // CCTK_ARRAY or CCTK_SCALAR

        const auto &restrict arraygroupdata =
            *ghext->globaldata.arraygroupdata.at(rd.gi);
        const nan_handling_t nan_handling = arraygroupdata.do_checkpoint
                                                ? nan_handling_t::forbid_nans
                                                : nan_handling_t::allow_nans;
        const valid_t &need = rd.valid;
        error_if_invalid(
            arraygroupdata, rd.vi, rd.tl, need, [attribute, cctkGH]() {
              std::ostringstream buf;
              buf << "CallFunction iteration " << cctkGH->cctk_iteration << " "
                  << attribute->where << ": " << attribute->thorn
                  << "::" << attribute->routine << " checking input";
              return buf.str();
            });
        check_valid_ga(
            rd.gi, rd.vi, rd.tl, nan_handling, [attribute, cctkGH]() {
              std::ostringstream buf;
              buf << "CallFunction iteration " << cctkGH->cctk_iteration << " "
                  << attribute->where << ": " << attribute->thorn
                  << "::" << attribute->routine << " checking input";
              return buf.str();
            });
      }
    }
  }

  // Poison those output variables that are not input variables
  if (poison_undefined_values) {
    std::map<clause_t, valid_t> isread;
    const std::vector<clause_t> &reads =
        decode_clauses(attribute, rdwr_t::read);
    for (const auto &rd : reads) {
      clause_t cl = rd;
      cl.valid = valid_t();
      assert(isread.count(cl) == 0);
      isread[cl] = rd.valid;
    }
    const std::vector<clause_t> &writes =
        decode_clauses(attribute, rdwr_t::write);
    for (const auto &wr : writes) {
      clause_t cl = wr;
      cl.valid = valid_t();
      valid_t need;
      if (isread.count(cl) > 0)
        need = isread[cl];

      if (CCTK_GroupTypeI(wr.gi) == CCTK_GF) {

        active_levels->loop_serially([&](auto &restrict leveldata) {
          auto &restrict groupdata = *leveldata.groupdata.at(wr.gi);
          const valid_t &provided = wr.valid;
          // The flesh can accidentally describe timelevels that do
          // not exist.
          if (wr.tl > int(groupdata.valid.size()))
            CCTK_VERROR("Accessing non-existent timelevel %d of variable %s",
                        wr.tl, groupdata.groupname.c_str());
          groupdata.valid.at(wr.tl).at(wr.vi).set_invalid(
              provided & ~need,
              [iteration = cctkGH->cctk_iteration, where = attribute->where,
               thorn = attribute->thorn, routine = attribute->routine] {
                std::ostringstream buf;
                buf << "CallFunction iteration " << iteration << " " << where
                    << ": " << thorn << "::" << routine
                    << ": Poison output variables that are not input variables";
                return buf.str();
              });
        });
        poison_invalid_gf(*active_levels, wr.gi, wr.vi, wr.tl);
      } else { // CCTK_ARRAY or CCTK_SCALAR
        auto &restrict arraygroupdata =
            *ghext->globaldata.arraygroupdata.at(wr.gi);
        const valid_t &provided = wr.valid;
        arraygroupdata.valid.at(wr.tl).at(wr.vi).set_invalid(
            provided & ~need,
            [iteration = cctkGH->cctk_iteration, where = attribute->where,
             thorn = attribute->thorn, routine = attribute->routine] {
              std::ostringstream buf;
              buf << "CallFunction iteration " << iteration << " " << where
                  << ": " << thorn << "::" << routine
                  << ": Poison output variables that are not input variables";
              return buf.str();
            });
        poison_invalid_ga(wr.gi, wr.vi, wr.tl);
      }
    }
  }

  // Calculate checksums over variables that are not written
  checksums_t checksums;
  if (poison_undefined_values) {
    const std::vector<clause_t> &writes =
        decode_clauses(attribute, rdwr_t::write);
    const int numgroups = CCTK_NumGroups();
    std::vector<std::vector<std::vector<valid_t> > > gfs(numgroups);
    for (int gi = 0; gi < numgroups; ++gi) {
      const int numvars = CCTK_NumVarsInGroupI(gi);
      gfs.at(gi).resize(numvars);
      for (int vi = 0; vi < numvars; ++vi) {
        const int numtimelevels = 1; // is expanded later if necessary
        gfs.at(gi).at(vi).resize(numtimelevels);
      }
    }
    for (const auto &wr : writes) {
      if (wr.tl >= int(gfs.at(wr.gi).at(wr.vi).size()))
        gfs.at(wr.gi).at(wr.vi).resize(wr.tl + 1);
      gfs.at(wr.gi).at(wr.vi).at(wr.tl) |= wr.valid;
    }

    checksums = calculate_checksums(gfs);
  }

  const mode_t mode = decode_mode(attribute);
  switch (mode) {
  case mode_t::local:
    // Call function once per tile
    active_levels->loop_parallel([&](int patch, int level, int index,
                                     int component, const cGH *local_cctkGH) {
      update_cctkGH(const_cast<cGH *>(local_cctkGH), cctkGH);
      CCTK_CallFunction(function, attribute, const_cast<cGH *>(local_cctkGH));
    });
    synchronize();
    break;

  case mode_t::meta:
  case mode_t::global:
  case mode_t::level:
    // Call function just once
    // Note: meta mode scheduling must continue to work even after we
    // shut down ourselves!
    CCTK_CallFunction(function, attribute, cctkGH);
    break;

  default:
    assert(0);
  }

  // Check checksums
  if (poison_undefined_values)
    check_checksums(checksums, [attribute, cctkGH]() {
      std::ostringstream buf;
      buf << "CallFunction iteration " << cctkGH->cctk_iteration << " "
          << attribute->where << ": " << attribute->thorn
          << "::" << attribute->routine << " checking output";
      return buf.str();
    });

  // Mark output variables as having valid data
  {
    const std::vector<clause_t> &writes =
        decode_clauses(attribute, rdwr_t::write);
    for (const auto &wr : writes) {
      if (CCTK_GroupTypeI(wr.gi) == CCTK_GF) {
        const auto &patchdata0 = ghext->patchdata.at(0);
        const auto &leveldata0 = patchdata0.leveldata.at(0);
        const auto &groupdata0 = *leveldata0.groupdata.at(wr.gi);
        const nan_handling_t nan_handling = groupdata0.do_checkpoint
                                                ? nan_handling_t::forbid_nans
                                                : nan_handling_t::allow_nans;

        active_levels->loop_serially([&](auto &restrict leveldata) {
          auto &restrict groupdata = *leveldata.groupdata.at(wr.gi);
          const valid_t &provided = wr.valid;
          groupdata.valid.at(wr.tl).at(wr.vi).set_valid(
              provided,
              [iteration = cctkGH->cctk_iteration, where = attribute->where,
               thorn = attribute->thorn, routine = attribute->routine] {
                std::ostringstream buf;
                buf << "CallFunction iteration " << iteration << " " << where
                    << ": " << thorn << "::" << routine
                    << ": Mark output variables as valid";
                return buf.str();
              });
        });
        check_valid_gf(*active_levels, wr.gi, wr.vi, wr.tl, nan_handling,
                       [attribute, cctkGH]() {
                         std::ostringstream buf;
                         buf << "CallFunction iteration "
                             << cctkGH->cctk_iteration << " "
                             << attribute->where << ": " << attribute->thorn
                             << "::" << attribute->routine
                             << " checking output";
                         return buf.str();
                       });
      } else { // CCTK_ARRAY or CCTK_SCALAR
        auto &restrict arraygroupdata =
            *ghext->globaldata.arraygroupdata.at(wr.gi);
        const nan_handling_t nan_handling = arraygroupdata.do_checkpoint
                                                ? nan_handling_t::forbid_nans
                                                : nan_handling_t::allow_nans;
        const valid_t &provided = wr.valid;
        arraygroupdata.valid.at(wr.tl).at(wr.vi).set_valid(
            provided,
            [iteration = cctkGH->cctk_iteration, where = attribute->where,
             thorn = attribute->thorn, routine = attribute->routine] {
              std::ostringstream buf;
              buf << "CallFunction iteration " << iteration << " " << where
                  << ": " << thorn << "::" << routine
                  << ": Mark output variables as valid";
              return buf.str();
            });
        check_valid_ga(
            wr.gi, wr.vi, wr.tl, nan_handling, [attribute, cctkGH]() {
              std::ostringstream buf;
              buf << "CallFunction iteration " << cctkGH->cctk_iteration << " "
                  << attribute->where << ": " << attribute->thorn
                  << "::" << attribute->routine << " checking output";
              return buf.str();
            });
      }
    }
  }

  // Mark invalid variables as having invalid data
  {
    const std::vector<clause_t> &invalids =
        decode_clauses(attribute, rdwr_t::invalid);
    for (const auto &inv : invalids) {
      if (CCTK_GroupTypeI(inv.gi) == CCTK_GF) {
        const auto &patchdata0 = ghext->patchdata.at(0);
        const auto &leveldata0 = patchdata0.leveldata.at(0);
        const auto &groupdata0 = *leveldata0.groupdata.at(inv.gi);
        const nan_handling_t nan_handling = groupdata0.do_checkpoint
                                                ? nan_handling_t::forbid_nans
                                                : nan_handling_t::allow_nans;

        active_levels->loop_serially([&](auto &restrict leveldata) {
          auto &restrict groupdata = *leveldata.groupdata.at(inv.gi);
          const valid_t &invalidated = inv.valid;
          groupdata.valid.at(inv.tl).at(inv.vi).set_invalid(
              invalidated,
              [iteration = cctkGH->cctk_iteration, where = attribute->where,
               thorn = attribute->thorn, routine = attribute->routine] {
                std::ostringstream buf;
                buf << "CallFunction iteration " << iteration << " " << where
                    << ": " << thorn << "::" << routine
                    << ": Mark invalid variables as invalid";
                return buf.str();
              });
        });
        check_valid_gf(*active_levels, inv.gi, inv.vi, inv.tl, nan_handling,
                       [attribute, cctkGH]() {
                         std::ostringstream buf;
                         buf << "CallFunction iteration "
                             << cctkGH->cctk_iteration << " "
                             << attribute->where << ": " << attribute->thorn
                             << "::" << attribute->routine
                             << " checking output";
                         return buf.str();
                       });
      } else { // CCTK_ARRAY or CCTK_SCALAR
        auto &restrict arraygroupdata =
            *ghext->globaldata.arraygroupdata.at(inv.gi);
        const nan_handling_t nan_handling = arraygroupdata.do_checkpoint
                                                ? nan_handling_t::forbid_nans
                                                : nan_handling_t::allow_nans;
        const valid_t &invalidated = inv.valid;
        arraygroupdata.valid.at(inv.tl).at(inv.vi).set_invalid(
            invalidated,
            [iteration = cctkGH->cctk_iteration, where = attribute->where,
             thorn = attribute->thorn, routine = attribute->routine] {
              std::ostringstream buf;
              buf << "CallFunction iteration " << iteration << " " << where
                  << ": " << thorn << "::" << routine
                  << ": Mark invalid variables as invalid";
              return buf.str();
            });
        check_valid_ga(
            inv.gi, inv.vi, inv.tl, nan_handling, [attribute, cctkGH]() {
              std::ostringstream buf;
              buf << "CallFunction iteration " << cctkGH->cctk_iteration << " "
                  << attribute->where << ": " << attribute->thorn
                  << "::" << attribute->routine << " checking output";
              return buf.str();
            });
      }
    }
  }

  constexpr int didsync = 0;
  return didsync;
}

bool sync_active = false; // Catch recursive calls

struct mark_sync_active {
  mark_sync_active() {
    if (sync_active)
      CCTK_ERROR(
          "Recursive call to SyncGroupsByDirI. Maybe you are syncing grid "
          "functions in the \"restrict\" bin while the parameter "
          "\"restrict_during_sync\" is true?");
    sync_active = true;
  }
  ~mark_sync_active() { sync_active = false; }
};

int SyncGroupsByDirI(const cGH *restrict cctkGH, int numgroups,
                     const int *groups0, const int *directions) {
  DECLARE_CCTK_PARAMETERS;

  assert(in_global_mode(cctkGH));

  mark_sync_active marked;

  static Timer timer("Sync");
  Interval interval(timer);

  assert(cctkGH);
  assert(numgroups >= 0);
  assert(groups0);

  if (verbose) {
    std::ostringstream buf;
    for (int n = 0; n < numgroups; ++n) {
      if (n != 0)
        buf << ", ";
      buf << CCTK_FullGroupName(groups0[n]);
    }
#pragma omp critical
    CCTK_VINFO("SyncGroups %s", buf.str().c_str());
  }

  const int gi_regrid_error = CCTK_GroupIndex("CarpetXRegrid::regrid_error");
  assert(gi_regrid_error >= 0);

  std::vector<int> groups;
  for (int n = 0; n < numgroups; ++n) {
    const int gi = groups0[n];
    if (CCTK_GroupTypeI(gi) != CCTK_GF)
      continue;
    // Don't restrict the regridding error
    if (gi == gi_regrid_error)
      continue;
    groups.push_back(gi);
  }

  // AMR-C7.  The sync is the scope the ordering claim is about, so it is marked
  // at entry and at exit and everything between the two markers belongs to it.
  if (log_sched_on()) {
    log_sched_flux_census(cctkGH);
    std::ostringstream fields;
    fields << "site=sync/enter numgroups=" << numgroups
           << " gf_groups=" << groups.size() << " nlevels="
           << ghext->num_levels() << " npatches=" << ghext->num_patches()
           << " restrict_during_sync=" << int(bool(restrict_during_sync))
           << " presync_mode=" << presync_mode << " levels=L["
           << (active_levels ? active_levels->min_level : -1) << ","
           << (active_levels ? active_levels->max_level : -1) << ")P["
           << (active_levels ? active_levels->min_patch : -1) << ","
           << (active_levels ? active_levels->max_patch : -1) << ")";
    log_sched(cctkGH, fields.str());
  }

  // Skip groups that have valid ghosts and boundaries
  const int n_groups_before_presync = groups.size();
  if (CCTK_EQUALS(presync_mode, "presync-only")) {
    active_levels->loop_serially([&](auto &restrict leveldata) {
      std::vector<int> new_groups;
      for (const int gi : groups) {
        auto &restrict groupdata = *leveldata.groupdata.at(gi);
        bool need_sync = false;
        for (int tl = 0; tl < int(groupdata.valid.size()); tl++) {
          if (need_sync)
            break;
          auto &timeleveldata = groupdata.valid.at(tl);
          for (int vi = 0; vi < int(timeleveldata.size()); vi++) {
            if (need_sync)
              break;
            valid_t have = groupdata.valid.at(tl).at(vi).get();
            if (!have.valid_ghosts || !have.valid_outer) {
              need_sync = true;
            }
          }
        }
        if (need_sync) {
          new_groups.push_back(gi);
        }
      }
      groups = new_groups;
    });
    // AMR-C7.  THIS FILTER SITS ABOVE THE RESTRICTION, AND THAT IS NOT
    // OBVIOUSLY RIGHT.  It drops a group whose ghosts and outer boundaries are
    // already valid -- a question about the GHOST zones -- and what it drops it
    // drops from the `Restrict` call below as well, which is a question about
    // the coarse INTERIOR under a refined level.  When it empties the list the
    // whole sync returns here and no level is restricted at all.  Every
    // parameter file in this project sets `presync-only`, so this is the live
    // path; the counts are measured in `evidence/amr/c7/`.  This step reports
    // it and does not change it (R1).
    if (log_sched_on() && groups.size() != size_t(n_groups_before_presync)) {
      std::ostringstream fields;
      fields << "site=sync/presync-filter groups_before="
             << n_groups_before_presync << " groups_after=" << groups.size()
             << " empties_the_sync=" << int(groups.size() == 0);
      log_sched(cctkGH, fields.str());
    }
    if (groups.size() == 0) {
      if (log_sched_on()) {
        std::ostringstream fields;
        fields << "site=sync/presync-skip groups_before="
               << n_groups_before_presync
               << " restricted_levels=0 mpinterp_calls=0";
        log_sched(cctkGH, fields.str());
        // AND THE BRACKET CLOSES HERE TOO.  An `enter` that is not always
        // matched by an `exit` makes a reader pair this sync's `enter` with the
        // NEXT sync's `exit` and reason about a window that is two syncs wide.
        // That is `[P422]`'s defect one level over -- a scope that merges two
        // traversals -- and it is the defect this whole instrument exists to
        // avoid, so every return path from here on emits the closing line.
        std::ostringstream ex;
        ex << "site=sync/exit gf_groups=0 via=presync-skip";
        log_sched(cctkGH, ex.str());
      }
      return 0;
    }
  }

  if (restrict_during_sync) {
    long n_restrict_calls = 0;
    active_levels->loop_fine_to_coarse([&](const auto &leveldata) {
      if (leveldata.level < ghext->num_levels() - 1) {
        ++n_restrict_calls;
        Restrict(cctkGH, leveldata.level, groups, "sync");
      }
    });
    // `[P135]`: a run in which nothing was restricted has to SAY so, or it
    // cannot be told from a run in which the instrument was not compiled in.
    if (log_sched_on()) {
      std::ostringstream fields;
      fields << "site=sync/restrict-loop calls=" << n_restrict_calls
             << " nlevels=" << ghext->num_levels();
      log_sched(cctkGH, fields.str());
    }
    // FIXME: cannot call POSTRESTRICT since this could contain a SYNC leading
    // to an infinite loop. This means that outer boundaries will be left
    // invalid after an implicit restrict
    // CCTK_Traverse(cctkGH, "CCTK_POSTRESTRICT");
  }

  static const bool have_multipatch_boundaries =
      CCTK_IsFunctionAliased("MultiPatch_Interpolate");

  // mp_corners_7.md section 5 point 4: always-on (not CCTK_DEBUG-gated) log,
  // tagged with group name and a monotonic call counter, to confirm the
  // relative ordering of this per-group BC sweep against
  // MultiPatch1_Interpolate's own per-call log (see interpolate.cxx).
  {
    // BUGFIX_TODO.md B10: the increment used to sit OUTSIDE the `verbose` test,
    // so an unsynchronised mutable static was written on every sync of every
    // run while being read on none of them.  Moved inside.  The PRINTED numbers
    // do not move: the only path that reads the counter is the one that also
    // increments it, so any evidence quoting a call number (`[E3]`'s 32 NaNs at
    // "SyncGroupsByDirI call #2") was taken with `verbose` on and still reads
    // the same.
    if (verbose) {
      static long call_counter = 0;
      ++call_counter;
      for (const int gi : groups) {
#pragma omp critical
        CCTK_VINFO("SyncGroupsByDirI call #%ld: group %s", call_counter,
                   CCTK_FullGroupName(gi));
      }
    }
  }

  // Check preconditions
  for (const int gi : groups) {
    const auto &patchdata0 = ghext->patchdata.at(0);
    const auto &leveldata0 = patchdata0.leveldata.at(0);
    const auto &groupdata0 = *leveldata0.groupdata.at(gi);
    const nan_handling_t nan_handling = groupdata0.do_checkpoint
                                            ? nan_handling_t::forbid_nans
                                            : nan_handling_t::allow_nans;
    // We always sync all directions.
    // If there is more than one time level, then we don't sync the
    // oldest.
    // TODO: during evolution, sync only one time level
    const int ntls0 = groupdata0.mfab.size();
    const int sync_tl0 = ntls0 > 1 ? ntls0 - 1 : ntls0;

    // BUGFIX_TODO.md step B8b ([O4]): A MULTI-PATCH SYNC OF A GROUP WITH MORE
    // THAN ONE SYNCED TIME LEVEL IS REFUSED.
    //
    // THE HOLE.  `MultiPatch_Interpolate` writes the interpatch ghost zones at
    // `tl = 0` only.  An interpatch face carries no stored outer boundary
    // condition, so the boundary passes write nothing there either, at any
    // time level.  Yet the postcondition below marks ghosts and outer
    // boundaries VALID for every `tl < sync_tl0`.  At `sync_tl0 > 1` that
    // certifies, as valid, cells this sync has just poisoned and nothing has
    // written -- which is precisely the charge this project is making against
    // the code it is fixing.  Neither state we could reach here is one we are
    // willing to certify: before an interpatch face stopped carrying an outer
    // BC, `tl >= 1` was written with an outer-boundary echo across a seam
    // instead.
    //
    // WHY `sync_tl0 > 1` AND NOT `ntls0 > 1`.  Every time-level loop in this
    // function -- the preconditions, the two `FillPatch_*` call sites, the
    // postcondition, the debug NaN sweep -- is bounded by
    // `sync_tl = ntls > 1 ? ntls - 1 : ntls`, the expression two lines above.
    // A TWO-time-level group therefore syncs `tl = 0` alone: `tl >= 1` is
    // never poisoned here, never filled here, and never marked valid here, so
    // it has no hole in it and refusing it would refuse a configuration that
    // works.  The predicate is the number of time levels this function
    // actually touches, not the number the group declares.
    //
    // WHY THE PATCH COUNT AND NOT `have_multipatch_boundaries`.  The hole needs
    // an interpatch face.  A `patch_system = "Cartesian"` grid is a single
    // patch with six physical outer faces and none, so the `all` pass fills its
    // boundaries at every time level exactly as on a non-patch grid; the alias
    // is nonetheless aliased there.  Same reason as the PARAMCHECK guard in
    // `driver.cxx`.
    //
    // UNEXERCISED, AND SAID SO RATHER THAN LEFT TO BE NOTICED.  A group's time
    // level count is compile time (`driver.cxx`, `mfab.resize`), and no thorn
    // in this thorn list declares a `CCTK_GF` group with more than two.  The
    // guard was reached and its message read on a purpose-built configuration
    // (`evidence/fix/b8/tl3/`), not on any rig that ships.  Note what that
    // witness also shows: without this check the same configuration already
    // stops one loop later, at the `error_if_invalid` below, complaining that
    // `tl = 1` is invalid.  This moves the refusal and names the reason; it
    // does not turn a silently wrong run into a refused one.
    if (ghext->num_patches() > 1 && sync_tl0 > 1)
      CCTK_VERROR(
          "SyncGroupsByDirI: group %s has %d time levels, of which %d are "
          "synchronised, on a multi-patch system (%d patches). Only the "
          "current time level of a synchronised group is filled across patch "
          "boundaries; the older ones would be marked valid without anything "
          "having written their interpatch ghost zones. Refusing rather than "
          "certifying them. Declare the group with at most 2 time levels, or "
          "use a single-patch grid.",
          groupdata0.groupname.c_str(), ntls0, sync_tl0, ghext->num_patches());

    active_levels->loop_serially([&](auto &restrict leveldata) {
      auto &restrict groupdata = *leveldata.groupdata.at(gi);

      if (leveldata.level > 0) {

        const int level = leveldata.level;
        const auto &restrict coarseleveldata =
            ghext->patchdata.at(leveldata.patch).leveldata.at(level - 1);
        auto &restrict coarsegroupdata = *coarseleveldata.groupdata.at(gi);
        assert(coarsegroupdata.numvars == groupdata.numvars);

        for (int tl = 0; tl < sync_tl0; ++tl) {
          for (int vi = 0; vi < groupdata.numvars; ++vi) {
            error_if_invalid(coarsegroupdata, vi, tl, make_valid_int(), []() {
              return "SyncGroupsByDirI on coarse level before prolongation";
            });
          }
        } // for tl

      } // if leveldata.level > 0

      for (int tl = 0; tl < sync_tl0; ++tl) {
        for (int vi = 0; vi < groupdata.numvars; ++vi) {
          // Synchronization only uses the interior
          error_if_invalid(groupdata, vi, tl, make_valid_int(),
                           []() { return "SyncGroupsByDirI before syncing"; });
          groupdata.valid.at(tl).at(vi).set_invalid(make_valid_ghosts(), []() {
            return "SyncGroupsByDirI before syncing: "
                   "Mark ghost zones as invalid";
          });
        }
      } // for tl
    });

    active_levels_t active_fine_levels = *active_levels;
    using std::max;
    active_fine_levels.min_level = max(active_fine_levels.min_level, 1);
    for (int tl = 0; tl < sync_tl0; ++tl) {
      for (int vi = 0; vi < groupdata0.numvars; ++vi) {
        check_valid_gf(active_fine_levels, gi, vi, tl, nan_handling, []() {
          return "SyncGroupsByDirI on coarse level before prolongation";
        });
        poison_invalid_gf(*active_levels, gi, vi, tl);
        check_valid_gf(*active_levels, gi, vi, tl, nan_handling,
                       []() { return "SyncGroupsByDirI before syncing"; });
      }
    } // for tl
  } // for gi

  // BUGFIX_TODO.md step B3 deleted the bootstrap MultiPatch_Interpolate that
  // used to run here, together with the level-0 FillBoundary pre-pass that fed
  // it and the 14-line mp_corners_7.md comment that justified it. Step B2 made
  // pass 1 skip the interpatch x outer corners, so it no longer reads an
  // interpatch ghost at all, and making those ghosts finite ahead of it -- the
  // bootstrap's only purpose -- is no longer anything. The surviving call is the
  // one below, after tasks1/2/3 have run: by then AMReX has filled the
  // intra-patch (box-split) ghosts itself, which is what the FillBoundary
  // pre-pass was standing in for, and pass 1 has written the outer ghosts.
  // There is now exactly ONE interpolate call per sync, so no call can read
  // another call's output within a sync -- the two-pass non-idempotency behind
  // the 696 + 312 corruption is absent by construction rather than gated off.
  // We need to loop over groups, patches, and levels in a definite
  // order so that AMReX's communication pattern does not get
  // confused. Therefore all the loops here are serial. The only
  // parallelization happens within AMReX and within our boundary
  // conditions. This is not efficient.

  task_manager tasks1;
  task_manager tasks2;
  task_manager tasks3;

  for (const int gi : groups) {
    active_levels->loop_serially([&](auto &restrict leveldata) {
      auto &restrict groupdata = *leveldata.groupdata.at(gi);

      // We always sync all directions.
      // If there is more than one time level, then we don't sync the
      // oldest.
      // TODO: during evolution, sync only one time level
      const int ntls = groupdata.mfab.size();
      const int sync_tl = ntls > 1 ? ntls - 1 : ntls;

      if (leveldata.level == 0) {
        // Copy from adjacent boxes on same level

        for (int tl = 0; tl < sync_tl; ++tl) {
          // BUGFIX_TODO.md step B2(c). The pass is chosen PER TIME LEVEL, not
          // per call site. Pass 1 runs once per synced time level; the
          // interpolator and pass 2 both run at `tl = 0` only. A flat
          // `skip_interpatch_corners` here would leave every `tl >= 1`
          // interpatch corner written by nobody and then marked valid.
          //
          // Gating on `have_multipatch_boundaries` as well is not decoration:
          // it makes a single-patch run pass `bc_pass_t::all` literally, so
          // the "this commit is inert on a single patch" claim is a statement
          // about which code path executes rather than about a boolean.
          //
          // What this ternary does NOT do, and B8 is what does: step B7 has
          // now removed the outer BC from `groupdata.boundaries` on interpatch
          // faces, so `all` no longer writes the interpatch ghost at
          // `tl >= 1` while the validity marks below still claim it -- the
          // `tl` axis is closed by B8's refusal of more than one time level,
          // not by this line. Do not read a correctness into it that it does
          // not have. (No group in any rig here declares more than one time
          // level, so this is latent rather than live.)
          const bc_pass_t bc_pass = (have_multipatch_boundaries && tl == 0)
                                        ? bc_pass_t::skip_interpatch_corners
                                        : bc_pass_t::all;
          tasks1.submit_serially([&tasks2, &leveldata, &groupdata, tl,
                                  bc_pass]() {
            FillPatch_Sync(tasks2, groupdata, *groupdata.mfab.at(tl),
                           ghext->patchdata.at(leveldata.patch)
                               .amrcore->Geom(leveldata.level),
                           bc_pass);
          });
        } // for tl

      } else { // if leveldata.level > 0
        // Copy from adjacent boxes on same level, and interpolate
        // from next coarser level

        const int level = leveldata.level;
        const auto &restrict coarseleveldata =
            ghext->patchdata.at(leveldata.patch).leveldata.at(level - 1);
        auto &restrict coarsegroupdata = *coarseleveldata.groupdata.at(gi);
        assert(coarsegroupdata.numvars == groupdata.numvars);

        amrex::Interpolater *const interpolator = groupdata.interpolator;

        for (int tl = 0; tl < sync_tl; ++tl) {

          // B2(c), same selection as the level-0 branch above. It reaches the
          // two calls in `FillPatch_ProlongateGhosts` that fill the real
          // `mfab`; the coarse temporary inside it is pinned to
          // `bc_pass_t::all` there, because nothing ever runs a second pass
          // over a temporary that `FillPatchInterp` is about to read.
          //
          // After B7 that `all` no longer covers the temporary's INTERPATCH
          // faces -- nothing does -- so the prolongation reads what
          // `mf_set_domain_bndry` left there ([P22]). That is now refused at
          // the fill by `check_camr_contract`, under the C-AMR contract stated
          // at the top of `fillpatch.cxx`, rather than by refusing multipatch
          // + AMR outright; see also the dispatch comment in
          // `boundaries_impl.hxx`.
          const bc_pass_t bc_pass = (have_multipatch_boundaries && tl == 0)
                                        ? bc_pass_t::skip_interpatch_corners
                                        : bc_pass_t::all;
          tasks1.submit_serially([&tasks2, &tasks3, &leveldata, &groupdata,
                                  &coarsegroupdata, interpolator, tl,
                                  bc_pass]() {
            FillPatch_ProlongateGhosts(tasks2, tasks3, groupdata,
                                       coarsegroupdata, *groupdata.mfab.at(tl),
                                       *coarsegroupdata.mfab.at(tl),
                                       ghext->patchdata.at(leveldata.patch)
                                           .amrcore->Geom(leveldata.level),
                                       ghext->patchdata.at(leveldata.patch)
                                           .amrcore->Geom(leveldata.level - 1),
                                       interpolator, groupdata.bcrecs, bc_pass);
          });

        } // for tl

      } // if leveldata.level > 0
    });
  } // for gi

  tasks1.run_tasks_serially();
  synchronize();
  tasks2.run_tasks_serially();
  synchronize();
  tasks3.run_tasks_serially();
  synchronize();

  // Check postconditions
  for (const int gi : groups) {
    const auto &patchdata0 = ghext->patchdata.at(0);
    const auto &leveldata0 = patchdata0.leveldata.at(0);
    const auto &groupdata0 = *leveldata0.groupdata.at(gi);
    const nan_handling_t nan_handling = groupdata0.do_checkpoint
                                            ? nan_handling_t::forbid_nans
                                            : nan_handling_t::allow_nans;
    // We always sync all directions.
    // If there is more than one time level, then we don't sync the
    // oldest.
    // TODO: during evolution, sync only one time level
    const int ntls0 = groupdata0.mfab.size();
    const int sync_tl0 = ntls0 > 1 ? ntls0 - 1 : ntls0;

    // For the multipatch case, corner ghost cells at the outer+interpatch face
    // intersection are not yet written at this point.  Since BUGFIX_TODO.md
    // step B2 the first BC pass (inside FillPatch_Sync /
    // FillPatch_ProlongateGhosts) runs `bc_pass_t::skip_interpatch_corners` at
    // tl = 0 and skips them outright, rather than writing them from an
    // unpopulated interpatch ghost zone as it used to.  MultiPatch_Interpolate
    // fills the interpatch ghosts and the corners-only second BC pass (below)
    // then writes the corners, once.  Defer the validity marks and poison/check
    // calls to after that second pass so that validity flags never claim corner
    // cells are valid while they are still unwritten.
    if (!have_multipatch_boundaries) {
      active_levels->loop_serially([&](auto &restrict leveldata) {
        auto &restrict groupdata = *leveldata.groupdata.at(gi);

        for (int tl = 0; tl < sync_tl0; ++tl) {
          for (int vi = 0; vi < groupdata.numvars; ++vi) {
            groupdata.valid.at(tl).at(vi).set_ghosts(true, []() {
              return "SyncGroupsByDirI after syncing: "
                     "Mark ghost zones as valid";
            });
            if (groupdata.all_faces_have_symmetries_or_boundaries())
              groupdata.valid.at(tl).at(vi).set_outer(true, []() {
                return "SyncGroupsByDirI after syncing: "
                       "Mark outer boundaries as valid";
              });
          }
        } // for tl
      });

      for (int tl = 0; tl < sync_tl0; ++tl) {
        for (int vi = 0; vi < groupdata0.numvars; ++vi) {
          poison_invalid_gf(*active_levels, gi, vi, tl);
          check_valid_gf(*active_levels, gi, vi, tl, nan_handling,
                         []() { return "SyncGroupsByDirI after syncing"; });
        }
      } // for tl
    }
  } // for gi

  if (have_multipatch_boundaries) {
    std::vector<CCTK_INT> cactusvarinds;
    for (int group : groups) {
      const auto &groupdata =
          *ghext->patchdata.at(0).leveldata.at(0).groupdata.at(group);
      for (int var = 0; var < groupdata.numvars; ++var)
        cactusvarinds.push_back(groupdata.firstvarindex + var);
    }
#ifdef CCTK_DEBUG
    log_mp_interpolate_call("SyncGroupsByDirI main", cactusvarinds);
#endif
    // AMR-C7.  THE INTERPATCH FILL, bracketed.  `CapyrX_MultiPatch`'s own
    // `MPLEVELS` line (`CAPYRX_LOG_LEVELS`) lands between these two markers and
    // is the independent, second-repo reading of the same ordering.
    if (log_sched_on()) {
      std::ostringstream fields;
      fields << "site=sync/mpinterp-pre nvars=" << cactusvarinds.size()
             << " ngroups=" << groups.size();
      log_sched(cctkGH, fields.str());
    }
    // The sync's only interpolation pass (BUGFIX_TODO.md step B3 deleted the
    // bootstrap one above), so slave_overlap's write-back is unconditional
    // again: nothing later in this sync reads these interior cells as a donor,
    // and no earlier call in it can have read them either.
    MultiPatch_Interpolate(cctkGH, cactusvarinds.size(),
                           cactusvarinds.data());
    if (log_sched_on()) {
      std::ostringstream fields;
      fields << "site=sync/mpinterp-post nvars=" << cactusvarinds.size();
      log_sched(cctkGH, fields.str());
    }

    // AMR-B4b (`[P123]`): the interpolator has just overwritten interior
    // cells; their inter-box ghost copies are one write stale until this runs.
    // See `reshare_interior_writes` above for why this is safe against `I` and
    // `O`, and why it is here rather than after the second BC pass.
    reshare_interior_writes(groups);

    // Second BC pass: correct corner ghost cells at the outer+interpatch face
    // intersection.
    //
    // Background: MultiPatch_Interpolate fills interpatch ghost cells by
    // interpolating from neighbouring patches.  However it skips any ghost
    // cell where p.NI[d] != 0 in a direction that is an outer boundary
    // (see loop_bnd skip logic in CapyrX_MultiPatch/src/interpolate.cxx).
    // Those "corner" cells — simultaneously in an interpatch ghost zone in
    // one direction and on an outer-BC face in another — are therefore never
    // touched by MultiPatch_Interpolate.
    //
    // Until BUGFIX_TODO.md step B2, the first BC pass (inside FillPatch_Sync /
    // FillPatch_ProlongateGhosts) DID write these corner cells, using
    // interpatch ghost sources that were not yet populated (MultiPatch had not
    // run).  For non-Dirichlet BCs (Neumann, Robin, linear extrapolation) that
    // produces wrong values: the stencil source src[d] = dst[d] is in the
    // interpatch ghost zone and was NaN/stale.  This pass then overwrote them
    // -- a second write of every corner, and of every pure outer ghost too.
    //
    // Now that MultiPatch_Interpolate has filled all pure interpatch ghost
    // cells, their values are valid sources, and this pass computes the corner
    // cell values from them.
    //
    // NOTE: This second pass must happen AFTER MultiPatch_Interpolate and
    // BEFORE the validity checks below.
    //
    // BUGFIX_TODO.md step B2(d): this pass is now CORNERS-ONLY, and runs at
    // `tl = 0` only.
    //
    //   - corners-only, because it used to be a full `all` pass: pass 1 had
    //     already written every pure outer ghost `O`, and this pass wrote all
    //     of them again.  That is the double write the objection names.  Pass 1
    //     now writes `O` and skips the corners; this pass writes the corners
    //     and skips everything else; the two are disjoint and every cell is
    //     written exactly once after the data motion.
    //
    //     MEASURED COST OF REMOVING THE SECOND `O` WRITE, and it is not zero
    //     (evidence/fix/b2/README.md, finding [P45]).  The second write was
    //     redundant only while nothing changed the interior between the two
    //     passes.  With `CapyrX_MultiPatch::slave_overlap = yes`,
    //     `MultiPatch_Interpolate` writes INTERIOR overlap cells, so an outer
    //     BC that reads the interior -- neumann, linear_extrapolation, robin,
    //     but not dirichlet, which reads nothing -- gave a different answer the
    //     second time, and that answer was the one consistent with the final
    //     interior.  After this commit `O` is computed from the PRE-slave
    //     interior and is stale by one interpolation error: on the battery's
    //     analytic smooth field, max |delta| falls 1.0e-1 -> 4.8e-3 -> 3.8e-4
    //     -> 3.9e-5 -> 1.4e-5 as interpolation_order goes 0 -> 4, and a Neumann
    //     leg goes from 0 to 1792 outer ghosts (of 7500) that no longer equal
    //     the interior plane they copy.  Every slave-OFF leg is bit-identical.
    //     This is a property of the target ordering in BUGFIX_TODO.md section
    //     0, whose DAG has the interpolator writing only `I`; with
    //     `slave_overlap` the real graph has an interior -> `O` edge and is
    //     cyclic through Channel 1.  It is recorded, not fixed here: removing
    //     the second write is exactly what this step is for, and re-adding it
    //     conditionally would be a new mechanism.
    //
    //   - `tl = 0` only, because `MultiPatch_Interpolate` writes `tl = 0`. A
    //     corners-only pass at `tl >= 1` would read an interpatch ghost that
    //     nobody has filled -- exactly the read this partition exists to
    //     remove. At `tl >= 1` pass 1 runs `all` and keeps doing what `main`
    //     does; see the ternary at the `FillPatch_*` call sites above, and B8
    //     for why that is refused rather than relied on.
    //
    // Nothing else moves. The `set_ghosts` / `set_outer` / `poison_invalid_gf`
    // trio is a unit (moving the `set_*` calls without the poison call
    // silently destroys the just-communicated inter-box ghost data); on this
    // path it has already been split, and repairing that is NOT part of this
    // commit.
    active_levels->loop_serially([&](auto &restrict leveldata) {
      for (const int gi : groups) {
        auto &restrict groupdata = *leveldata.groupdata.at(gi);
        groupdata.apply_boundary_conditions(*groupdata.mfab.at(0),
                                            bc_pass_t::interpatch_corners_only);
      }
    });

    // Corner cells are now written.  Mark ghost zones and outer boundaries
    // valid for all groups across all patches and levels.  This is deferred
    // from the postcondition loop above because the first BC pass leaves the
    // corner cells at outer+interpatch junctions UNWRITTEN until this point
    // (before BUGFIX_TODO.md step B2 it left them written from an unpopulated
    // source, which is why this comment used to say "stale" rather than
    // "unwritten").
    for (const int gi : groups) {
      const auto &groupdata0 =
          *ghext->patchdata.at(0).leveldata.at(0).groupdata.at(gi);
      const int ntls0 = groupdata0.mfab.size();
      const int sync_tl0 = ntls0 > 1 ? ntls0 - 1 : ntls0;
      active_levels->loop_serially([&](auto &restrict leveldata) {
        auto &restrict groupdata = *leveldata.groupdata.at(gi);
        for (int tl = 0; tl < sync_tl0; ++tl) {
          for (int vi = 0; vi < groupdata.numvars; ++vi) {
            groupdata.valid.at(tl).at(vi).set_ghosts(true, []() {
              return "SyncGroupsByDirI after 2nd BC pass: "
                     "Mark ghost zones as valid";
            });
            if (groupdata.all_faces_have_symmetries_or_boundaries())
              groupdata.valid.at(tl).at(vi).set_outer(true, []() {
                return "SyncGroupsByDirI after 2nd BC pass: "
                       "Mark outer boundaries as valid";
              });
          }
        }
      });
    }

#ifdef CCTK_DEBUG
    // Verify that the 2nd BC pass zeroed out NaN values from ghost zones for
    // groups that forbid NaNs (do_checkpoint=yes).  A non-zero count after the
    // 2nd pass indicates the corner-cell fix is incomplete.
    active_levels->loop_serially([&](auto &restrict leveldata) {
      for (const int gi : groups) {
        const auto &restrict groupdata = *leveldata.groupdata.at(gi);
        if (!groupdata.do_checkpoint)
          continue; // allow NaNs in non-checkpointed groups
        const int ntls = groupdata.mfab.size();
        const int sync_tl = ntls > 1 ? ntls - 1 : ntls;
        for (int tl = 0; tl < sync_tl; ++tl) {
          const amrex::MultiFab &mf = *groupdata.mfab.at(tl);
          // MultiFab::contains_nan() scans all cells including ghost zones.
          if (mf.contains_nan()) {
#pragma omp critical
            CCTK_VERROR(
                "CCTK_DEBUG SyncGroupsByDirI: After 2nd BC pass + "
                "MultiPatch_Interpolate, group '%s' patch %d level %d tl=%d "
                "still contains NaN. The 2nd BC pass should have cleared all "
                "corner-cell NaNs sourced from valid interpatch ghost cells. "
                "Remaining NaN indicates an unresolved ghost-zone bug.",
                groupdata.groupname.c_str(), leveldata.patch, leveldata.level,
                tl);
          }
        }
      }
    });
#endif // CCTK_DEBUG

    // WHY THERE IS NO `poison_invalid_gf` HERE, AND WHY ADDING ONE WOULD BE A
    // BUG (BUGFIX_TODO.md step B5, withdrawn; `multipatch_case.md` [C8]/[N2]).
    //
    // The non-multipatch postcondition above runs `set_* -> poison_invalid_gf
    // -> check_valid_gf`; this path runs `set_* -> check_valid_gf` and skips
    // the poison call. That asymmetry looks like an oversight and has been
    // proposed as a symmetry repair twice. It is not one: on a patch system
    // the added call is a no-op where it is safe and a data destroyer where it
    // is not.
    //
    //   - `poison_invalid_gf` poisons `where_t::boundary` whenever
    //     `valid_outer` is false (`valid.cxx:198-202`).
    //   - `where_t::boundary` is not the physical outer boundary here. It is
    //     whatever `cctk_bbox` marks, and `cctk_bbox` comes from
    //     `GridDesc::GridDesc(leveldata, mfp)`, which sets
    //     `bbox[f][d] = vbx[...] == domain[...]` with NO symmetry test (`:178`,
    //     reaching `cctk_bbox` at `:698`). On a patch system every interpatch
    //     face IS a patch-domain face, so `where_t::boundary` is the whole of
    //     outer ghosts + interpatch ghosts + their corners, and
    //     `where_t::ghosts` is only the inter-box ghosts. (The one bbox
    //     definition in the tree that does test
    //     `symmetries[f][d] != symmetry_t::none` is the sibling constructor
    //     `GridDesc::GridDesc(leveldata, global_component)` -- `:235`, its
    //     bbox at `:285` -- and it has no callers, in this thorn or any
    //     other.)
    //   - `set_ghosts(true)` therefore does not cover the interpatch ghosts.
    //     The *outer* bit governs them, and `set_outer(true)` above is
    //     conditional on `all_faces_have_symmetries_or_boundaries()`
    //     (`driver.cxx:995-1005`), which is per patch AND per group.
    //
    // So with an outer boundary condition set, `set_outer(true)` has just run
    // and a poison call here does nothing. With `boundary_* = "none"` the
    // predicate is false on any patch that owns a physical outer face, and the
    // call NaN-fills the interpatch ghosts `MultiPatch_Interpolate` wrote
    // correctly a few lines above -- between that write and its next reader.
    //
    // The real defect is that CarpetX has no `O`-only validity bit: one bit
    // covers the outer boundary, the interpatch ghosts and their corners, so
    // the interpolator's own output lives inside the region the poison call
    // would clear. That is `main`'s defect, not this branch's, and splitting
    // the bit is separate work with its own failing-before test. Until it is
    // split, this path deliberately marks validity and checks, and does not
    // poison.

    for (const int gi : groups) {
      const auto &patchdata0 = ghext->patchdata.at(0);
      const auto &leveldata0 = patchdata0.leveldata.at(0);
      const auto &groupdata0 = *leveldata0.groupdata.at(gi);
      const nan_handling_t nan_handling = groupdata0.do_checkpoint
                                              ? nan_handling_t::forbid_nans
                                              : nan_handling_t::allow_nans;
      // We always sync all directions.
      // If there is more than one time level, then we don't sync the
      // oldest.
      // TODO: during evolution, sync only one time level
      const int ntls0 = groupdata0.mfab.size();
      const int sync_tl0 = ntls0 > 1 ? ntls0 - 1 : ntls0;

      for (int tl = 0; tl < sync_tl0; ++tl)
        for (int vi = 0; vi < groupdata0.numvars; ++vi)
          check_valid_gf(*active_levels, gi, vi, tl, nan_handling,
                         []() { return "SyncGroupsByDirI after syncing"; });

    } // for gi

  } else {
    assert(ghext->num_patches() == 1);
  }

  // AMR-B4a: the sync's postcondition that nobody was checking.  It runs after
  // the second BC pass and after the validity marks, i.e. on exactly the state
  // the next reader -- and the TSV/Silo writer -- will see, which is what makes
  // its `bad_cells` column comparable with `rowdiff.py`'s multi-valued cell
  // count on the same run's output.  It is placed outside the multipatch
  // branch on purpose: the invariant is about AMReX boxes, not about patches,
  // and a single-patch AMR rig is a legitimate subject for it.
  interbox_check("SyncGroupsByDirI", groups);

  // AMR-C7: closes the scope opened by `site=sync/enter`.  Everything between
  // the two markers happened inside this sync.
  if (log_sched_on()) {
    std::ostringstream fields;
    fields << "site=sync/exit gf_groups=" << groups.size() << " via=full";
    log_sched(cctkGH, fields.str());
  }

  assert(sync_active);

  return numgroups; // number of groups synchronized
}

void Reflux(const cGH *cctkGH, int level) {
  DECLARE_CCTK_PARAMETERS;

  // AMR-C7.  The census runs on the first logged event of the run, whichever
  // it is; putting the call here as well as in `SyncGroupsByDirI` means a run
  // that never syncs still reports it.
  if (log_sched_on())
    log_sched_flux_census(cctkGH);

  if (!do_reflux) {
    // `[N13]`: "nothing was refluxed" has two causes and they are not the same
    // measurement.  Say which one this is.
    if (log_sched_on()) {
      std::ostringstream fields;
      fields << "site=evol/reflux-off level=" << level << " do_reflux=0";
      log_sched(cctkGH, fields.str());
    }
    return;
  }

  static Timer timer("Reflux");
  Interval interval(timer);

  // AMR-C7's counters.  `examined` is the (patch, GF group) pairs this call
  // walked and `with_freg` the ones that carry a flux register, which is the
  // only kind this function can move.  Both are printed, so a zero is a
  // measured zero and not a silence (`[P135]`, `[P184]`).
  long n_patches = 0, n_examined = 0, n_with_freg = 0;

  for (const auto &patchdata : ghext->patchdata) {
    if (level + 1 < int(patchdata.leveldata.size())) {
      ++n_patches;
      auto &leveldata = patchdata.leveldata.at(level);
      const auto &fineleveldata = patchdata.leveldata.at(level + 1);
      for (int gi = 0; gi < int(leveldata.groupdata.size()); ++gi) {
        const int tl = 0;
        cGroup group;
        int ierr = CCTK_GroupData(gi, &group);
        assert(!ierr);

        if (group.grouptype != CCTK_GF)
          continue;
        ++n_examined;

        auto &groupdata = *leveldata.groupdata.at(gi);
        const auto &finegroupdata = *fineleveldata.groupdata.at(gi);
        const nan_handling_t nan_handling = groupdata.do_checkpoint
                                                ? nan_handling_t::forbid_nans
                                                : nan_handling_t::allow_nans;

        // If the group has associated fluxes
        if (finegroupdata.freg) {
          ++n_with_freg;

          // Check coarse and fine data and fluxes are valid
          for (int vi = 0; vi < finegroupdata.numvars; ++vi) {
            error_if_invalid(finegroupdata, vi, tl, make_valid_int(), []() {
              return "Reflux before refluxing: Fine level data";
            });
            error_if_invalid(groupdata, vi, tl, make_valid_int(), []() {
              return "Reflux before refluxing: Coarse level data";
            });
          }
          for (int d = 0; d < dim; ++d) {
            const int flux_gi = finegroupdata.fluxes.at(d);
            const auto &flux_finegroupdata =
                *fineleveldata.groupdata.at(flux_gi);
            const auto &flux_groupdata = *leveldata.groupdata.at(flux_gi);
            for (int vi = 0; vi < finegroupdata.numvars; ++vi) {
              error_if_invalid(
                  flux_finegroupdata, vi, tl, make_valid_int(), [&]() {
                    std::ostringstream buf;
                    buf << "Reflux: Fine level flux in direction " << d;
                    return buf.str();
                  });
              error_if_invalid(flux_groupdata, vi, tl, make_valid_int(), [&]() {
                std::ostringstream buf;
                buf << "Reflux: Coarse level flux in direction " << d;
                return buf.str();
              });
            }
          }

          for (int d = 0; d < dim; ++d) {
            const int flux_gi = finegroupdata.fluxes.at(d);
            const auto &flux_finegroupdata =
                *fineleveldata.groupdata.at(flux_gi);
            const auto &flux_groupdata = *leveldata.groupdata.at(flux_gi);
            finegroupdata.freg->CrseInit(*flux_groupdata.mfab.at(tl), d, 0, 0,
                                         flux_groupdata.numvars, -1);
            finegroupdata.freg->FineAdd(*flux_finegroupdata.mfab.at(tl), d, 0,
                                        0, flux_finegroupdata.numvars, 1);
          }
          const amrex::Geometry &geom = patchdata.amrcore->Geom(level);
          finegroupdata.freg->Reflux(*groupdata.mfab.at(tl), 1.0, 0, 0,
                                     groupdata.numvars, geom);

          const active_levels_t active_levels(level, level + 1);
          for (int vi = 0; vi < finegroupdata.numvars; ++vi)
            check_valid_gf(active_levels, gi, vi, tl, nan_handling, []() {
              return "Reflux after refluxing: Fine level data";
            });
        }
      } // for gi
    } // if level exists
  } // for patchdata

  if (log_sched_on()) {
    std::ostringstream fields;
    fields << "site=evol/reflux level=" << level << " do_reflux=1"
           << " patches=" << n_patches << " examined=" << n_examined
           << " with_freg=" << n_with_freg;
    log_sched(cctkGH, fields.str());
  }
}

void Restrict(const cGH *cctkGH, int level, const std::vector<int> &groups,
              const char *const site) {
  DECLARE_CCTK_PARAMETERS;

  // AMR-C7.  Emitted BEFORE the assert below, deliberately: `do_restrict = no`
  // is a shipped parameter and this assert is live in both configurations
  // (`[P199]`: `NDEBUG` is defined in neither), so the run aborts here with no
  // Cactus-level message at all.  With the instrument on, the last line before
  // the abort names the site, the level and the parameter value.  This step
  // does not change that behaviour -- one mechanism per commit (R1).
  if (log_sched_on()) {
    std::ostringstream fields;
    fields << "site=" << site << "/restrict-enter level=" << level
           << " do_restrict=" << int(bool(do_restrict))
           << " ngroups=" << groups.size();
    log_sched(cctkGH, fields.str());
  }

#warning "TODO"
  assert(do_restrict);
  if (!do_restrict)
    return;

  static Timer timer("Restrict");
  Interval interval(timer);

  const int gi_regrid_error = CCTK_GroupIndex("CarpetXRegrid::regrid_error");
  assert(gi_regrid_error >= 0);

  // AMR-C7's counters.  `examined` is (patch, group) pairs this call looked at,
  // `restricted` is the ones on which an `average_down*` actually ran, and the
  // two skip counters say WHY the difference, so that a zero has one cause and
  // not three (`[N13]`).
  long n_patches = 0, n_examined = 0, n_restricted = 0;
  long n_skip_regrid_error = 0, n_skip_no_do_restrict = 0;

  for (const auto &patchdata : ghext->patchdata) {
    const int patch = patchdata.patch;
    if (level + 1 < int(patchdata.leveldata.size())) {
      ++n_patches;
      auto &leveldata = patchdata.leveldata.at(level);
      const auto &fineleveldata = patchdata.leveldata.at(level + 1);
      const active_levels_t active_levels(level, level + 1, patch, patch + 1);
      const active_levels_t active_fine_levels(level + 1, level + 2, patch,
                                               patch + 1);

      for (const int gi : groups) {
        cGroup group;
        int ierr = CCTK_GroupData(gi, &group);
        assert(!ierr);

        assert(group.grouptype == CCTK_GF);

        auto &groupdata = *leveldata.groupdata.at(gi);
        const auto &finegroupdata = *fineleveldata.groupdata.at(gi);
        const amrex::IntVect reffact{2, 2, 2};
        const nan_handling_t nan_handling = groupdata.do_checkpoint
                                                ? nan_handling_t::forbid_nans
                                                : nan_handling_t::allow_nans;

        ++n_examined;

        // Don't restrict the regridding error
        if (gi == gi_regrid_error) {
          ++n_skip_regrid_error;
          continue;
        }
        // Don't restrict groups that have restriction disabled
        if (!groupdata.do_restrict) {
          ++n_skip_no_do_restrict;
          continue;
        }
        ++n_restricted;

        // If there is more than one time level, then we don't restrict the
        // oldest.
        // TODO: during evolution, restrict only one time level
        int ntls = groupdata.mfab.size();
        int restrict_tl = ntls > 1 ? ntls - 1 : ntls;
        for (int tl = 0; tl < restrict_tl; ++tl) {

          for (int vi = 0; vi < groupdata.numvars; ++vi) {

            // Restriction only uses the interior
            error_if_invalid(finegroupdata, vi, tl, make_valid_int(), []() {
              return "Restrict on fine level before restricting";
            });
            poison_invalid_gf(active_fine_levels, gi, vi, tl);
            check_valid_gf(active_fine_levels, gi, vi, tl, nan_handling, []() {
              return "Restrict on fine level before restricting";
            });
            error_if_invalid(groupdata, vi, tl, make_valid_int(), []() {
              return "Restrict on coarse level before restricting";
            });
            poison_invalid_gf(active_levels, gi, vi, tl);
            check_valid_gf(active_levels, gi, vi, tl, nan_handling, []() {
              return "Restrict on coarse level before restricting";
            });
          }

#if 1
          {
            static Timer timer("Restrict::average_down");
            Interval interval(timer);
#warning                                                                       \
    "TODO: Allow different restriction operators, and ensure this is conservative"
            // rank: 0: vertex, 1: edge, 2: face, 3: volume
            int rank = 0;
            for (int d = 0; d < dim; ++d)
              rank += groupdata.indextype.at(d);
            switch (rank) {
            case 0:
              average_down_nodal(*finegroupdata.mfab.at(tl),
                                 *groupdata.mfab.at(tl), reffact);
              break;
            case 1:
              average_down_edges(*finegroupdata.mfab.at(tl),
                                 *groupdata.mfab.at(tl), reffact);
              break;
            case 2:
              average_down_faces(*finegroupdata.mfab.at(tl),
                                 *groupdata.mfab.at(tl), reffact);
              break;
            case 3:
              average_down(*finegroupdata.mfab.at(tl), *groupdata.mfab.at(tl),
                           0, groupdata.numvars, reffact);
              break;
            default:
              assert(0);
            }
          }
#endif

          // TODO: Also remember old why_valid for interior?
          for (int vi = 0; vi < groupdata.numvars; ++vi) {
            // Should we mark ghosts and maybe outer boundaries as
            // valid as well?
            groupdata.valid.at(tl).at(vi).set_invalid(
                make_valid_outer() | make_valid_ghosts(),
                []() { return "Restrict"; });
            poison_invalid_gf(active_levels, gi, vi, tl);
            check_valid_gf(active_levels, gi, vi, tl, nan_handling, []() {
              return "Restrict on coarse level after restricting";
            });
          }

        } // for tl
      } // for gi
    } // if level exists
  } // for patchdata

  if (log_sched_on()) {
    std::ostringstream fields;
    fields << "site=" << site << "/restrict-exit level=" << level
           << " patches=" << n_patches << " examined=" << n_examined
           << " restricted=" << n_restricted
           << " skip_regrid_error=" << n_skip_regrid_error
           << " skip_no_do_restrict=" << n_skip_no_do_restrict;
    log_sched(cctkGH, fields.str());
  }
}

void Restrict(const cGH *cctkGH, int level, const char *const site) {
  const int numgroups = CCTK_NumGroups();
  std::vector<int> groups;
  groups.reserve(numgroups);
  const auto &patchdata0 = ghext->patchdata.at(0);
  const auto &leveldata0 = patchdata0.leveldata.at(0);
  for (const auto &groupdataptr : leveldata0.groupdata) {
    // Restrict only grid functions
    if (groupdataptr) {
      auto &restrict groupdata = *groupdataptr;
      // Restrict only evolved grid functions
      if (groupdata.do_checkpoint)
        groups.push_back(groupdata.groupindex);
    }
  }
  Restrict(cctkGH, level, groups, site);
}

// storage handling
namespace {
int GroupStorageCrease(const cGH *cctkGH, int n_groups, const int *groups,
                       const int *requested_tls, int *status, const bool inc) {
  DECLARE_CCTK_PARAMETERS;

  assert(cctkGH);
  assert(n_groups >= 0);
  assert(groups);
  assert(requested_tls);
  for (int n = 0; n < n_groups; ++n) {
    if (groups[n] < 0 or groups[n] >= CCTK_NumGroups()) {
      CCTK_VWARN(CCTK_WARN_ALERT, "Group index %d is illegal", groups[n]);
      return -1;
    }
    assert(groups[n] >= 0 and groups[n] < CCTK_NumGroups());
    assert(requested_tls[n] >= 0 or requested_tls[n] == -1);
  }

  // sanitize list of requested timelevels
  std::vector<int> tls(n_groups);
  for (int n = 0; n < n_groups; ++n) {
    int ntls = requested_tls[n];
    int const declared_tls = CCTK_DeclaredTimeLevelsGI(groups[n]);
    if (inc and declared_tls < 2 and ntls > declared_tls) {
      CCTK_VWARN(CCTK_WARN_ALERT,
                 "Attempting to activate %d timelevels for group '%s' which "
                 "only has a single timelevel declared in interface.ccl. "
                 "Please declared at least 2 timelevels in interface.ccl to "
                 "allow more timelevels to be created at runtime.",
                 ntls, CCTK_FullGroupName(groups[n]));
      ntls = declared_tls;
    }
    if (ntls == -1) {
      ntls = declared_tls;
    }
    tls.at(n) = ntls;
  }

  // TODO: actually do something
  int min_num_timelevels = INT_MAX;
  for (int n = 0; n < n_groups; ++n) {
    int const gid = groups[n];

    cGroup group;
    int ierr = CCTK_GroupData(gid, &group);
    assert(not ierr);

    // Record previous number of allocated time levels
    if (status) {
      // Note: This remembers only the last level
      status[n] = group.numtimelevels;
    }

    // Record (minimum of) current number of time levels
    using std::min;
    min_num_timelevels = min(min_num_timelevels, group.numtimelevels);
  } // for n
  if (min_num_timelevels == INT_MAX) {
    min_num_timelevels = 0;
  }

  return min_num_timelevels;
}
} // namespace

int GroupStorageIncrease(const cGH *cctkGH, int n_groups, const int *groups,
                         const int *tls, int *status) {
  DECLARE_CCTK_PARAMETERS;

  return GroupStorageCrease(cctkGH, n_groups, groups, tls, status, true);
}

int GroupStorageDecrease(const cGH *cctkGH, int n_groups, const int *groups,
                         const int *tls, int *status) {
  DECLARE_CCTK_PARAMETERS;

  return GroupStorageCrease(cctkGH, n_groups, groups, tls, status, false);
}

int EnableGroupStorage(const cGH *cctkGH, const char *groupname) {
  const int group = CCTK_GroupIndex(groupname);
  assert(group >= 0 and group < CCTK_NumGroups());
  // TODO: decide whether to use CCTK_MaxActiveTimeLevelsGI
  const int tls = CCTK_DeclaredTimeLevelsGI(group);
  int status;
  GroupStorageIncrease(cctkGH, 1, &group, &tls, &status);
  // Return whether storage was allocated previously
  return status;
}

int DisableGroupStorage(const cGH *cctkGH, const char *groupname) {
  const int group = CCTK_GroupIndex(groupname);
  assert(group >= 0 and group < CCTK_NumGroups());
  const int tls = 0;
  int status;
  GroupStorageDecrease(cctkGH, 1, &group, &tls, &status);
  // Return whether storage was allocated previously
  return status;
}

} // namespace CarpetX
