#ifndef CARPETX_CARPETX_IO_SLICE_HXX
#define CARPETX_CARPETX_IO_SLICE_HXX

#include <cctk.h>

#include <vect.hxx>

#include <array>
#include <cmath>
#include <optional>
#include <utility>

namespace CarpetX {

// A 2D planar slice: one normal axis pinned to a physical coordinate.
// "Exactly one normal axis" is structural, so 0D/1D/3D "slices" are
// unrepresentable.
struct slice_t {
  int normal_dir;  // 0=x, 1=y, 2=z
  CCTK_REAL coord; // physical coordinate of the plane along normal_dir

  // The two in-plane axes (d != normal_dir), ascending. Used for the
  // output_directions metadata: xy->{0,1}, xz->{0,2}, yz->{1,2}.
  std::array<int, 2> inplane_dirs() const {
    std::array<int, 2> r{};
    int k = 0;
    for (int d = 0; d < 3; ++d)
      if (d != normal_dir)
        r[k++] = d;
    return r;
  }

  // Restrict a half-open index box [lo, hi) to a thickness-1 slab at the
  // plane. Callers pass a centering-aware origin x0[d] (the physical coordinate
  // of index 0; ProbLo(d) + cellCentered(d)*dx[d]/2, matching io_tsv.cxx:274),
  // so idx = lrint((coord - x0)/dx) selects the correct plane for both vertex-
  // and cell-centred data. Returns [lo, hi) with lo[n]=idx, hi[n]=idx+1, or
  // nullopt when the plane lies outside.
  std::optional<std::pair<Arith::vect<int, 3>, Arith::vect<int, 3> > >
  restrict_box(Arith::vect<int, 3> lo, Arith::vect<int, 3> hi,
               const Arith::vect<CCTK_REAL, 3> &x0,
               const Arith::vect<CCTK_REAL, 3> &dx) const {
    const int n = normal_dir;
    const int idx = std::lrint((coord - x0[n]) / dx[n]);
    if (idx < lo[n] || idx >= hi[n])
      return std::nullopt;
    lo[n] = idx;
    hi[n] = idx + 1;
    return std::make_pair(lo, hi);
  }

  // Membership uses a centering-independent cell-canonical frame (origin cx0,
  // valid box [cval_lo, cval_hi), half-open so a boundary plane belongs to one
  // box); the slab uses the centering-aware index (origin sx0) so it matches
  // the variable's own grid, keeping every multivar aligned with the multimesh.
  // Returns nullopt when the plane lies outside this box.
  std::optional<std::pair<Arith::vect<int, 3>, Arith::vect<int, 3> > >
  restrict_box_interior(Arith::vect<int, 3> ext_lo, Arith::vect<int, 3> ext_hi,
                        const Arith::vect<int, 3> &cval_lo,
                        const Arith::vect<int, 3> &cval_hi,
                        const Arith::vect<CCTK_REAL, 3> &cx0,
                        const Arith::vect<CCTK_REAL, 3> &sx0,
                        const Arith::vect<CCTK_REAL, 3> &dx) const {
    const int n = normal_dir;
    const int idx_cell = std::lrint((coord - cx0[n]) / dx[n]);
    if (idx_cell < cval_lo[n] || idx_cell >= cval_hi[n])
      return std::nullopt;
    const int idx_var = std::lrint((coord - sx0[n]) / dx[n]);
    ext_lo[n] = idx_var;
    ext_hi[n] = idx_var + 1;
    return std::make_pair(ext_lo, ext_hi);
  }
};

// Copy the sub-box [sub_lo, sub_hi) out of a Fortran-ordered (stride-1 in
// dir 0) source array spanning [src_lo, src_hi) into a contiguous
// Fortran-ordered destination buffer of shape (sub_hi - sub_lo). `dst` must
// hold at least prod(sub_hi - sub_lo) elements. Lifted verbatim from the
// openPMD ghost-strip loop (io_openpmd.cxx:1836-1860).
inline void extract_subbox(CCTK_REAL *restrict dst,
                           const CCTK_REAL *restrict src,
                           const Arith::vect<int, 3> &src_lo,
                           const Arith::vect<int, 3> &src_hi,
                           const Arith::vect<int, 3> &sub_lo,
                           const Arith::vect<int, 3> &sub_hi) {
  const Arith::vect<int, 3> src_shape = src_hi - src_lo;
  const Arith::vect<int, 3> sub_shape = sub_hi - sub_lo;
  const Arith::vect<int, 3> offset = sub_lo - src_lo;
  constexpr int sdi = 1;
  const int sdj = sdi * src_shape[0];
  const int sdk = sdj * src_shape[1];
  const CCTK_REAL *restrict const src0 =
      src + sdi * offset[0] + sdj * offset[1] + sdk * offset[2];
  constexpr int ddi = 1;
  const int ddj = ddi * sub_shape[0];
  const int ddk = ddj * sub_shape[1];
  for (int k = 0; k < sub_shape[2]; ++k)
    for (int j = 0; j < sub_shape[1]; ++j)
#pragma omp simd
      for (int i = 0; i < sub_shape[0]; ++i)
        dst[ddi * i + ddj * j + ddk * k] = src0[sdi * i + sdj * j + sdk * k];
}

} // namespace CarpetX

#endif // #ifndef CARPETX_CARPETX_IO_SLICE_HXX
