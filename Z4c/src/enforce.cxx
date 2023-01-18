#include "physics.hxx"

#include <loop_device.hxx>
#include <mat.hxx>
#include <simd.hxx>
#include <vec.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#ifdef __CUDACC__
#include <nvToolsExt.h>
#endif

#include <cmath>
#include <sstream>

namespace Z4c {
using namespace Arith;
using namespace Loop;
using namespace std;

extern "C" void Z4c_Enforce(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_Z4c_Enforce;
  DECLARE_CCTK_PARAMETERS;

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout1(cctkGH, indextype);

  const GF3D2<CCTK_REAL> &gf_chi = chi;

  const GF3D2<CCTK_REAL> &gf_alphaG = alphaG;

  const smat<GF3D2<CCTK_REAL>, 3> gf_gammat{
      gammatxx, gammatxy, gammatxz, gammatyy, gammatyz, gammatzz,
  };

  const smat<GF3D2<CCTK_REAL>, 3> gf_At{
      Atxx, Atxy, Atxz, Atyy, Atyz, Atzz,
  };

  typedef simd<CCTK_REAL> vreal;
  typedef simdl<CCTK_REAL> vbool;
  constexpr size_t vsize = tuple_size_v<vreal>;

#ifdef __CUDACC__
  const nvtxRangeId_t range = nvtxRangeStartA("Z4c_Enforce::enforce");
#endif
  grid.loop_all_device<0, 0, 0, vsize>(
      grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
        const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
        const GF3D2index index1(layout1, p.I);

        // Load
        const vreal chi_old = gf_chi(mask, index1, 1);
        const vreal alphaG_old = gf_alphaG(mask, index1);

        const smat<vreal, 3> gammat_old =
            gf_gammat(mask, index1, one<smat<int, 3> >()());
        const smat<vreal, 3> At_old = gf_At(mask, index1);

        // Enforce floors

        const vreal chi = fmax(vreal(chi_floor), chi_old);
        const vreal alphaG = fmax(vreal(alphaG_floor - 1), alphaG_old);

        // Enforce algebraic constraints
        // See arXiv:1212.2901 [gr-qc].

        const vreal detgammat_old = calc_det(gammat_old);
        const vreal chi1_old = 1 / cbrt(detgammat_old);
        const smat<vreal, 3> gammat([&](int a, int b) ARITH_INLINE {
          return chi1_old * gammat_old(a, b);
        });
#ifdef CCTK_DEBUG
        const vreal detgammat = calc_det(gammat);
        const vreal gammat_norm = maxabs(gammat);
        const vreal gammat_scale = gammat_norm;
#ifndef __CUDACC__
        if (!(all(fabs(detgammat - 1) <= 1.0e-12 * gammat_scale))) {
          ostringstream buf;
          buf << "det gammat is not one: gammat=" << gammat
              << " det(gammat)=" << detgammat;
          CCTK_VERROR("%s", buf.str().c_str());
        }
#endif
        assert(all(fabs(detgammat - 1) <= 1.0e-12 * gammat_scale));
#endif

        const smat<vreal, 3> gammatu = calc_inv(gammat, vreal(1));

        const vreal traceAt_old = sum_symm<3>([&](int x, int y) ARITH_INLINE {
          return gammatu(x, y) * At_old(x, y);
        });
        const smat<vreal, 3> At([&](int a, int b) ARITH_INLINE {
          return At_old(a, b) - traceAt_old / 3 * gammat(a, b);
        });
#ifdef CCTK_DEBUG
        const vreal traceAt = sum_symm<3>([&](int x, int y) ARITH_INLINE {
          return gammatu(x, y) * At(x, y);
        });
        const vreal gammatu_norm = maxabs(gammatu);
        const vreal At_norm = maxabs(At);
        const vreal At_scale = fmax(fmax(gammat_norm, gammatu_norm), At_norm);
#ifndef __CUDACC__
        if (!(all(fabs(traceAt) <= 1.0e-12 * At_scale))) {
          ostringstream buf;
          buf << "tr At: At=" << At << " tr(At)=" << traceAt;
          CCTK_VERROR("%s", buf.str().c_str());
        }
#endif
        assert(all(fabs(traceAt) <= 1.0e-12 * At_scale));
#endif

        // Store
        gf_chi.store(mask, index1, chi);
        gf_gammat.store(mask, index1, gammat);
        gf_At.store(mask, index1, At);
        gf_alphaG.store(mask, index1, alphaG);
      });
#ifdef __CUDACC__
  nvtxRangeEnd(range);
#endif
}

} // namespace Z4c
