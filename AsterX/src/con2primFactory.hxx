#ifndef CON2PRIMFACTORY_HXX
#define CON2PRIMFACTORY_HXX
#include "AMReX_GpuQualifiers.H"
#define CCTK_DEVICE AMREX_GPU_DEVICE
#define CCTK_HOST AMREX_GPU_HOST
#include "utils.hxx"

constexpr int CONS_D = 0;
constexpr int S1_COV = 1;
constexpr int S2_COV = 2;
constexpr int S3_COV = 3;
constexpr int TAU = 4;
constexpr int B1 = 5;
constexpr int B2 = 6;
constexpr int B3 = 7;

constexpr int RHO = 0;
constexpr int V1_CON = 1;
constexpr int V2_CON = 2;
constexpr int V3_CON = 3;
constexpr int V1_COV = 1;
constexpr int V2_COV = 2;
constexpr int V3_COV = 3;
constexpr int EPS = 4;
constexpr int PRESS = 5;

constexpr int TT = 0;
constexpr int TX = 1;
constexpr int TY = 2;
constexpr int TZ = 3;
constexpr int XX = 4;
constexpr int XY = 5;
constexpr int XZ = 6;
constexpr int YY = 7;
constexpr int YZ = 8;
constexpr int ZZ = 9;

#include <math.h>

namespace AsterX
{

  constexpr int NCONS = 8;
  constexpr int NPRIMS = 8;

  /* Abstract class con2primFactory */
  class con2primFactory
  {
  public:
    /* The constructor must initialize the following variables */
    CCTK_REAL ConservedVars[NCONS];      // Conserved to solve
    CCTK_REAL PrimitiveVarsSeed[NPRIMS]; // Primitive seeds
    CCTK_REAL PrimitiveVars[NPRIMS];     // Primitive solution
    CCTK_INT maxIterations;
    CCTK_REAL tolerance;
    CCTK_REAL gcov[10], gcon[10];
    CCTK_REAL Bsq, BiSi, Ssq;
    CCTK_REAL WLorentz_Seed, bsq_Seed;
    CCTK_HOST CCTK_DEVICE void get_Ssq_Exact();  // From cons (exact)
    CCTK_HOST CCTK_DEVICE void get_Bsq_Exact();  // From cons (exact)
    CCTK_HOST CCTK_DEVICE void get_BiSi_Exact(); // From cons (exact)
    CCTK_HOST CCTK_DEVICE void get_WLorentz_bsq_Seeds();

    /* These must be set for 2DNRNoble scheme */
    CCTK_INT Failed_2DNRNoble;
    CCTK_INT Nit_2DNRNoble;
    CCTK_REAL vsq_Sol, Press_Seed, Z_Seed, Z_Sol, vsq_Seed;
    CCTK_HOST CCTK_DEVICE void get_Press_Seed(); // From seed prims and cons
    CCTK_HOST CCTK_DEVICE void get_Z_Seed();     // From seed prims and cons
    CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_f0(CCTK_REAL Z, CCTK_REAL Vsq);
    CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_f1(CCTK_REAL Z, CCTK_REAL Vsq);
    CCTK_HOST CCTK_DEVICE CCTK_REAL get_Press_funcZVsq(CCTK_REAL Z,
                                                       CCTK_REAL Vsq);
    CCTK_HOST CCTK_DEVICE CCTK_REAL get_dPdZ_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq);
    CCTK_HOST CCTK_DEVICE CCTK_REAL get_dPdVsq_funcZVsq(CCTK_REAL Z,
                                                        CCTK_REAL Vsq);
    CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_df0dZ(CCTK_REAL Z,
                                                        CCTK_REAL Vsq);
    CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_df0dVsq(CCTK_REAL Z,
                                                          CCTK_REAL Vsq);
    CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_df1dZ(CCTK_REAL Z,
                                                        CCTK_REAL Vsq);
    CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_df1dVsq(CCTK_REAL Z,
                                                          CCTK_REAL Vsq);
    CCTK_HOST CCTK_DEVICE void WZ2Prim();

    /* Called by 1DBrentPalenzuela */
    CCTK_INT Failed_1DBrentPalenzuela;
    CCTK_REAL xPalenzuela_Sol;
    CCTK_HOST CCTK_DEVICE CCTK_REAL get_Press_funcRhoEps(CCTK_REAL &rho_loc, CCTK_REAL &eps_loc);
    CCTK_HOST CCTK_DEVICE void xPalenzuelaToPrim();
  };

  template <typename typeEoS>
  CCTK_HOST CCTK_DEVICE void Con2Prim_2DNRNoble(CCTK_INT max_iter, CCTK_REAL tolf,
                                                typeEoS &plasma);

  template <typename T>
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  set_to_atmosphere(const T &rho_abs_min, const T &poly_K, const T &gamma,
                    const T sqrtg, const T Bsq, const vec<GF3D2<T>, 6> &prims,
                    const vec<GF3D2<T>, 5> &cons, const PointDesc &p)
  {

    prims(RHO)(p.I) = rho_abs_min;
    prims(V1_CON)(p.I) = 0.0;
    prims(V2_CON)(p.I) = 0.0;
    prims(V3_CON)(p.I) = 0.0;
    prims(PRESS)(p.I) = poly_K * pow(prims(RHO)(p.I), gamma);
    prims(EPS)(p.I) = prims(PRESS)(p.I) / ((gamma - 1.0) * prims(RHO)(p.I));

    cons(CONS_D)(p.I) = sqrtg * prims(RHO)(p.I);
    cons(S1_COV)(p.I) = 0.0;
    cons(S2_COV)(p.I) = 0.0;
    cons(S3_COV)(p.I) = 0.0;
    cons(TAU)(p.I) =
        sqrtg * ((prims(RHO)(p.I) * (1 + prims(EPS)(p.I)) + prims(PRESS)(p.I)) +
                 Bsq - (prims(PRESS)(p.I) + 0.5 * Bsq)) -
        cons(CONS_D)(p.I);
  }

} // namespace AsterX
#endif // #ifndef CON2PRIMFACTORY_HXX
