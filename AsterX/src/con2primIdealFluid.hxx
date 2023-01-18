#include <con2primFactory.hxx>
#ifndef CON2PRIM_IDEALFLUID_HXX
#define CON2PRIM_IDEALFLUID_HXX
namespace AsterX
{

    class idealFluid : public con2primFactory
    {
    public:
        /* Some attributes */
        CCTK_REAL GammaIdealFluid;
        /* Constructor */
        CCTK_HOST CCTK_DEVICE idealFluid(CCTK_REAL gamma, CCTK_REAL (&cons)[NCONS], CCTK_REAL (&prim)[NPRIMS], CCTK_REAL (&gcov)[4][4], CCTK_REAL (&gcon)[4][4]);
        CCTK_HOST CCTK_DEVICE void get_Ssq_Exact();  // From cons (exact)
        CCTK_HOST CCTK_DEVICE void get_Bsq_Exact();  // From cons (exact)
        CCTK_HOST CCTK_DEVICE void get_BiSi_Exact(); // From cons (exact)
        CCTK_HOST CCTK_DEVICE void get_WLorentz_bsq_Seeds();

        /* Called by 2DNRNoble */
        CCTK_HOST CCTK_DEVICE void get_Press_Seed();
        CCTK_HOST CCTK_DEVICE void get_Z_Seed();
        CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_f0(CCTK_REAL Z, CCTK_REAL Vsq);
        CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_f1(CCTK_REAL Z, CCTK_REAL Vsq);
        CCTK_HOST CCTK_DEVICE CCTK_REAL get_Press_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq);
        CCTK_HOST CCTK_DEVICE CCTK_REAL get_dPdZ_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq);
        CCTK_HOST CCTK_DEVICE CCTK_REAL get_dPdVsq_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq);
        CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_df0dZ(CCTK_REAL Z, CCTK_REAL Vsq);
        CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_df0dVsq(CCTK_REAL Z, CCTK_REAL Vsq);
        CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_df1dZ(CCTK_REAL Z, CCTK_REAL Vsq);
        CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_df1dVsq(CCTK_REAL Z, CCTK_REAL Vsq);
        CCTK_HOST CCTK_DEVICE void WZ2Prim();

        /* Called by 1DBrentPalenzuela */
        CCTK_INT Failed_1DBrentPalenzuela;
        CCTK_REAL xPalenzuela_Sol;
        CCTK_HOST CCTK_DEVICE CCTK_REAL get_Press_funcRhoEps(CCTK_REAL &rho_loc, CCTK_REAL &eps_loc);
        CCTK_HOST CCTK_DEVICE void xPalenzuelaToPrim();

        /* Destructor */
        CCTK_HOST CCTK_DEVICE ~idealFluid();
    };

    /* Constructor */
    CCTK_HOST CCTK_DEVICE inline idealFluid::idealFluid(CCTK_REAL gamma,
                                                        CCTK_REAL (&cons)[NCONS],
                                                        CCTK_REAL (&prim)[NPRIMS],
                                                        CCTK_REAL (&g_lo)[4][4],
                                                        CCTK_REAL (&g_up)[4][4])
    {

        GammaIdealFluid = gamma;

        PrimitiveVarsSeed[RHO] = prim[RHO];
        PrimitiveVarsSeed[V1_CON] = prim[V1_CON];
        PrimitiveVarsSeed[V2_CON] = prim[V2_CON];
        PrimitiveVarsSeed[V3_CON] = prim[V3_CON];
        PrimitiveVarsSeed[EPS] = prim[EPS];
        PrimitiveVarsSeed[B1] = prim[B1];
        PrimitiveVarsSeed[B2] = prim[B2];
        PrimitiveVarsSeed[B3] = prim[B3];

        ConservedVars[CONS_D] = cons[CONS_D];
        ConservedVars[S1_COV] = cons[S1_COV];
        ConservedVars[S2_COV] = cons[S2_COV];
        ConservedVars[S3_COV] = cons[S3_COV];
        ConservedVars[EPS] = cons[EPS];
        ConservedVars[B1] = cons[B1];
        ConservedVars[B2] = cons[B2];
        ConservedVars[B3] = cons[B3];

        gcov[TT] = g_lo[0][0];
        gcov[TX] = g_lo[0][1];
        gcov[TY] = g_lo[0][2];
        gcov[TZ] = g_lo[0][3];
        gcov[XX] = g_lo[1][1];
        gcov[XY] = g_lo[1][2];
        gcov[XZ] = g_lo[1][3];
        gcov[YY] = g_lo[2][2];
        gcov[YZ] = g_lo[2][3];
        gcov[ZZ] = g_lo[3][3];

        gcon[TT] = g_up[0][0];
        gcon[TX] = g_up[0][1];
        gcon[TY] = g_up[0][2];
        gcon[TZ] = g_up[0][3];
        gcon[XX] = g_up[1][1];
        gcon[XY] = g_up[1][2];
        gcon[XZ] = g_up[1][3];
        gcon[YY] = g_up[2][2];
        gcon[YZ] = g_up[2][3];
        gcon[ZZ] = g_up[3][3];

        get_Ssq_Exact();
        get_Bsq_Exact();
        get_BiSi_Exact();
        get_WLorentz_bsq_Seeds();

        /* update rho seed from CONS_D and gamma */
        // rho consistent with con[CONS_D] should be better guess than rho from last
        // timestep
        PrimitiveVarsSeed[RHO] = ConservedVars[CONS_D] / WLorentz_Seed;
    }

    CCTK_HOST CCTK_DEVICE inline void idealFluid::get_Ssq_Exact()
    {

        /* calculate S_squared */
        Ssq = ConservedVars[S1_COV] *
              (gcon[XX] * ConservedVars[S1_COV] + gcon[XY] * ConservedVars[S2_COV] +
               gcon[XZ] * ConservedVars[S3_COV]);
        Ssq += ConservedVars[S2_COV] *
               (gcon[XY] * ConservedVars[S1_COV] + gcon[YY] * ConservedVars[S2_COV] +
                gcon[YZ] * ConservedVars[S3_COV]);
        Ssq += ConservedVars[S3_COV] *
               (gcon[XZ] * ConservedVars[S1_COV] + gcon[YZ] * ConservedVars[S2_COV] +
                gcon[ZZ] * ConservedVars[S3_COV]);
        if ((Ssq < 0.) && (fabs(Ssq) < 1.0e-13))
        {
            Ssq = fabs(Ssq);
        }
    }

    CCTK_HOST CCTK_DEVICE inline void idealFluid::get_Bsq_Exact()
    {
        CCTK_REAL B1_cov = gcov[XX] * PrimitiveVarsSeed[B1] + gcov[XY] * PrimitiveVarsSeed[B2] + gcov[XZ] * PrimitiveVarsSeed[B3];
        CCTK_REAL B2_cov = gcov[XY] * PrimitiveVarsSeed[B1] + gcov[YY] * PrimitiveVarsSeed[B2] + gcov[YZ] * PrimitiveVarsSeed[B3];
        CCTK_REAL B3_cov = gcov[XZ] * PrimitiveVarsSeed[B1] + gcov[YZ] * PrimitiveVarsSeed[B2] + gcov[ZZ] * PrimitiveVarsSeed[B3];
        Bsq = B1_cov * PrimitiveVarsSeed[B1] + B2_cov * PrimitiveVarsSeed[B2] + B3_cov * PrimitiveVarsSeed[B3];
    }

    CCTK_HOST CCTK_DEVICE inline void idealFluid::get_BiSi_Exact()
    {
        BiSi = PrimitiveVarsSeed[B1] * ConservedVars[S1_COV] + PrimitiveVarsSeed[B2] * ConservedVars[S2_COV] + PrimitiveVarsSeed[B3] * ConservedVars[S3_COV];
    }

    CCTK_HOST CCTK_DEVICE inline void idealFluid::get_WLorentz_bsq_Seeds()
    {
        // covariant Valencia velocity:
        CCTK_REAL v1_cov = gcov[XX] * PrimitiveVarsSeed[V1_CON] +
                           gcov[XY] * PrimitiveVarsSeed[V2_CON] +
                           gcov[XZ] * PrimitiveVarsSeed[V3_CON];
        CCTK_REAL v2_cov = gcov[XY] * PrimitiveVarsSeed[V1_CON] +
                           gcov[YY] * PrimitiveVarsSeed[V2_CON] +
                           gcov[YZ] * PrimitiveVarsSeed[V3_CON];
        CCTK_REAL v3_cov = gcov[XZ] * PrimitiveVarsSeed[V1_CON] +
                           gcov[YZ] * PrimitiveVarsSeed[V2_CON] +
                           gcov[ZZ] * PrimitiveVarsSeed[V3_CON];

        CCTK_REAL vsq = v1_cov * PrimitiveVarsSeed[V1_CON] +
                        v2_cov * PrimitiveVarsSeed[V2_CON] +
                        v3_cov * PrimitiveVarsSeed[V3_CON];

        if ((vsq < 0.) && (fabs(vsq) < 1.0e-13))
        {
            vsq = fabs(vsq);
        }

        WLorentz_Seed = 1. / sqrt(1. - vsq);

        CCTK_REAL alp = sqrt(-1. / gcon[TT]);
        CCTK_REAL bt = WLorentz_Seed / alp *
                       (ConservedVars[B1] * v1_cov + ConservedVars[B2] * v2_cov +
                        ConservedVars[B3] * v3_cov);

        bsq_Seed = (Bsq + (alp * bt) * (alp * bt)) / (WLorentz_Seed * WLorentz_Seed);
    }

    /* Called by 2dNRNoble */
    CCTK_HOST CCTK_DEVICE inline void idealFluid::get_Press_Seed()
    {
        Press_Seed =
            PrimitiveVarsSeed[RHO] * PrimitiveVarsSeed[EPS] * (GammaIdealFluid - 1.0);
    }

    CCTK_HOST CCTK_DEVICE inline void idealFluid::get_Z_Seed()
    {
        Z_Seed = (PrimitiveVarsSeed[RHO] +
                  PrimitiveVarsSeed[EPS] * PrimitiveVarsSeed[RHO] + Press_Seed) *
                 WLorentz_Seed * WLorentz_Seed;
    }

    CCTK_HOST CCTK_DEVICE inline CCTK_REAL idealFluid::get_2DNRNoble_f0(CCTK_REAL Z,
                                                                        CCTK_REAL Vsq)
    {
        return (Vsq * (Bsq + Z) * (Bsq + Z) -
                (BiSi * BiSi * (Bsq + 2.0 * Z)) / (Z * Z) - Ssq);
    }

    CCTK_HOST CCTK_DEVICE inline CCTK_REAL idealFluid::get_2DNRNoble_f1(CCTK_REAL Z,
                                                                        CCTK_REAL Vsq)
    {
        CCTK_REAL Press = get_Press_funcZVsq(Z, Vsq);
        return ConservedVars[TAU] + ConservedVars[CONS_D] - Bsq / 2.0 * (1 + Vsq) +
               BiSi * BiSi / 2.0 / (Z * Z) - Z + Press;
    }

    CCTK_HOST CCTK_DEVICE inline CCTK_REAL idealFluid::get_2DNRNoble_df0dZ(CCTK_REAL Z,
                                                                           CCTK_REAL Vsq)
    {
        return (2.0 * Vsq * (Bsq + Z) - 2.0 * BiSi * BiSi / (Z * Z) +
                2.0 * BiSi * BiSi * (Bsq + 2.0 * Z) / (Z * Z * Z));
    }

    CCTK_HOST CCTK_DEVICE inline CCTK_REAL
    idealFluid::get_2DNRNoble_df0dVsq(CCTK_REAL Z, CCTK_REAL Vsq)
    {
        return (Bsq + Z) * (Bsq + Z);
    }

    CCTK_HOST CCTK_DEVICE inline CCTK_REAL idealFluid::get_2DNRNoble_df1dZ(CCTK_REAL Z,
                                                                           CCTK_REAL Vsq)
    {
        return (-BiSi * BiSi / (Z * Z * Z) - 1.0 + get_dPdZ_funcZVsq(Z, Vsq));
    }

    CCTK_HOST CCTK_DEVICE inline CCTK_REAL
    idealFluid::get_2DNRNoble_df1dVsq(CCTK_REAL Z, CCTK_REAL Vsq)
    {
        return (-Bsq / 2.0 + get_dPdVsq_funcZVsq(Z, Vsq));
    }

    CCTK_HOST CCTK_DEVICE inline CCTK_REAL idealFluid::get_Press_funcZVsq(CCTK_REAL Z,
                                                                          CCTK_REAL Vsq)
    {
        return ((Z * (1.0 - Vsq) - ConservedVars[CONS_D] * sqrt(1.0 - Vsq)) *
                (GammaIdealFluid - 1.0) / (GammaIdealFluid));
    }

    CCTK_HOST CCTK_DEVICE inline CCTK_REAL idealFluid::get_dPdZ_funcZVsq(CCTK_REAL Z,
                                                                         CCTK_REAL Vsq)
    {
        return ((1.0 - Vsq) * (GammaIdealFluid - 1.0) / GammaIdealFluid);
    }

    CCTK_HOST CCTK_DEVICE inline CCTK_REAL idealFluid::get_dPdVsq_funcZVsq(CCTK_REAL Z,
                                                                           CCTK_REAL Vsq)
    {
        return ((-Z + ConservedVars[CONS_D] / (2.0 * sqrt(1.0 - Vsq))) *
                (GammaIdealFluid - 1.0) / GammaIdealFluid);
    }

    CCTK_HOST CCTK_DEVICE inline void idealFluid::WZ2Prim()
    {
        CCTK_REAL W_Sol = 1.0 / sqrt(1.0 - vsq_Sol);

        PrimitiveVars[RHO] = ConservedVars[CONS_D] / W_Sol;

        PrimitiveVars[V1_CON] =
            (gcon[XX] * ConservedVars[S1_COV] + gcon[XY] * ConservedVars[S2_COV] +
             gcon[XZ] * ConservedVars[S3_COV]) /
            (Z_Sol + Bsq);
        PrimitiveVars[V1_CON] += BiSi * ConservedVars[B1] / (Z_Sol * (Z_Sol + Bsq));

        PrimitiveVars[V2_CON] =
            (gcon[XY] * ConservedVars[S1_COV] + gcon[YY] * ConservedVars[S2_COV] +
             gcon[YZ] * ConservedVars[S3_COV]) /
            (Z_Sol + Bsq);
        PrimitiveVars[V2_CON] += BiSi * ConservedVars[B2] / (Z_Sol * (Z_Sol + Bsq));

        PrimitiveVars[V3_CON] =
            (gcon[XZ] * ConservedVars[S1_COV] + gcon[YZ] * ConservedVars[S2_COV] +
             gcon[ZZ] * ConservedVars[S3_COV]) /
            (Z_Sol + Bsq);
        PrimitiveVars[V3_CON] += BiSi * ConservedVars[B3] / (Z_Sol * (Z_Sol + Bsq));

        PrimitiveVars[EPS] =
            (Z_Sol * (1. - vsq_Sol) / PrimitiveVars[RHO] - 1.0) / GammaIdealFluid;
    }

    /* Called by 1DBrentPalenzuela */

    CCTK_HOST CCTK_DEVICE inline CCTK_REAL idealFluid::get_Press_funcRhoEps(CCTK_REAL &rho_loc, CCTK_REAL &eps_loc)
    {
        return rho_loc * eps_loc * (GammaIdealFluid - 1.0);
    }

    CCTK_HOST CCTK_DEVICE inline void idealFluid::xPalenzuelaToPrim()
    {

        const CCTK_REAL qPalenzuela = ConservedVars[TAU] / ConservedVars[CONS_D];
        const CCTK_REAL rPalenzuela = Ssq / pow(ConservedVars[CONS_D], 2);
        const CCTK_REAL sPalenzuela = Bsq / ConservedVars[CONS_D];
        const CCTK_REAL tPalenzuela = BiSi / pow(ConservedVars[CONS_D], 3. / 2.);

        // (i)
        CCTK_REAL Wminus2 = 1.0 - (xPalenzuela_Sol * xPalenzuela_Sol * rPalenzuela + (2.0 * xPalenzuela_Sol + sPalenzuela) * tPalenzuela * tPalenzuela) / (xPalenzuela_Sol * xPalenzuela_Sol * (xPalenzuela_Sol + sPalenzuela) * (xPalenzuela_Sol + sPalenzuela));
        Wminus2 = fmin(fmax(Wminus2, 1e-10), 1 - 1e-10);
        const CCTK_REAL W_loc = pow(Wminus2, -0.5);

        // (ii)
        PrimitiveVars[RHO] = ConservedVars[CONS_D] / W_loc;

        // (iii)
        PrimitiveVars[EPS] = W_loc - 1.0 + (1.0 - W_loc * W_loc) * xPalenzuela_Sol / W_loc + W_loc * (qPalenzuela - sPalenzuela + tPalenzuela * tPalenzuela / (2 * xPalenzuela_Sol * xPalenzuela_Sol) + sPalenzuela / (2.0 * W_loc * W_loc));

        // (iv)
        //CCTK_REAL P_loc = get_Press_funcRhoEps(rho_loc, eps_loc);

        // Taken from WZ2Prim (2DNRNoble)
        Z_Sol = xPalenzuela_Sol * PrimitiveVars[RHO] * W_loc;
        CCTK_REAL alp = sqrt(-1. / gcon[TT]);

        CCTK_REAL v1_Valencia = (gcon[XX] * ConservedVars[S1_COV] + gcon[XY] * ConservedVars[S2_COV] + gcon[XZ] * ConservedVars[S3_COV]) / (Z_Sol + Bsq);
        v1_Valencia += BiSi * PrimitiveVars[B1] / (Z_Sol * (Z_Sol + Bsq));

        CCTK_REAL v2_Valencia = (gcon[XY] * ConservedVars[S1_COV] + gcon[YY] * ConservedVars[S2_COV] + gcon[YZ] * ConservedVars[S3_COV]) / (Z_Sol + Bsq);
        v2_Valencia += BiSi * PrimitiveVars[B2] / (Z_Sol * (Z_Sol + Bsq));

        CCTK_REAL v3_Valencia = (gcon[XZ] * ConservedVars[S1_COV] + gcon[YZ] * ConservedVars[S2_COV] + gcon[ZZ] * ConservedVars[S3_COV]) / (Z_Sol + Bsq);
        v3_Valencia += BiSi * PrimitiveVars[B3] / (Z_Sol * (Z_Sol + Bsq));

        PrimitiveVars[V1_CON] = alp * v1_Valencia - alp * alp * gcon[TX];
        PrimitiveVars[V2_CON] = alp * v2_Valencia - alp * alp * gcon[TY];
        PrimitiveVars[V3_CON] = alp * v3_Valencia - alp * alp * gcon[TZ];

        PrimitiveVars[B1] = ConservedVars[B1];
        PrimitiveVars[B2] = ConservedVars[B2];
        PrimitiveVars[B3] = ConservedVars[B3];
    }

    /* Destructor */
    CCTK_HOST CCTK_DEVICE inline idealFluid::~idealFluid()
    {
        // How to destruct properly a vector?
    }
}

#endif // #ifndef CON2PRIM_HXX
