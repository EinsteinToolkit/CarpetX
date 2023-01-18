#ifndef WEYL_HXX
#define WEYL_HXX

#include "derivs.hxx"
#include "physics.hxx"

#include <cplx.hxx>
#include <defs.hxx>
#include <mat.hxx>
#include <rten.hxx>
#include <simd.hxx>
#include <vec.hxx>

#include <cmath>

namespace Weyl {

template <typename T> struct weyl_vars_noderivs {

  // Position
  const vec<T, 4> coord;

  // ADM variables
  const smat<T, 3> gamma;
  const T alpha;
  const vec<T, 3> beta;

  // Intermediate quantities
  const vec<T, 3> betal;
  const T gtt;

  // 4-metric
  const smat<T, 4> g;

  // Inverse 4-metric
  const T detg;
  const smat<T, 4> gu;

  // Tetrad
  const vec<T, 4> et, ephi, etheta, er;
  const vec<T, 4> l, n;
  const vec<cplx<T>, 4> m;

  inline ARITH_INLINE ARITH_DEVICE ARITH_HOST weyl_vars_noderivs(
      const T time, const vec<T, 3> &coord3, const smat<T, 3> &gamma,
      const T &alpha, const vec<T, 3> &beta)
      : coord{time, coord3(0), coord3(1), coord3(2)},
        //
        gamma(gamma), alpha(alpha), beta(beta),
        //
        betal([&](int a) ARITH_INLINE {
          return sum<3>([&](int x)
                            ARITH_INLINE { return gamma(a, x) * beta(x); });
        }),
        //
        gtt(-pow2(alpha) //
            + sum<3>([&](int x) ARITH_INLINE { return betal(x) * beta(x); })),
        //
        g([&](int a, int b) ARITH_INLINE {
          if (a == 0 && b == 0)
            return gtt;
          if (a == 0)
            return betal(b - 1);
          if (b == 0)
            return betal(a - 1);
          return gamma(a - 1, b - 1);
        }),
        detg(calc_det(g)),                   //
        gu(calc_inv(g, detg)),               //
        et(calc_et(gu)),                     //
        ephi(calc_ephi(coord, g)),           //
        etheta(calc_etheta(coord, g, ephi)), //
        er(calc_er(coord, g, etheta, ephi)), //
        l([&](int a) ARITH_INLINE { return (et(a) + er(a)) / sqrt(T(2)); }),
        n([&](int a) ARITH_INLINE { return (et(a) - er(a)) / sqrt(T(2)); }),
        m([&](int a)
              ARITH_INLINE { return cplx<T>(etheta(a), ephi(a)) / sqrt(T(2)); })
  //
  {}

  inline ARITH_INLINE ARITH_DEVICE ARITH_HOST
  weyl_vars_noderivs(const vec<T, 4> &coord, const smat<T, 4> &g)
      : coord(coord),                        //
        gamma(), alpha(), beta(), betal(),   //
        g(g), detg(calc_det(g)),             //
        gu(calc_inv(g, detg)),               //
        et(calc_et(gu)),                     //
        ephi(calc_ephi(coord, g)),           //
        etheta(calc_etheta(coord, g, ephi)), //
        er(calc_er(coord, g, etheta, ephi)), //
        l([&](int a) ARITH_INLINE { return (et(a) + er(a)) / sqrt(T(2)); }),
        n([&](int a) ARITH_INLINE { return (et(a) - er(a)) / sqrt(T(2)); }),
        m([&](int a)
              ARITH_INLINE { return cplx<T>(etheta(a), ephi(a)) / sqrt(T(2)); })
  //
  {}
};

template <typename T> struct weyl_vars : weyl_vars_noderivs<T> {

  // C++ is tedious:

  // Position
  using weyl_vars_noderivs<T>::coord;

  // ADM variables
  using weyl_vars_noderivs<T>::alpha;
  using weyl_vars_noderivs<T>::beta;
  using weyl_vars_noderivs<T>::gamma;

  // Intermediate quantities
  using weyl_vars_noderivs<T>::betal;

  // 4-metric
  using weyl_vars_noderivs<T>::g;

  // Inverse 4-metric
  using weyl_vars_noderivs<T>::detg;
  using weyl_vars_noderivs<T>::gu;

  // Tetrad
  using weyl_vars_noderivs<T>::et;
  using weyl_vars_noderivs<T>::ephi;
  using weyl_vars_noderivs<T>::etheta;
  using weyl_vars_noderivs<T>::er;
  using weyl_vars_noderivs<T>::l;
  using weyl_vars_noderivs<T>::n;
  using weyl_vars_noderivs<T>::m;

  // Time derivatives of ADM variables
  const smat<T, 3> k;
  const T dtalpha;
  const vec<T, 3> dtbeta;

  // Spatial derivatives of ADM variables
  const smat<vec<T, 3>, 3> dgamma;
  const smat<smat<T, 3>, 3> ddgamma;
  const vec<T, 3> dalpha;

  // Second time derivatives of ADM variables
  const smat<T, 3> dtk;
  const T dt2alpha;
  const vec<T, 3> dt2beta;

  // Space-time derivatives of ADM variables
  const smat<vec<T, 3>, 3> dk;
  const vec<T, 3> ddtalpha;
  const vec<vec<T, 3>, 3> ddtbeta;

  // Second spatial derivatives of ADM variables
  const smat<T, 3> ddalpha;
  const vec<vec<T, 3>, 3> dbeta;
  const vec<smat<T, 3>, 3> ddbeta;

  // Intermediate quantities
  const mat<T, 3> dtgamma;
  const vec<T, 3> dtbetal;
  const vec<vec<T, 3>, 3> dbetal;
  const T dtgtt;
  const vec<T, 3> dgtt;

  const smat<T, 3> dt2gamma;
  const smat<vec<T, 3>, 3> ddtgamma;
  const vec<T, 3> dt2betal;
  const vec<vec<T, 3>, 3> ddtbetal;
  const vec<smat<T, 3>, 3> ddbetal;
  const T dt2gtt;
  const vec<T, 3> ddtgtt;
  const smat<T, 3> ddgtt;

  // Derivatives of 4-metric
  const smat<vec<T, 4>, 4> dg;
  const smat<smat<T, 4>, 4> ddg;
  const smat<vec<T, 4>, 4> dgu;

  // Christoffel symbol
  const vec<smat<T, 4>, 4> Gammal;
  const vec<smat<T, 4>, 4> Gamma;

  const vec<smat<vec<T, 4>, 4>, 4> dGammal;
  const vec<smat<vec<T, 4>, 4>, 4> dGamma;

  // Riemann, Ricci, Weyl
  const rten<T, 4> Rm;
  const smat<T, 4> R;
  const T Rsc;
  // TODO: C is trace free
  const rten<T, 4> C;

  // Ricci and Weyl scalars
  const T Lambda;
  const T Phi00, Phi11, Phi22;
  const cplx<T> Phi10, Phi20, Phi21;
  const cplx<T> Psi0, Psi1, Psi2, Psi3, Psi4;

  // Gradient of tetrad
  const vec<vec<T, 4>, 4> det, dephi, detheta, der;
  const vec<vec<T, 4>, 4> dl, dn;
  const vec<vec<cplx<T>, 4>, 4> dm;

  // Newman-Penrose spin coefficients
  const cplx<T> npkappa, npsigma, nprho, nptau, npepsilon, npbeta, npalpha,
      npgamma, nppi, npmu, nplambda, npnu;

  inline ARITH_INLINE ARITH_DEVICE ARITH_HOST
  weyl_vars(const T time, const vec<T, 3> &coord3,
            //
            const smat<T, 3> &gamma, const T &alpha, const vec<T, 3> &beta,
            //
            const mat<T, 3> &k, const T &dtalpha, const vec<T, 3> &dtbeta,
            //
            const mat<vec<T, 3>, 3> &dgamma, const vec<T, 3> &dalpha,
            const vec<vec<T, 3>, 3> &dbeta,
            //
            const mat<T, 3> &dtk, const T &dt2alpha, const vec<T, 3> &dt2beta,
            //
            const mat<vec<T, 3>, 3> &dk, const vec<T, 3> &ddtalpha,
            const vec<vec<T, 3>, 3> &ddtbeta,
            //
            const mat<mat<T, 3>, 3> &ddgamma, const mat<T, 3> &ddalpha,
            const vec<mat<T, 3>, 3> &ddbeta)
      : weyl_vars_noderivs<T>(time, coord3, gamma, alpha, beta),
        // Time derivatives of ADM variables
        k(k), dtalpha(dtalpha), dtbeta(dtbeta),
        // Spatial derivatives of ADM variables
        dgamma(dgamma), ddgamma(ddgamma), dalpha(dalpha),
        // Second time derivatives of ADM variables
        dtk(dtk), dt2alpha(dt2alpha), dt2beta(dt2beta),
        // Space-time derivatives of ADM variables
        dk(dk), ddtalpha(ddtalpha), ddtbeta(ddtbeta),
        // Second spatial derivatives of ADM variables
        ddalpha(ddalpha), dbeta(dbeta), ddbeta(ddbeta),
        //
        // dt gamma_ij = -2 alpha K_ij
        //               + gamma_kj beta^k,i + gamma_ik beta^k,j
        //               + beta^k gamma_ij,k
        dtgamma([&](int a, int b) ARITH_INLINE {
          return -2 * alpha * k(a, b) //
                 + sum<3>([&](int x) ARITH_INLINE {
                     return gamma(x, b) * dbeta(x)(a);
                   }) //
                 + sum<3>([&](int x) ARITH_INLINE {
                     return gamma(a, x) * dbeta(x)(b);
                   }) //
                 + sum<3>([&](int x) ARITH_INLINE {
                     return beta(x) * dgamma(a, b)(x);
                   });
        }),
        dtbetal([&](int a) ARITH_INLINE {
          return sum<3>([&](int x) ARITH_INLINE {
            return dtgamma(a, x) * beta(x) + gamma(a, x) * dtbeta(x);
          });
        }),
        dbetal([&](int a) ARITH_INLINE {
          return vec<T, 3>([&](int b) ARITH_INLINE {
            return sum<3>([&](int x) ARITH_INLINE {
              return dgamma(a, x)(b) * beta(x) + gamma(a, x) * dbeta(x)(b);
            });
          });
        }),
        dtgtt(
            -2 * alpha * dtalpha //
            +
            sum<3>([&](int x) ARITH_INLINE { return dtbetal(x) * beta(x); }) //
            + sum<3>([&](int x) ARITH_INLINE { return betal(x) * dtbeta(x); })),
        dgtt([&](int a) {
          return -2 * alpha * dalpha(a) //
                 + sum<3>([&](int x) ARITH_INLINE {
                     return dbetal(x)(a) * beta(x);
                   }) //
                 + sum<3>([&](int x)
                              ARITH_INLINE { return betal(x) * dbeta(x)(a); });
        }),
        //
        dt2gamma([&](int a, int b) ARITH_INLINE {
          return -2 * dtalpha * k(a, b)  //
                 - 2 * alpha * dtk(a, b) //
                 + sum<3>([&](int x) ARITH_INLINE {
                     return dtgamma(x, b) * dbeta(x)(a) +
                            gamma(x, b) * ddtbeta(x)(a);
                   }) //
                 + sum<3>([&](int x) ARITH_INLINE {
                     return dtgamma(a, x) * dbeta(x)(b) +
                            gamma(a, x) * ddtbeta(x)(b);
                   }) //
                 + sum<3>([&](int x) ARITH_INLINE {
                     return dtbeta(x) * dgamma(a, b)(x) +
                            beta(x) * ddtgamma(a, b)(x);
                   });
        }),
        ddtgamma([&](int a, int b) ARITH_INLINE {
          return vec<T, 3>([&](int c) ARITH_INLINE {
            return -2 * dalpha(c) * k(a, b)  //
                   - 2 * alpha * dk(a, b)(c) //
                   + sum<3>([&](int x) ARITH_INLINE {
                       return dgamma(x, b)(c) * dbeta(x)(a) +
                              gamma(x, b) * ddbeta(x)(a, c);
                     }) //
                   + sum<3>([&](int x) ARITH_INLINE {
                       return dgamma(a, x)(c) * dbeta(x)(b) +
                              gamma(a, x) * ddbeta(x)(b, c);
                     }) //
                   + sum<3>([&](int x) ARITH_INLINE {
                       return dbeta(x)(c) * dgamma(a, b)(x) +
                              beta(x) * ddgamma(a, b)(x, c);
                     });
          });
        }),
        dt2betal([&](int a) ARITH_INLINE {
          return sum<3>([&](int x) ARITH_INLINE {
            return dt2gamma(a, x) * beta(x)        //
                   + 2 * dtgamma(a, x) * dtbeta(x) //
                   + dtgamma(a, x) * dt2beta(x);
          });
        }),
        ddtbetal([&](int a) ARITH_INLINE {
          return vec<T, 3>([&](int b) ARITH_INLINE {
            return sum<3>([&](int x) ARITH_INLINE {
              return ddtgamma(a, x)(b) * beta(x)   //
                     + dgamma(a, x)(b) * dtbeta(x) //
                     + dtgamma(a, x) * dbeta(x)(b) //
                     + gamma(a, x) * ddtbeta(x)(b);
            });
          });
        }),
        ddbetal([&](int a) ARITH_INLINE {
          return mat<T, 3>([&](int b, int c) ARITH_INLINE {
            return sum<3>([&](int x) ARITH_INLINE {
              return ddgamma(a, x)(b, c) * beta(x)   //
                     + dgamma(a, x)(b) * dbeta(x)(c) //
                     + dgamma(a, x)(c) * dbeta(x)(b) //
                     + gamma(a, x) * ddbeta(x)(b, c);
            });
          });
        }),
        dt2gtt(-2 * pow2(dtalpha)     //
               - 2 * alpha * dt2alpha //
               + sum<3>([&](int x)
                            ARITH_INLINE { return dt2betal(x) * beta(x); }) //
               + 2 * sum<3>([&](int x) ARITH_INLINE {
                   return dtbetal(x) * dtbeta(x);
                 }) //
               + sum<3>([&](int x)
                            ARITH_INLINE { return betal(x) * dt2beta(x); })),
        ddtgtt([&](int a) {
          return -2 * dalpha(a) * dtalpha  //
                 - 2 * alpha * ddtalpha(a) //
                 + sum<3>([&](int x) ARITH_INLINE {
                     return ddtbetal(x)(a) * beta(x);
                   }) //
                 + sum<3>([&](int x) ARITH_INLINE {
                     return dtbetal(x) * dbeta(x)(a);
                   }) //
                 + sum<3>([&](int x) ARITH_INLINE {
                     return dbetal(x)(a) * dtbeta(x);
                   }) //
                 + sum<3>([&](int x) ARITH_INLINE {
                     return betal(x) * ddtbeta(x)(a);
                   });
        }),
        ddgtt([&](int a, int b) {
          return -2 * dalpha(a) * dalpha(b)  //
                 - 2 * alpha * ddalpha(a, b) //
                 + sum<3>([&](int x) ARITH_INLINE {
                     return ddbetal(x)(a, b) * beta(x);
                   }) //
                 + sum<3>([&](int x) ARITH_INLINE {
                     return dbetal(x)(a) * dbeta(x)(b);
                   }) //
                 + sum<3>([&](int x) ARITH_INLINE {
                     return dbetal(x)(b) * dbeta(x)(a);
                   }) //
                 + sum<3>([&](int x) ARITH_INLINE {
                     return betal(x) * ddbeta(x)(a, b);
                   });
        }),
        //
        dg([&](int a, int b) ARITH_INLINE {
          return vec<T, 4>([&](int c) ARITH_INLINE {
            if (c == 0) {
              if (a == 0 && b == 0)
                return dtgtt;
              if (a == 0)
                return dtbetal(b - 1);
              if (b == 0)
                return dtbetal(a - 1);
              return dtgamma(a - 1, b - 1);
            }
            if (a == 0 && b == 0)
              return dgtt(c - 1);
            if (a == 0)
              return dbetal(b - 1)(c - 1);
            if (b == 0)
              return dbetal(a - 1)(c - 1);
            return dgamma(a - 1, b - 1)(c - 1);
          });
        }),
        //
        ddg([&](int a, int b) ARITH_INLINE {
          return mat<T, 4>([&](int c, int d) ARITH_INLINE {
            if (c == 0 && d == 0) {
              if (a == 0 && b == 0)
                return dt2gtt;
              if (a == 0)
                return dt2betal(b - 1);
              if (b == 0)
                return dt2betal(a - 1);
              return dt2gamma(a - 1, b - 1);
            }
            if (c == 0) {
              if (a == 0 && b == 0)
                return ddtgtt(d - 1);
              if (a == 0)
                return ddtbetal(b - 1)(d - 1);
              if (b == 0)
                return ddtbetal(a - 1)(d - 1);
              return ddtgamma(a - 1, b - 1)(d - 1);
            }
            if (d == 0) {
              if (a == 0 && b == 0)
                return ddtgtt(c - 1);
              if (a == 0)
                return ddtbetal(b - 1)(c - 1);
              if (b == 0)
                return ddtbetal(a - 1)(c - 1);
              return ddtgamma(a - 1, b - 1)(c - 1);
            }
            if (a == 0 && b == 0)
              return ddgtt(c - 1, d - 1);
            if (a == 0)
              return ddbetal(b - 1)(c - 1, d - 1);
            if (b == 0)
              return ddbetal(a - 1)(c - 1, d - 1);
            return ddgamma(a - 1, b - 1)(c - 1, d - 1);
          });
        }),
        //
        dgu(calc_dgu(gu, dg)),
        //
        Gammal(calc_gammal(dg)),       //
        Gamma(calc_gamma(gu, Gammal)), //
        //
        dGammal(calc_dgammal(ddg)),                    //
        dGamma(calc_dgamma(gu, dgu, Gammal, dGammal)), //
        //
        Rm(calc_riemann(g, Gamma, dGamma)), //
        R(calc_ricci(gu, Rm)),              //
        Rsc(calc_trace(R, gu)),             //
        C(calc_weyl(g, Rm, R, Rsc)),
        //
        // Badri Krishnan's PhD thesis, appendix A
        Lambda(Rsc / 24), //
        Phi00(sum_symm<4>([&](int a, int b) ARITH_INLINE {
          return R(a, b) * l(a) * l(b) / 2;
        })),
        Phi11(sum<4>([&](int a, int b) ARITH_INLINE {
          return R(a, b) * (l(a) * n(b) + real(m(a) * conj(m(b)))) / 4;
        })),
        Phi22(sum_symm<4>([&](int a, int b) ARITH_INLINE {
          return R(a, b) * n(a) * n(b) / 2;
        })),
        Phi10(sum<4>([&](int a, int b) ARITH_INLINE {
          return R(a, b) * l(a) * conj(m(b)) / T(2);
        })),
        Phi20(sum_symm<4>([&](int a, int b) ARITH_INLINE {
          return R(a, b) * conj(m(a)) * conj(m(b)) / T(2);
        })),
        Phi21(sum<4>([&](int a, int b) ARITH_INLINE {
          return R(a, b) * conj(m(a)) * n(b) / T(2);
        })),
        Psi0(sum<4>([&](int b) ARITH_INLINE {
          // sum<4>([&](int a, int b, int c, int
          // d) ARITH_INLINE {
          //   return C(a, b, c, d) * l(a) * m(b) * l(c) * m(d);
          // })
          return m(b) * sum<4>([&](int d) ARITH_INLINE {
                   return m(d) * sum<4>([&](int a) ARITH_INLINE {
                            return l(a) * sum<4>([&](int c) ARITH_INLINE {
                                     return C(a, b, c, d) * l(c);
                                   });
                          });
                 });
        })),
        Psi1(
            // sum<4>([&](int a, int b, int c, int
            // d) ARITH_INLINE {
            //   return C(a, b, c, d) * l(a) * m(b) * l(c) * n(d);
            // })
            sum<4>([&](int b) ARITH_INLINE {
              return m(b) * sum<4>([&](int a) ARITH_INLINE {
                       return l(a) * sum<4>([&](int c) ARITH_INLINE {
                                return l(c) * sum<4>([&](int d) ARITH_INLINE {
                                         return C(a, b, c, d) * n(d);
                                       });
                              });
                     });
            })),
        Psi2(
            // sum<4>([&](int a, int b, int c, int
            // d) ARITH_INLINE {
            //   return C(a, b, c, d) * l(a) * m(b) * conj(m(c)) * n(d);
            // })
            sum<4>([&](int b) ARITH_INLINE {
              return m(b) * sum<4>([&](int c) ARITH_INLINE {
                       return conj(m(c)) * sum<4>([&](int a) ARITH_INLINE {
                                return l(a) * sum<4>([&](int d) ARITH_INLINE {
                                         return C(a, b, c, d) * n(d);
                                       });
                              });
                     });
            })),
        Psi3(
            // sum<4>([&](int a, int b, int c, int
            // d) ARITH_INLINE {
            //   return C(a, b, c, d) * l(a) * n(b) * conj(m(c)) * n(d);
            // })
            sum<4>([&](int c) ARITH_INLINE {
              return conj(m(c)) * sum<4>([&](int a) ARITH_INLINE {
                       return l(a) * sum<4>([&](int b) ARITH_INLINE {
                                return n(b) * sum<4>([&](int d) ARITH_INLINE {
                                         return C(a, b, c, d) * n(d);
                                       });
                              });
                     });
            })),
        Psi4(
            // sum<4>([&](int a, int b, int c, int
            // d) ARITH_INLINE {
            //   return C(a, b, c, d) * conj(m(a)) * n(b) * conj(m(c)) * n(d);
            // })
            sum<4>([&](int a) ARITH_INLINE {
              return conj(m(a)) * sum<4>([&](int c) ARITH_INLINE {
                       return conj(m(c)) * sum<4>([&](int b) ARITH_INLINE {
                                return n(b) * sum<4>([&](int d) ARITH_INLINE {
                                         return C(a, b, c, d) * n(d);
                                       });
                              });
                     });
            })),
        det(calc_det(gu, dgu, et, Gamma)),                                    //
        dephi(calc_dephi(coord, g, dg, ephi, Gamma)),                         //
        detheta(calc_detheta(coord, g, dg, ephi, dephi, etheta, Gamma)),      //
        der(calc_der(coord, g, dg, ephi, dephi, etheta, detheta, er, Gamma)), //
        dl([&](int a) ARITH_INLINE { return (det(a) + der(a)) / sqrt(T(2)); }),
        dn([&](int a) ARITH_INLINE { return (det(a) - der(a)) / sqrt(T(2)); }),
        dm([&](int a) ARITH_INLINE {
          return vec<cplx<T>, 4>([&](int b) ARITH_INLINE {
            return cplx<T>(detheta(a)(b), dephi(a)(b)) / sqrt(T(2));
          });
        }),
        npkappa(sum<4>([&](int a, int b)
                           ARITH_INLINE { return -m(a) * l(b) * dl(a)(b); })),
        npsigma(sum<4>([&](int a, int b)
                           ARITH_INLINE { return -m(a) * m(b) * dl(a)(b); })),
        nprho(sum<4>([&](int a, int b) ARITH_INLINE {
          return -m(a) * conj(m(b)) * dl(a)(b);
        })),
        nptau(sum<4>([&](int a, int b)
                         ARITH_INLINE { return -m(a) * n(b) * dl(a)(b); })),
        npepsilon(sum<4>([&](int a, int b) ARITH_INLINE {
          return (conj(m(a)) * l(b) * dm(a)(b) - n(a) * l(b) * dl(a)(b)) / T(2);
        })),
        npbeta(sum<4>([&](int a, int b) ARITH_INLINE {
          return (conj(m(a)) * m(b) * dm(a)(b) - n(a) * m(b) * dl(a)(b)) / T(2);
        })),
        npalpha(sum<4>([&](int a, int b) ARITH_INLINE {
          return (conj(m(a)) * conj(m(b)) * dm(a)(b) -
                  n(a) * conj(m(b)) * dl(a)(b)) /
                 T(2);
        })),
        npgamma(sum<4>([&](int a, int b) ARITH_INLINE {
          return (conj(m(a)) * n(b) * dm(a)(b) - n(a) * n(b) * dl(a)(b)) / T(2);
        })),
        nppi(sum<4>([&](int a, int b)
                        ARITH_INLINE { return conj(m(a)) * l(b) * dn(a)(b); })),
        npmu(sum<4>([&](int a, int b)
                        ARITH_INLINE { return conj(m(a)) * m(b) * dn(a)(b); })),
        nplambda(sum<4>([&](int a, int b) ARITH_INLINE {
          return conj(m(a)) * conj(m(b)) * dn(a)(b);
        })),
        npnu(sum<4>([&](int a, int b)
                        ARITH_INLINE { return conj(m(a)) * n(b) * dn(a)(b); }))
  //
  {}

  inline ARITH_INLINE ARITH_DEVICE ARITH_HOST
  weyl_vars(const vec<T, 4> &coord, //
            const mat<T, 4> &g, const mat<vec<T, 4>, 4> &dg,
            const mat<mat<T, 4>, 4> &ddg)
      : weyl_vars_noderivs<T>(coord, g),
        //
        k(), dtalpha(), dtbeta(), dgamma(), ddgamma(), dalpha(), dtk(),
        dt2alpha(), dt2beta(), dk(), ddtalpha(), ddtbeta(), ddalpha(), dbeta(),
        ddbeta(), dtgamma(), dtbetal(), dbetal(), dt2gamma(), ddtgamma(),
        dt2betal(), ddtbetal(), ddbetal(),
        //
        dg(dg), ddg(ddg),
        //
        dgu(calc_dgu(gu, dg)),
        //
        Gammal(calc_gammal(dg)),       //
        Gamma(calc_gamma(gu, Gammal)), //
        //
        dGammal(calc_dgammal(ddg)),                    //
        dGamma(calc_dgamma(gu, dgu, Gammal, dGammal)), //
        //
        Rm(calc_riemann(g, Gamma, dGamma)), //
        R(calc_ricci(gu, Rm)),              //
        Rsc(R.trace(gu)),                   //
        C(calc_weyl(g, Rm, R, Rsc)),
        //
        // Badri Krishnan's PhD thesis, appendix A
        Lambda(Rsc / 24), //
        Phi00(sum_symm<4>([&](int a, int b) ARITH_INLINE {
          return R(a, b) * l(a) * l(b) / 2;
        })),
        Phi11(sum<4>([&](int a, int b) ARITH_INLINE {
          return R(a, b) * (l(a) * n(b) + real(m(a) * conj(m(b)))) / 4;
        })),
        Phi22(sum_symm<4>([&](int a, int b) ARITH_INLINE {
          return R(a, b) * n(a) * n(b) / 2;
        })),
        Phi10(sum<4>([&](int a, int b) ARITH_INLINE {
          return R(a, b) * l(a) * conj(m(b)) / T(2);
        })),
        Phi20(sum_symm<4>([&](int a, int b) ARITH_INLINE {
          return R(a, b) * conj(m(a)) * conj(m(b)) / T(2);
        })),
        Phi21(sum<4>([&](int a, int b) ARITH_INLINE {
          return R(a, b) * conj(m(a)) * n(b) / T(2);
        })),
        Psi0(sum<4>([&](int b) ARITH_INLINE {
          // sum<4>([&](int a, int b, int c, int
          // d) ARITH_INLINE {
          //   return C(a, b, c, d) * l(a) * m(b) * l(c) * m(d);
          // })
          return m(b) * sum<4>([&](int d) ARITH_INLINE {
                   return m(d) * sum<4>([&](int a) ARITH_INLINE {
                            return l(a) * sum<4>([&](int c) ARITH_INLINE {
                                     return C(a, b, c, d) * l(c);
                                   });
                          });
                 });
        })),
        Psi1(
            // sum<4>([&](int a, int b, int c, int
            // d) ARITH_INLINE {
            //   return C(a, b, c, d) * l(a) * m(b) * l(c) * n(d);
            // })
            sum<4>([&](int b) ARITH_INLINE {
              return m(b) * sum<4>([&](int a) ARITH_INLINE {
                       return l(a) * sum<4>([&](int c) ARITH_INLINE {
                                return l(c) * sum<4>([&](int d) ARITH_INLINE {
                                         return C(a, b, c, d) * n(d);
                                       });
                              });
                     });
            })),
        Psi2(
            // sum<4>([&](int a, int b, int c, int
            // d) ARITH_INLINE {
            //   return C(a, b, c, d) * l(a) * m(b) * conj(m(c)) * n(d);
            // })
            sum<4>([&](int b) ARITH_INLINE {
              return m(b) * sum<4>([&](int c) ARITH_INLINE {
                       return conj(m(c)) * sum<4>([&](int a) ARITH_INLINE {
                                return l(a) * sum<4>([&](int d) ARITH_INLINE {
                                         return C(a, b, c, d) * n(d);
                                       });
                              });
                     });
            })),
        Psi3(
            // sum<4>([&](int a, int b, int c, int
            // d) ARITH_INLINE {
            //   return C(a, b, c, d) * l(a) * n(b) * conj(m(c)) * n(d);
            // })
            sum<4>([&](int c) ARITH_INLINE {
              return conj(m(c)) * sum<4>([&](int a) ARITH_INLINE {
                       return l(a) * sum<4>([&](int b) ARITH_INLINE {
                                return n(b) * sum<4>([&](int d) ARITH_INLINE {
                                         return C(a, b, c, d) * n(d);
                                       });
                              });
                     });
            })),
        Psi4(
            // sum<4>([&](int a, int b, int c, int
            // d) ARITH_INLINE {
            //   return C(a, b, c, d) * conj(m(a)) * n(b) * conj(m(c)) * n(d);
            // })
            sum<4>([&](int a) ARITH_INLINE {
              return conj(m(a)) * sum<4>([&](int c) ARITH_INLINE {
                       return conj(m(c)) * sum<4>([&](int b) ARITH_INLINE {
                                return n(b) * sum<4>([&](int d) ARITH_INLINE {
                                         return C(a, b, c, d) * n(d);
                                       });
                              });
                     });
            })),
        det(calc_det(gu, dgu, et, Gamma)),                                    //
        dephi(calc_dephi(coord, g, dg, ephi, Gamma)),                         //
        detheta(calc_detheta(coord, g, dg, ephi, dephi, etheta, Gamma)),      //
        der(calc_der(coord, g, dg, ephi, dephi, etheta, detheta, er, Gamma)), //
        dl([&](int a) ARITH_INLINE { return (det(a) + der(a)) / sqrt(T(2)); }),
        dn([&](int a) ARITH_INLINE { return (det(a) - der(a)) / sqrt(T(2)); }),
        dm([&](int a) ARITH_INLINE {
          return vec<cplx<T>, 4>([&](int b) ARITH_INLINE {
            return cplx<T>(detheta(a)(b), dephi(a)(b)) / sqrt(T(2));
          });
        }),
        npkappa(sum<4>([&](int a, int b)
                           ARITH_INLINE { return -m(a) * l(b) * dl(a)(b); })),
        npsigma(sum<4>([&](int a, int b)
                           ARITH_INLINE { return -m(a) * m(b) * dl(a)(b); })),
        nprho(sum<4>([&](int a, int b) ARITH_INLINE {
          return -m(a) * conj(m(b)) * dl(a)(b);
        })),
        nptau(sum<4>([&](int a, int b)
                         ARITH_INLINE { return -m(a) * n(b) * dl(a)(b); })),
        npepsilon(sum<4>([&](int a, int b) ARITH_INLINE {
          return (conj(m(a)) * l(b) * dm(a)(b) - n(a) * l(b) * dl(a)(b)) / T(2);
        })),
        npbeta(sum<4>([&](int a, int b) ARITH_INLINE {
          return (conj(m(a)) * m(b) * dm(a)(b) - n(a) * m(b) * dl(a)(b)) / T(2);
        })),
        npalpha(sum<4>([&](int a, int b) ARITH_INLINE {
          return (conj(m(a)) * conj(m(b)) * dm(a)(b) -
                  n(a) * conj(m(b)) * dl(a)(b)) /
                 T(2);
        })),
        npgamma(sum<4>([&](int a, int b) ARITH_INLINE {
          return (conj(m(a)) * n(b) * dm(a)(b) - n(a) * n(b) * dl(a)(b)) / T(2);
        })),
        nppi(sum<4>([&](int a, int b)
                        ARITH_INLINE { return conj(m(a)) * l(b) * dn(a)(b); })),
        npmu(sum<4>([&](int a, int b)
                        ARITH_INLINE { return conj(m(a)) * m(b) * dn(a)(b); })),
        nplambda(sum<4>([&](int a, int b) ARITH_INLINE {
          return conj(m(a)) * conj(m(b)) * dn(a)(b);
        })),
        npnu(sum<4>([&](int a, int b)
                        ARITH_INLINE { return conj(m(a)) * n(b) * dn(a)(b); }))
  //
  {}
};

} // namespace Weyl

#endif // #ifndef WEYL_HXX
