#ifndef CARPETX_CARPETX_TYPEPROPS_HH
#define CARPETX_CARPETX_TYPEPROPS_HH

#include <cctk.h>

namespace CarpetLib {

template <typename T> struct typeprops {
  typedef T complex;
  typedef T real;
  // Create a complex number from a real number
  static inline complex const &fromreal(real const &x) { return x; }
};

#ifdef HAVE_CCTK_COMPLEX8
template <> struct typeprops<CCTK_COMPLEX8> {
  typedef CCTK_COMPLEX8 complex;
  typedef CCTK_REAL4 real;
  static inline complex fromreal(real const x) { return CCTK_Cmplx8(x, 0); }
};
#endif

#ifdef HAVE_CCTK_COMPLEX16
template <> struct typeprops<CCTK_COMPLEX16> {
  typedef CCTK_COMPLEX16 complex;
  typedef CCTK_REAL8 real;
  static inline complex fromreal(real const x) { return CCTK_Cmplx16(x, 0); }
};
#endif

#ifdef HAVE_CCTK_COMPLEX32
template <> struct typeprops<CCTK_COMPLEX32> {
  typedef CCTK_COMPLEX32 complex;
  typedef CCTK_REAL16 real;
  static inline complex fromreal(real const x) { return CCTK_Cmplx32(x, 0); }
};
#endif

// Return the specific Cactus variable type for a Cactus variable type
static inline int specific_cactus_type(int const vartype) {
  switch (vartype) {
  case CCTK_VARIABLE_INT:
#ifdef CCTK_INTEGER_PRECISION_1
    return CCTK_VARIABLE_INT1;
#endif
#ifdef CCTK_INTEGER_PRECISION_2
    return CCTK_VARIABLE_INT2;
#endif
#ifdef CCTK_INTEGER_PRECISION_4
    return CCTK_VARIABLE_INT4;
#endif
#ifdef CCTK_INTEGER_PRECISION_8
    return CCTK_VARIABLE_INT8;
#endif
#ifdef CCTK_INTEGER_PRECISION_16
    return CCTK_VARIABLE_INT16;
#endif
    return -1;
  case CCTK_VARIABLE_REAL:
#ifdef CCTK_REAL_PRECISION_4
    return CCTK_VARIABLE_REAL4;
#endif
#ifdef CCTK_REAL_PRECISION_8
    return CCTK_VARIABLE_REAL8;
#endif
#ifdef CCTK_REAL_PRECISION_16
    return CCTK_VARIABLE_REAL16;
#endif
    return -1;
  case CCTK_VARIABLE_COMPLEX:
#ifdef CCTK_REAL_PRECISION_4
    return CCTK_VARIABLE_COMPLEX8;
#endif
#ifdef CCTK_REAL_PRECISION_8
    return CCTK_VARIABLE_COMPLEX16;
#endif
#ifdef CCTK_REAL_PRECISION_16
    return CCTK_VARIABLE_COMPLEX32;
#endif
    return -1;
  }
  return vartype;
}
} // namespace CarpetLib

#endif // #ifndef CARPETX_CARPETX_TYPEPROPS_HH
