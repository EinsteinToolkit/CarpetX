#ifndef FIXMATH_HXX
#define FIXMATH_HXX

#if 0

// Include this file BEFORE including <cctk.h>, best before including
// any other include files

#ifdef _CCTK_MATH_H_
#error                                                                         \
    "The Cactus include file <cctk_Math.h> has already been included. However, this file <fixmath.hxx> needs to be included first. It is best to include <fixmath.hxx> before any Cactus include files."
#endif

// This provides broken `#defines` for `isnan` etc.
#include <cctk.h>

// Disable the broken definitions
#undef copysign
#undef fpclassify
#undef isfinite
#undef isinf
#undef isnan
#undef isnormal
#undef signbit

#endif

#endif // #ifndef FIXMATH_HXX
