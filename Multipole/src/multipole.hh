#include "cctk.h"
#include "cctk_Arguments.h"

#include <string>
#include <vector>
#include <cassert>

// Multipole_Calc
//      This is the main scheduling file.  Because we are completely local here
//      and do not use cactus arrays etc, we schedule only one function and then
//      like program like one would in C, C++ with this function taking the 
//      place of int main(void).
//
//      This function calls functions to accomplish 3 things:
//        1) Interpolate psi4 onto a sphere
//        2) Integrate psi4 with the ylm's over that sphere
//        2) Output the mode decomposed psi4
extern "C" void Multipole_Calc(CCTK_ARGUMENTS);

namespace Multipole {
// information about variable which we decompose
struct variable_desc
{
  int index;
  int imag_index;
  CCTK_INT spin_weight;
  std::string name;
};

// a simple array class to hold complex modes abs(m) <= l, l <= lmax, for
// nradii radii for nvars variables
class mode_array {
  public:
    mode_array(int nvars_, int nradii_, int lmax_) : nvars(nvars_),
      nradii(nradii_), lmax(lmax_),
      modes(size_t(nvars * nradii * (lmax+1)*(lmax+1) * 2)) {}
    virtual ~mode_array() {}
    // default copy and assignment is ok

    CCTK_REAL& operator()(int v, int ri, int l, int m, bool is_im) {
      return modes.at(mode_idx(v, ri, l, m, is_im));
    }

    const CCTK_REAL& operator()(int v, int ri, int l, int m, bool is_im) const {
      return modes.at(mode_idx(v, ri, l, m, is_im));
    }

    int get_nvars() const { return nvars; }
    int get_nradii() const { return nradii; }
    int get_lmax() const { return lmax; }
  private:
    size_t mode_idx(int v, int ri, int l, int m, int is_im) const {
      assert(v >= 0 && v < nvars);
      assert(ri >= 0 && ri < nradii);
      assert(l >= 0 && l <= lmax);
      assert(m <= l && -m <= l);
      return size_t(v * nradii * (lmax+1)*(lmax+1) * 2 +
                    ri * (lmax+1)*(lmax+1) * 2 +
                    (l*l + (m+l)) * 2 + is_im);
    }

    const int nvars;
    const int nradii;
    const int lmax;
    std::vector<CCTK_REAL> modes;
};
}
