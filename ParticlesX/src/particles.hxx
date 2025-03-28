#ifndef PARTICLESX_HXX
#define PARTICLESX_HXX

#include <vect.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>

#include <vector>

namespace ParticlesX {

constexpr int dim = 3;

struct Particle {
  CCTK_REAL id;                    // id (integer)
  CCTK_REAL tau;                   // eigentime
  Arith::vect<CCTK_REAL, dim> pos; // position
  Arith::vect<CCTK_REAL, dim> vel; // velocity
  Arith::vect<CCTK_REAL, dim> acc; // acceleration
};

extern "C" void ParticlesX_output(CCTK_ARGUMENTS);

class Particles {
  std::vector<CCTK_REAL> id;
  std::vector<CCTK_REAL> tau;
  Arith::vect<std::vector<CCTK_REAL>, dim> pos;
  Arith::vect<std::vector<CCTK_REAL>, dim> vel;
  Arith::vect<std::vector<CCTK_REAL>, dim> acc;

  friend void ParticlesX_output(CCTK_ARGUMENTS);

public:
  Particles() = delete;

  Particles(const Particles &) = delete;
  Particles(Particles &&) = delete;
  Particles &operator=(const Particles &) = delete;
  Particles &operator=(Particles &&) = delete;

  Particles(int n);
  int size() const;

  Particle get(int n) const;
  void set(const Particle &p, int n);
};

Particles &get_particles();

} // namespace ParticlesX

#endif // #ifndef PARTICLESX_HXX
