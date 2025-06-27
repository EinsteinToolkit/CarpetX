#include "particles.hxx"

#include <vect.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#ifndef HAVE_CAPABILITY_ADIOS2
#error "Particles require ADIOS2 support for I/O"
#endif
#include <adios2.h>
#if !defined ADIOS2_USE_MPI || !ADIOS2_USE_MPI
#error                                                                         \
    "ParticlesX requires an MPI-enabled ADIOS2 library. Please compile ADIOS2 with MPI support and enable its use by defining ADIOS2_USE_MPI."
#endif

#ifndef HAVE_CAPABILITY_MPI
#error "Particles requires MPI support for I/O"
#endif
#include <mpi.h>

#include <iomanip>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace ParticlesX {

Particles::Particles(int n)
    : id(n), tau(n), pos{std::vector<CCTK_REAL>(n), std::vector<CCTK_REAL>(n),
                         std::vector<CCTK_REAL>(n)},
      vel{std::vector<CCTK_REAL>(n), std::vector<CCTK_REAL>(n),
          std::vector<CCTK_REAL>(n)},
      acc{std::vector<CCTK_REAL>(n), std::vector<CCTK_REAL>(n),
          std::vector<CCTK_REAL>(n)} {}

int Particles::size() const { return id.size(); }

Particle Particles::get(int n) const {
  return Particle{
      id.at(n),
      tau.at(n),
      {pos[0].at(n), pos[1].at(n), pos[2].at(n)},
      {vel[0].at(n), vel[1].at(n), vel[2].at(n)},
      {acc[0].at(n), acc[1].at(n), acc[2].at(n)},
  };
}

void Particles::set(const Particle &p, int n) {
  id.at(n) = p.id;
  tau.at(n) = p.tau;
  for (int d = 0; d < dim; ++d)
    pos[d].at(n) = p.pos[d];
  for (int d = 0; d < dim; ++d)
    vel[d].at(n) = p.vel[d];
  for (int d = 0; d < dim; ++d)
    acc[d].at(n) = p.acc[d];
}

std::unique_ptr<Particles> particles_ptr;

Particles &get_particles() { return *particles_ptr; }

////////////////////////////////////////////////////////////////////////////////

extern "C" void ParticlesX_setup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ParticlesX_setup;
  DECLARE_CCTK_PARAMETERS;

  const int myproc = CCTK_MyProc(cctkGH);
  // const int nprocs = CCTK_nProcs(cctkGH);
  const int nparticles = myproc == 0 ? num_particles : 0;

  particles_ptr = std::make_unique<Particles>(nparticles);
}

extern "C" void ParticlesX_output(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ParticlesX_output;
  DECLARE_CCTK_PARAMETERS;

  if (out_every <= 0 || cctk_iteration % out_every != 0)
    return;

  Particles &particles = get_particles();
  const int nparticles = particles.size();

  std::ostringstream buf;
  buf << out_dir << "/"                                               //
      << "particles"                                                  //
      << ".it" << std::setw(8) << std::setfill('0') << cctk_iteration //
      << ".bp5";
  const std::string filename = buf.str();

  adios2::ADIOS adios = adios2::ADIOS("", MPI_COMM_WORLD, "Fortran");
  adios2::IO io = adios.DeclareIO("IO");
  adios2::Engine engine = io.Open(filename, adios2::Mode::Write);

#ifdef ADIOS2_HAVE_BLOSC2
  const adios2::Operator compressor =
      adios.DefineOperator("Blosc2Compressor", adios2::ops::LosslessBlosc);
  // Use a high compression rate and a byteshuffle filter
  const adios2::Params compressor_options{
      {adios2::ops::blosc::key::clevel, adios2::ops::blosc::value::clevel_9},
      {adios2::ops::blosc::key::doshuffle,
       adios2::ops::blosc::value::doshuffle_shuffle},
  };
#endif

  auto output_variable = [&](const std::string &name,
                             const std::vector<CCTK_REAL> &data) {
    adios2::Variable<CCTK_REAL> var =
        io.DefineVariable<CCTK_REAL>(name, {}, {}, {std::size_t(nparticles)});
#ifdef ADIOS2_HAVE_BLOSC2
    var.AddOperation(compressor, compressor_options);
#endif
    var.SetSelection({{}, {std::size_t(nparticles)}});
    engine.Put(var, data.data());
  };

  engine.BeginStep();

  io.DefineAttribute("cctk_iteration", &cctk_iteration, 1);
  io.DefineAttribute("cctk_time", &cctk_time, 1);

  output_variable("id", particles.id);
  output_variable("tau", particles.tau);
  output_variable("posx", particles.pos[0]);
  output_variable("posy", particles.pos[1]);
  output_variable("posz", particles.pos[2]);
  output_variable("velx", particles.vel[0]);
  output_variable("vely", particles.vel[1]);
  output_variable("velz", particles.vel[2]);
  output_variable("accx", particles.acc[0]);
  output_variable("accy", particles.acc[1]);
  output_variable("accz", particles.acc[2]);

  engine.PerformPuts();
  engine.EndStep();
  engine.Close();
}

} // namespace ParticlesX
