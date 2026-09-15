#include <defs.hxx>
#include <loop_device.hxx>
#include <vect.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <atomic>
#include <cmath>

namespace TestProlongateStar {
using namespace Arith;
using namespace std;

////////////////////////////////////////////////////////////////////////////////

// Let's test all polynomial powers in a single run! We choose one particular
// set of powers for each iteration.
void choose_powers(const int iteration, int &restrict px, int &restrict py,
                   int &restrict pz) {
  DECLARE_CCTK_PARAMETERS;
  px = iteration % (max_power_x + 1);
  py = iteration / (max_power_x + 1) % (max_power_y + 1);
  pz = iteration / (max_power_x + 1) / (max_power_y + 1) % (max_power_z + 1);
}

CCTK_DEVICE CCTK_HOST constexpr CCTK_REAL polyfun1d(const CCTK_REAL x,
                                                    const int px) {
  // Our vertex-centred test function is a polynomial in the coordinates
  return pown(x, px);
}

CCTK_DEVICE CCTK_HOST constexpr CCTK_REAL polyfun(const CCTK_REAL x,
                                                  const CCTK_REAL y,
                                                  const CCTK_REAL z,
                                                  const int px, const int py,
                                                  const int pz) {
  return polyfun1d(x, px) + polyfun1d(y, py) + polyfun1d(z, pz);
}

CCTK_DEVICE CCTK_HOST constexpr CCTK_REAL ipolyfun1d(const CCTK_REAL x,
                                                     const int px) {
  // Our cell-centred test function is the integral of a polynomial in the
  // coordinates
  return (px + 1) * pown(x, px + 1);
}

CCTK_DEVICE CCTK_HOST constexpr CCTK_REAL consfun1d(const CCTK_REAL x,
                                                    const CCTK_REAL dx,
                                                    const int px) {
  return (ipolyfun1d(x + dx / 2, px) - ipolyfun1d(x - dx / 2, px)) / dx;
}

CCTK_DEVICE CCTK_HOST constexpr CCTK_REAL
consfun(const CCTK_REAL x, const CCTK_REAL y, const CCTK_REAL z,
        const CCTK_REAL dx, const CCTK_REAL dy, const CCTK_REAL dz,
        const int px, const int py, const int pz) {
  return consfun1d(x, dx, px) + consfun1d(y, dy, py) + consfun1d(z, dz, pz);
}

CCTK_DEVICE CCTK_HOST constexpr CCTK_REAL
fluxxfun(const CCTK_REAL x, const CCTK_REAL y, const CCTK_REAL z,
         const CCTK_REAL dx, const CCTK_REAL dy, const CCTK_REAL dz,
         const int px, const int py, const int pz) {
  return polyfun1d(x, px) + consfun1d(y, dy, py) + consfun1d(z, dz, pz);
}

CCTK_DEVICE CCTK_HOST constexpr CCTK_REAL
fluxyfun(const CCTK_REAL x, const CCTK_REAL y, const CCTK_REAL z,
         const CCTK_REAL dx, const CCTK_REAL dy, const CCTK_REAL dz,
         const int px, const int py, const int pz) {
  return consfun1d(x, dx, px) + polyfun1d(y, py) + consfun1d(z, dz, pz);
}

CCTK_DEVICE CCTK_HOST constexpr CCTK_REAL
fluxzfun(const CCTK_REAL x, const CCTK_REAL y, const CCTK_REAL z,
         const CCTK_REAL dx, const CCTK_REAL dy, const CCTK_REAL dz,
         const int px, const int py, const int pz) {
  return consfun1d(x, dx, px) + consfun1d(y, dy, py) + polyfun1d(z, pz);
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void TestProlongateStar_Set(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestProlongateStar_Set;
  DECLARE_CCTK_PARAMETERS;

  if (cctk_component == 0 && cctk_tile_min[0] == 0 && cctk_tile_min[1] == 0 &&
      cctk_tile_min[2] == 0)
    CCTK_INFO("TestProlongateStar_Set");

  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);

  int px, py, pz;
  choose_powers(cctk_iteration, px, py, pz);

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones, [=] CCTK_DEVICE(const Loop::PointDesc &p) {
        const CCTK_REAL x = vcoordx(p.I);
        const CCTK_REAL y = vcoordy(p.I);
        const CCTK_REAL z = vcoordz(p.I);
        const CCTK_REAL goodval = polyfun(x, y, z, px, py, pz);
        alp(p.I) = goodval;
      });

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const Loop::PointDesc &p) {
        const CCTK_REAL x = ccoordx(p.I);
        const CCTK_REAL y = ccoordy(p.I);
        const CCTK_REAL z = ccoordz(p.I);
        const CCTK_REAL goodval = consfun(x, y, z, dx, dy, dz, px, py, pz);
        dens(p.I) = goodval;
      });

  grid.loop_int_device<0, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const Loop::PointDesc &p) {
        const CCTK_REAL x = p.X[0];
        const CCTK_REAL y = p.X[1];
        const CCTK_REAL z = p.X[2];
        const CCTK_REAL goodval = fluxxfun(x, y, z, dx, dy, dz, px, py, pz);
        flux_x(p.I) = goodval;
      });

  grid.loop_int_device<1, 0, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const Loop::PointDesc &p) {
        const CCTK_REAL x = p.X[0];
        const CCTK_REAL y = p.X[1];
        const CCTK_REAL z = p.X[2];
        const CCTK_REAL goodval = fluxyfun(x, y, z, dx, dy, dz, px, py, pz);
        flux_y(p.I) = goodval;
      });

  grid.loop_int_device<1, 1, 0>(
      grid.nghostzones, [=] CCTK_DEVICE(const Loop::PointDesc &p) {
        const CCTK_REAL x = p.X[0];
        const CCTK_REAL y = p.X[1];
        const CCTK_REAL z = p.X[2];
        const CCTK_REAL goodval = fluxzfun(x, y, z, dx, dy, dz, px, py, pz);
        flux_z(p.I) = goodval;
      });
}

////////////////////////////////////////////////////////////////////////////////

std::atomic<bool> found_error;

extern "C" void TestProlongateStar_CheckBegin(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestProlongateStar_CheckBegin;
  DECLARE_CCTK_PARAMETERS;

  found_error = false;
}

extern "C" void TestProlongateStar_Check(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestProlongateStar_Check;
  DECLARE_CCTK_PARAMETERS;

  using std::abs, std::max;

  if (cctk_component == 0 && cctk_tile_min[0] == 0 && cctk_tile_min[1] == 0 &&
      cctk_tile_min[2] == 0)
    CCTK_INFO("TestProlongateStar_Check");

  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);

  int px, py, pz;
  choose_powers(cctk_iteration, px, py, pz);

  const CCTK_REAL maxrelerr = max({px, py, pz}) < 9 ? 0 : 1.0e-15;

  grid.loop_all<0, 0, 0>(grid.nghostzones, [&](const Loop::PointDesc &p) {
    // Skip outer boundary points
    if (maxabs(p.NI) != 0)
      return;

    const CCTK_REAL x = vcoordx(p.I);
    const CCTK_REAL y = vcoordy(p.I);
    const CCTK_REAL z = vcoordz(p.I);
    const CCTK_REAL goodval = polyfun(x, y, z, px, py, pz);

    if (abs(alp(p.I) - goodval) > maxrelerr * abs(goodval)) {
      CCTK_VERROR("alp px=%d py=%d pz=%d expected=%.17g, result=%.17g", px, py,
                  pz, double(goodval), double(alp(p.I)));
      found_error = true;
    }
  });

  // Skip testing too-high prolongation orders for conserved variables
  if (max({px, py, pz}) <= dens_prolongation_order) {
    grid.loop_all<1, 1, 1>(grid.nghostzones, [&](const Loop::PointDesc &p) {
      // Skip outer boundary points
      if (maxabs(p.NI) != 0)
        return;

      const CCTK_REAL x = ccoordx(p.I);
      const CCTK_REAL y = ccoordy(p.I);
      const CCTK_REAL z = ccoordz(p.I);
      const CCTK_REAL goodval = consfun(x, y, z, dx, dy, dz, px, py, pz);

      if (abs(dens(p.I) - goodval) > maxrelerr * abs(goodval)) {
        CCTK_VERROR("dens p=[%d,%d,%d] i=[%d,%d,%d] x=[%.17g,%.17g,%.17g] "
                    "expected=%.17g result=%.17g abserr=%.17g relerr=%.17g",
                    px, py, pz, p.I[0], p.I[1], p.I[2], double(p.X[0]),
                    double(p.X[1]), double(p.X[2]), double(goodval),
                    double(dens(p.I)), double(abs(dens(p.I) - goodval)),
                    double(abs((dens(p.I) - goodval)) / goodval));
        found_error = true;
      }
    });
  }

  // Skip testing too-high prolongation orders for flux variables
  if (max({px, py, pz}) <= flux_prolongation_order) {
    grid.loop_all<0, 1, 1>(grid.nghostzones, [&](const Loop::PointDesc &p) {
      // Skip outer boundary points
      if (maxabs(p.NI) != 0)
        return;

      const CCTK_REAL x = p.X[0];
      const CCTK_REAL y = p.X[1];
      const CCTK_REAL z = p.X[2];
      const CCTK_REAL goodval = fluxxfun(x, y, z, dx, dy, dz, px, py, pz);

      if (abs(flux_x(p.I) - goodval) > maxrelerr * abs(goodval)) {
        CCTK_VERROR("flux_x p=[%d,%d,%d] i=[%d,%d,%d] x=[%.17g,%.17g,%.17g] "
                    "expected=%.17g result=%.17g abserr=%.17g relerr=%.17g",
                    px, py, pz, p.I[0], p.I[1], p.I[2], double(p.X[0]),
                    double(p.X[1]), double(p.X[2]), double(goodval),
                    double(flux_x(p.I)), double(abs(flux_x(p.I) - goodval)),
                    double(abs((flux_x(p.I) - goodval)) / goodval));
        found_error = true;
      }
    });

    grid.loop_all<1, 0, 1>(grid.nghostzones, [&](const Loop::PointDesc &p) {
      // Skip outer boundary points
      if (maxabs(p.NI) != 0)
        return;

      const CCTK_REAL x = p.X[0];
      const CCTK_REAL y = p.X[1];
      const CCTK_REAL z = p.X[2];
      const CCTK_REAL goodval = fluxyfun(x, y, z, dx, dy, dz, px, py, pz);

      if (abs(flux_y(p.I) - goodval) > maxrelerr * abs(goodval)) {
        CCTK_VERROR("flux_y px=%d py=%d pz=%d expected=%.17g, result=%.17g", px,
                    py, pz, double(goodval), double(flux_y(p.I)));
        found_error = true;
      }
    });

    grid.loop_all<1, 1, 0>(grid.nghostzones, [&](const Loop::PointDesc &p) {
      // Skip outer boundary points
      if (maxabs(p.NI) != 0)
        return;

      const CCTK_REAL x = p.X[0];
      const CCTK_REAL y = p.X[1];
      const CCTK_REAL z = p.X[2];
      const CCTK_REAL goodval = fluxzfun(x, y, z, dx, dy, dz, px, py, pz);

      if (abs(flux_z(p.I) - goodval) > maxrelerr * abs(goodval)) {
        CCTK_VERROR("flux_z px=%d py=%d pz=%d expected=%.17g, result=%.17g", px,
                    py, pz, double(goodval), double(flux_z(p.I)));
        found_error = true;
      }
    });
  }
}

extern "C" void TestProlongateStar_CheckEnd(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestProlongateStar_CheckEnd;
  DECLARE_CCTK_PARAMETERS;

  if (found_error)
    CCTK_ERROR("Found error");
}

} // namespace TestProlongateStar
