#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iomanip>
#include <sstream>
#include <string>
#include <map>
#include <sys/stat.h>
#include <errno.h>
#include <iostream>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#ifdef HAVE_CAPABILITY_HDF5
// We currently support the HDF5 1.6 API (and when using 1.8 the
// compatibility mode introduced by H5_USE_16_API).  Several machines
// in SimFactory use HDF5 1.6, so we cannot drop support for it.  It
// seems it is hard to support both the 1.6 and 1.8 API
// simultaneously; for example H5Fopen takes a different number of
// arguments in the two versions.
#define H5_USE_16_API
#include <hdf5.h>
#endif

#include "utils.hh"
#include "integrate.hh"
#include "multipole.hh"

// check return code of HDF5 call abort with an error message if there was an error.
// adapted from CarpetIOHDF5/src/CarpetIOHDF5.hh.
#define HDF5_ERROR(fn_call)                                                   \
        do {                                                                  \
          hid_t _error_code = fn_call;                                        \
                                                                              \
                                                                              \
          if (_error_code < 0)                                                \
          {                                                                   \
            CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,              \
                        "HDF5 call '%s' returned error code %d",              \
                        #fn_call, (int)_error_code);                          \
          }                                                                   \
        } while (0)


using namespace std;

////////////////////////////////////////////////////////////////////////
// File manipulation
////////////////////////////////////////////////////////////////////////

FILE *Multipole_OpenOutputFile(CCTK_ARGUMENTS, const string &name)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  bool first_time = cctk_iteration == 0;
  const char *mode = first_time ? "w" : "a";
  const char *my_out_dir = strcmp(out_dir, "") ? out_dir : io_out_dir;
  const int err = CCTK_CreateDirectory(0755, my_out_dir);
  if (err < 0)
    CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING, 
               "Multipole output directory %s could not be created (error code %d)",
               my_out_dir, err);

  string output_name(string(my_out_dir) + "/" + name);

  FILE *fp = fopen(output_name.c_str(), mode);

  if (fp == 0)
  {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING, "%s", (string("Could not open output file ") + output_name).c_str());
  }

  return fp;
}

////////////////////////////////////////////////////////////////////////
// Unused
////////////////////////////////////////////////////////////////////////

void Multipole_OutputArray(CCTK_ARGUMENTS, FILE *f, int array_size, 
                           CCTK_REAL const th[], CCTK_REAL const ph[], 
                           CCTK_REAL const xs[], CCTK_REAL const ys[], CCTK_REAL const zs[],
                           CCTK_REAL const data[])
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  CCTK_REAL last_ph = ph[0];

  for (int i = 0; i < array_size; i++)
  {
    if (ph[i] != last_ph) // Separate blocks for gnuplot
      fprintf(f, "\n");
    fprintf(f, "%f %f %f %f %f %f %.19g\n", cctk_time, th[i], ph[i], xs[i], ys[i], zs[i], data[i]);
    last_ph = ph[i];
  }
}


void Multipole_OutputArrayToFile(CCTK_ARGUMENTS, const string &name, int array_size,
                                 CCTK_REAL const th[], CCTK_REAL const ph[],
                                 CCTK_REAL const xs[], CCTK_REAL const ys[], CCTK_REAL const zs[],
                                 CCTK_REAL const data[])
{
  DECLARE_CCTK_ARGUMENTS;

  if (FILE *fp = Multipole_OpenOutputFile(CCTK_PASS_CTOC, name))
  {
    Multipole_OutputArray(CCTK_PASS_CTOC, fp, array_size, th, ph, xs, ys, zs, data);
    fclose(fp);
  }
}

////////////////////////////////////////////////////////////////////////
// Misc
////////////////////////////////////////////////////////////////////////

void Multipole_Output1D(CCTK_ARGUMENTS, const string &name, int array_size,
                        CCTK_REAL const th[], CCTK_REAL const ph[], mp_coord coord,
                        CCTK_REAL const data[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (FILE *f = Multipole_OpenOutputFile(CCTK_PASS_CTOC, name))
  {
    fprintf(f, "\"Time = %.19g\n", cctk_time);

    if (coord == mp_theta)
    {
      for (int i = 0; i <= ntheta; i++)
      {
        int idx = Multipole_Index(i, 0, ntheta);
        fprintf(f, "%f %.19g\n", th[idx], data[idx]);
      }
    }
    else if (coord == mp_phi)
    {
      for (int i = 0; i <= nphi; i++)
      {
        int idx = Multipole_Index(ntheta / 4, i, ntheta);
        fprintf(f, "%f %.19g\n", ph[idx], data[idx]);
      }
    }
    fprintf(f, "\n\n");
    fclose(f);
  }
}

////////////////////////////////////////////////////////////////////////
// Complex IO
////////////////////////////////////////////////////////////////////////

void Multipole_OutputComplex(CCTK_ARGUMENTS, FILE *fp, CCTK_REAL redata, CCTK_REAL imdata)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;
  fprintf(fp, "%f %.19g %.19g\n", cctk_time, redata, imdata);
}

void Multipole_OutputComplexToFile(CCTK_ARGUMENTS, const string &name, CCTK_REAL redata, CCTK_REAL imdata)
{
  DECLARE_CCTK_ARGUMENTS;

  if (FILE *fp = Multipole_OpenOutputFile(CCTK_PASS_CTOC, name))
  {
    Multipole_OutputComplex(CCTK_PASS_CTOC, fp, redata, imdata);
    fclose(fp);
  }
}

////////////////////////////////////////////////////////////////////////
// HDF5 complex output
////////////////////////////////////////////////////////////////////////

#ifdef HAVE_CAPABILITY_HDF5

static bool file_exists(const string &name)
{
  struct stat sts;
  return !(stat(name.c_str(), &sts) == -1 && errno == ENOENT);
}

static bool dataset_exists(hid_t file, const string &dataset_name)
{
  // To test whether a dataset exists, the recommended way in API 1.6
  // is to use H5Gget_objinfo, but this prints an error to stderr if
  // the dataset does not exist.  We explicitly avoid this by wrapping
  // the call in H5E_BEGIN_TRY/H5E_END_TRY statements.  In 1.8,
  // H5Gget_objinfo is deprecated, and H5Lexists does the job.  See
  // http://www.mail-archive.com/hdf-forum@hdfgroup.org/msg00125.html

  #if 1
  bool exists;
  H5E_BEGIN_TRY
  {
    exists = H5Gget_objinfo(file, dataset_name.c_str(), 1, NULL) >= 0;
  }
  H5E_END_TRY;
  return exists;
  #else
  return H5Lexists(file, dataset_name.c_str(), H5P_DEFAULT);
  #endif
}

void Multipole_OutputComplexToH5File(CCTK_ARGUMENTS, const Multipole::variable_desc vars[],
                                     const CCTK_REAL radii[], const Multipole::mode_array& modes)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const char *my_out_dir = strcmp(out_dir, "") ? out_dir : io_out_dir;
  if (CCTK_CreateDirectory(0755, my_out_dir) < 0)
    CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING, 
               "Multipole output directory %s could not be created",
               my_out_dir);

  static map<string,bool> checked; // Has the given file been checked
                                   // for truncation? map<*,bool>
                                   // defaults to false
  for (int v = 0; v < modes.get_nvars(); v++)
  {
    string basename = "mp_" + vars[v].name + ".h5";
    string output_name = my_out_dir + string("/") + basename;

    hid_t file;

    if (!file_exists(output_name) || (!checked[output_name] && IO_TruncateOutputFiles(cctkGH)))
    {
      file = H5Fcreate(output_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }
    else
    {
      file = H5Fopen(output_name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    }

    checked[output_name] = true;

    for(int i=0; i<modes.get_nradii(); i++)
    {
      const CCTK_REAL rad = radii[i];
      for (int l=0; l <= modes.get_lmax(); l++)
      {
        for (int m=-l; m <= l; m++)
        {
          ostringstream datasetname;
          datasetname << "l" << l << "_m" << m << "_r" << setiosflags(ios::fixed) << setprecision(2) << rad;

          hid_t dataset = -1;

          if (dataset_exists(file, datasetname.str()))
          {
            dataset = H5Dopen(file, datasetname.str().c_str());
          }
          else
          {
            hsize_t dims[2] = {0,3};
            hsize_t maxdims[2] = {H5S_UNLIMITED, 3};
            hid_t dataspace = H5Screate_simple(2, dims, maxdims);

            hid_t cparms = -1;
            hsize_t chunk_dims[2] = {hsize_t(hdf5_chunk_size),3};
            cparms = H5Pcreate (H5P_DATASET_CREATE);
            HDF5_ERROR(H5Pset_chunk(cparms, 2, chunk_dims));

            dataset = H5Dcreate(file, datasetname.str().c_str(), H5T_NATIVE_DOUBLE, dataspace, cparms);
            H5Pclose(cparms);
          }

          hid_t filespace = H5Dget_space (dataset);
          hsize_t filedims[2];
          hsize_t maxdims[2];
          HDF5_ERROR(H5Sget_simple_extent_dims(filespace, filedims, maxdims));

          filedims[0] += 1;
          hsize_t size[2] = {filedims[0], filedims[1]};
          HDF5_ERROR(H5Dextend (dataset, size));
          HDF5_ERROR(H5Sclose(filespace));

          /* Select a hyperslab  */
          hsize_t offset[2] = {filedims[0]-1, 0};
          hsize_t dims2[2] = {1,3};
          filespace = H5Dget_space (dataset);
          HDF5_ERROR(H5Sselect_hyperslab (filespace, H5S_SELECT_SET, offset, NULL,
                                          dims2, NULL));

          CCTK_REAL data[] = {cctk_time, modes(v, i, l, m, 0), modes(v, i, l, m, 1)};

          hid_t memdataspace = H5Screate_simple(2, dims2, NULL);

          /* Write the data to the hyperslab  */
          HDF5_ERROR(H5Dwrite (dataset, H5T_NATIVE_DOUBLE, memdataspace, filespace,
                               H5P_DEFAULT, data));

          HDF5_ERROR(H5Dclose(dataset));
          HDF5_ERROR(H5Sclose(filespace));
          HDF5_ERROR(H5Sclose(memdataspace));
        }
      }
    }
    HDF5_ERROR(H5Fclose(file));
  }
}

#else

void Multipole_OutputComplexToH5File(CCTK_ARGUMENTS, const Multipole::variable_desc vars[],
                                     const CCTK_REAL radii[], const Multipole::mode_array& modes)
{
  CCTK_WARN(0,"HDF5 output has been requested but Cactus has been compiled without HDF5 support");
}

#endif

////////////////////////////////////////////////////////////////////////
// Coordinates
////////////////////////////////////////////////////////////////////////

void Multipole_CoordSetup(CCTK_REAL xhat[], CCTK_REAL yhat[],
                          CCTK_REAL zhat[], CCTK_REAL th[],
                          CCTK_REAL ph[])
{
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL PI = acos(-1.0);

  // Add an offset for midpoint integration.
  CCTK_REAL is_midpoint = 0.0;
  if (CCTK_Equals(integration_method, "midpoint")) {
    is_midpoint = 1.0;
  }
  const CCTK_REAL dth = PI/(ntheta+is_midpoint);
  const CCTK_REAL dph = 2*PI/(nphi+is_midpoint);
  for (int it = 0; it <= ntheta; it++)
    {
      for (int ip = 0; ip <= nphi; ip++)
        {
          const int i = Multipole_Index(it, ip, ntheta);
          // Check for when midpoint enabled:
          //   dth = PI/(ntheta+1) -> ntheta = PI/dth - 1
          //   th[i] = i * dth + 0.5*dth
          //  Therefore:
          //   th[0]      = 0.5*dth      <- GOOD.
          //   th[ntheta] = ntheta*dth + 0.5*dth
          //              = ((PI/dth)-1)*dth + 0.5*dth
          //              = PI - dth + 0.5*dth
          //              = PI - 0.5*dth <- GOOD.
          //  Similarly for ph.

          // Check for when midpoint disabled:
          //   dth = PI/ntheta -> ntheta = PI/dth
          //   th[i] = i * dth
          //  Therefore:
          //   th[0]      = 0      <- GOOD.
          //   th[ntheta] = ntheta*dth
          //              = PI/dth*dth
          //              = PI     <- GOOD.
          // Similarly for ph.
          th[i] = it * dth + 0.5*dth*is_midpoint;
          ph[i] = ip * dph + 0.5*dph*is_midpoint;
          xhat[i] = cos(ph[i])*sin(th[i]);
          yhat[i] = sin(ph[i])*sin(th[i]);
          zhat[i] = cos(th[i]);
        }
    }
}

void Multipole_ScaleCartesian(int ntheta, int nphi, CCTK_REAL r,
                              CCTK_REAL const xhat[], CCTK_REAL const yhat[], CCTK_REAL const zhat[],
                              CCTK_REAL x[], CCTK_REAL y[], CCTK_REAL z[])
{
  for (int it = 0; it <= ntheta; it++)
  {
    for (int ip = 0; ip <= nphi; ip++)
    {
      const int i = Multipole_Index(it, ip, ntheta);

      x[i] = r * xhat[i];
      y[i] = r * yhat[i];
      z[i] = r * zhat[i];
    }
  }
}

////////////////////////////////////////////////////////////////////////
// Integration
////////////////////////////////////////////////////////////////////////

void Multipole_Integrate(int array_size, int nthetap,
    CCTK_REAL const array1r[], CCTK_REAL const array1i[],
    CCTK_REAL const array2r[], CCTK_REAL const array2i[],
    CCTK_REAL const th[], CCTK_REAL const ph[], 
    CCTK_REAL *outre, CCTK_REAL *outim)
{
  DECLARE_CCTK_PARAMETERS

  int il = Multipole_Index(0,0,ntheta);
  int iu = Multipole_Index(1,0,ntheta);
  CCTK_REAL dth = th[iu] - th[il];
  iu = Multipole_Index(0,1,ntheta);
  CCTK_REAL dph = ph[iu] - ph[il];

  static CCTK_REAL *fr = 0;
  static CCTK_REAL *fi = 0;
  static bool allocated_memory = false;

  // Construct an array for the real integrand
  if (!allocated_memory)
  {
    fr = new CCTK_REAL[array_size];
    fi = new CCTK_REAL[array_size];
    allocated_memory = true;
  }
    
  // the below calculations take the integral of conj(array1)*array2*sin(th)
  for (int i = 0; i < array_size; i++)
  {
    fr[i] = (array1r[i] * array2r[i] + 
             array1i[i] * array2i[i] ) * sin(th[i]);
    fi[i] = (array1r[i] * array2i[i] - 
             array1i[i] * array2r[i] ) * sin(th[i]);
  }
    
  if (CCTK_Equals(integration_method, "midpoint"))
  {
    *outre = Midpoint2DIntegral(fr, ntheta, nphi, dth, dph);
    *outim = Midpoint2DIntegral(fi, ntheta, nphi, dth, dph);
  }
  else if (CCTK_Equals(integration_method, "trapezoidal"))
  {
    *outre = Trapezoidal2DIntegral(fr, ntheta, nphi, dth, dph);
    *outim = Trapezoidal2DIntegral(fi, ntheta, nphi, dth, dph);
  }
  else if (CCTK_Equals(integration_method, "Simpson"))
  {
    if (nphi % 2 != 0 || ntheta % 2 != 0)
    {
      CCTK_WARN (CCTK_WARN_ABORT, "The Simpson integration method requires even ntheta and even nphi");
    }
    *outre = Simpson2DIntegral(fr, ntheta, nphi, dth, dph);
    *outim = Simpson2DIntegral(fi, ntheta, nphi, dth, dph);
  }
  else if (CCTK_Equals(integration_method, "DriscollHealy"))
  {
    if (ntheta % 2 != 0)
    {
      CCTK_WARN (CCTK_WARN_ABORT, "The Driscoll&Healy integration method requires even ntheta");
    }
    *outre = DriscollHealy2DIntegral(fr, ntheta, nphi, dth, dph);
    *outim = DriscollHealy2DIntegral(fi, ntheta, nphi, dth, dph);
  }
  else
  {
    CCTK_WARN(CCTK_WARN_ABORT, "internal error");
  }
} 
