#include "solve.hxx"

namespace ODESolvers {
using namespace std;

////////////////////////////////////////////////////////////////////////////////

// Initialize the temporary mfab mechanism
void statecomp_t::init_tmp_mfabs() {
  assert(CarpetX::active_levels);
  CarpetX::active_levels->loop_serially([&](const auto &leveldata) {
    for (const auto &groupdataptr : leveldata.groupdata) {
      if (groupdataptr == nullptr)
        continue;
      const auto &groupdata = *groupdataptr;
      groupdata.init_tmp_mfabs();
    }
  });
}

// Free all temporary mfabs that we might have allocated
void statecomp_t::free_tmp_mfabs() {
  assert(CarpetX::active_levels);
  CarpetX::active_levels->loop_serially([&](const auto &leveldata) {
    for (const auto &groupdataptr : leveldata.groupdata) {
      if (groupdataptr == nullptr)
        continue;
      const auto &groupdata = *groupdataptr;
      groupdata.free_tmp_mfabs();
    }
  });
}

// State that the state vector has valid data in the interior
void statecomp_t::set_valid(const valid_t valid) const {
  const int tl = this->timelevel;
  for (auto groupdata : groupdatas) {
    for (int vi = 0; vi < groupdata->numvars; ++vi) {
      groupdata->valid.at(tl).at(vi).set_int(valid.valid_int, [=]() {
        std::ostringstream buf;
        buf << "ODESolvers after lincomb: Mark interior as "
            << (valid.valid_int ? "valid" : "invalid");
        return buf.str();
      });
      groupdata->valid.at(tl).at(vi).set_outer(valid.valid_outer, [=]() {
        std::ostringstream buf;
        buf << "ODESolvers after lincomb: Mark outer boundary as "
            << (valid.valid_outer ? "valid" : "invalid");
        return buf.str();
      });
      groupdata->valid.at(tl).at(vi).set_ghosts(valid.valid_ghosts, [=]() {
        std::ostringstream buf;
        buf << "ODESolvers after lincomb: Mark ghosts as "
            << (valid.valid_int ? "valid" : "invalid");
        return buf.str();
      });
      // // TODO: Parallelize over patches, levels, group, variables, and
      // // timelevels
      // const active_levels_t active_levels(
      //     groupdata->level, groupdata->level + 1, groupdata->patch,
      //     groupdata->patch + 1);
      // CarpetX::poison_invalid_gf(active_levels, groupdata->groupindex, vi,
      // tl);
    }
  }
}

// Combine validity information from several sources
template <std::size_t N>
void statecomp_t::combine_valids(const statecomp_t &dst, const CCTK_REAL scale,
                                 const std::array<CCTK_REAL, N> &factors,
                                 const std::array<const statecomp_t *, N> &srcs,
                                 const CarpetX::valid_t where) {
  const int ngroups = dst.groupdatas.size();
  for (const auto &src : srcs)
    assert(int(src->groupdatas.size()) == ngroups);
  for (int group = 0; group < ngroups; ++group) {
    const auto &dstgroup = dst.groupdatas.at(group);
    const int nvars = dstgroup->numvars;
    for (const auto &src : srcs) {
      const auto &srcgroup = src->groupdatas.at(group);
      assert(srcgroup->numvars == nvars);
    }
  }

  const int dst_tl = dst.timelevel;
  for (int group = 0; group < ngroups; ++group) {
    const auto &dstgroup = dst.groupdatas.at(group);
    const int nvars = dstgroup->numvars;
    for (int vi = 0; vi < nvars; ++vi) {
      CarpetX::valid_t valid = where;
      bool did_set_valid = false;
      if (scale != 0) {
        valid &= dstgroup->valid.at(dst_tl).at(vi).get();
        did_set_valid = true;
      }
      for (std::size_t m = 0; m < srcs.size(); ++m) {
        if (factors.at(m) != 0) {
          const auto &src = srcs.at(m);
          const auto &srcgroup = src->groupdatas.at(group);
          valid &= srcgroup->valid.at(src->timelevel).at(vi).get();
          did_set_valid = true;
        }
      }
      if (!did_set_valid)
        valid = valid_t(false);
      dstgroup->valid.at(dst_tl).at(vi) =
          why_valid_t(valid, []() { return "Set from RHS in ODESolvers"; });
    }
  }
}

// Ensure a state vector has valid data in the interior
void statecomp_t::check_valid(const valid_t required,
                              const function<string()> &why) const {
  const int tl = this->timelevel;
  for (const auto groupdata : groupdatas) {
    for (int vi = 0; vi < groupdata->numvars; ++vi) {
      CarpetX::error_if_invalid(*groupdata, vi, tl, required, why);
      // TODO: Parallelize over pathces, levels, group, variables, and
      // timelevels
      const CarpetX::active_levels_t active_levels(
          groupdata->level, groupdata->level + 1, groupdata->patch,
          groupdata->patch + 1);
      CarpetX::check_valid_gf(active_levels, groupdata->groupindex, vi, tl,
                              CarpetX::nan_handling_t::forbid_nans, why);
    }
  }
}

// Copy state vector into newly allocated memory
statecomp_t statecomp_t::copy(const CarpetX::valid_t where) const {
  const std::size_t size = mfabs.size();
  statecomp_t result;
  result.timelevel = this->timelevel;
  result.groupdatas.reserve(size);
  result.mfabs.reserve(size);
  for (std::size_t n = 0; n < size; ++n) {
    const auto groupdata = groupdatas.at(n);
    // This global nan-check doesn't work since we don't care about the
    // boundaries
    // #ifdef CCTK_DEBUG
    //     const auto &x = mfabs.at(n);
    //     if (x->contains_nan())
    //       CCTK_VERROR("statecomp_t::copy.x: Group %s contains nans",
    //                   groupdata->groupname.c_str());
    // #endif
    auto y = groupdata->alloc_tmp_mfab();
    result.groupdatas.push_back(groupdata);
    result.mfabs.push_back(y);
  }
  lincomb(result, 0, make_array(CCTK_REAL(1)), make_array(this), where);
  // This global nan-check doesn't work since we don't care about the boundaries
  // #ifdef CCTK_DEBUG
  //   for (std::size_t n = 0; n < size; ++n) {
  //     const auto groupdata = result.groupdatas.at(n);
  //     const auto &y = result.mfabs.at(n);
  //     if (y->contains_nan())
  //       CCTK_VERROR("statecomp_t::copy.y: Group %s contains nans",
  //                   groupdata->groupname.c_str());
  //   }
  // #endif
  return result;
}

template <std::size_t N>
void statecomp_t::lincomb(const statecomp_t &dst, const CCTK_REAL scale,
                          const std::array<CCTK_REAL, N> &factors,
                          const std::array<const statecomp_t *, N> &srcs,
                          const CarpetX::valid_t where) {
  const std::size_t size = dst.mfabs.size();
  for (std::size_t n = 0; n < N; ++n)
    assert(srcs[n]->mfabs.size() == size);
  for (std::size_t m = 0; m < size; ++m) {
    const auto ncomp = dst.mfabs.at(m)->nComp();
    const auto ngrowvect = dst.mfabs.at(m)->nGrowVect();
    for (std::size_t n = 0; n < N; ++n) {
      assert(srcs[n]->mfabs.at(m)->nComp() == ncomp);
      assert(srcs[n]->mfabs.at(m)->nGrowVect() == ngrowvect);
    }
  }

  using std::isfinite;
  assert(isfinite(scale));
  const bool read_dst = scale != 0;
  for (std::size_t n = 0; n < N; ++n)
    assert(isfinite(factors[n]));

  statecomp_t::combine_valids(dst, scale, factors, srcs, where);

#ifndef AMREX_USE_GPU
  std::vector<std::function<void()> > tasks;
#endif

  // TODO: Poison ghosts/boundaries

  for (std::size_t m = 0; m < size; ++m) {
    const std::ptrdiff_t ncomps = dst.mfabs.at(m)->nComp();
    const auto mfitinfo = amrex::MFItInfo().DisableDeviceSync();
    for (amrex::MFIter mfi(*dst.mfabs.at(m), mfitinfo); mfi.isValid(); ++mfi) {
      const amrex::Array4<CCTK_REAL> dstvar = dst.mfabs.at(m)->array(mfi);
      std::array<amrex::Array4<const CCTK_REAL>, N> srcvars;
      for (std::size_t n = 0; n < N; ++n)
        srcvars[n] = srcs[n]->mfabs.at(m)->const_array(mfi);
      for (std::size_t n = 0; n < N; ++n) {
        assert(srcvars[n].jstride == dstvar.jstride);
        assert(srcvars[n].kstride == dstvar.kstride);
        assert(srcvars[n].nstride == dstvar.nstride);
      }
      const std::ptrdiff_t nstride = dstvar.nstride;
      const std::ptrdiff_t npoints = nstride * ncomps;

      CCTK_REAL *restrict const dstptr = dstvar.dataPtr();
      std::array<const CCTK_REAL *restrict, N> srcptrs;
      for (std::size_t n = 0; n < N; ++n)
        srcptrs[n] = srcvars[n].dataPtr();

#ifndef AMREX_USE_GPU
      // CPU

      const std::ptrdiff_t ntiles = omp_get_max_threads();
      const std::ptrdiff_t tile_size = Arith::align_ceil(
          Arith::div_ceil(npoints, ntiles), std::ptrdiff_t(64));

      for (std::ptrdiff_t imin = 0; imin < npoints; imin += tile_size) {
        using std::min;
        const std::ptrdiff_t imax = min(npoints, imin + tile_size);

        if (!read_dst && N == 1 && factors[0] == 1) {
          // Copy

          auto task = [=]() {
            std::memcpy(&dstptr[imin], &srcptrs[0][imin],
                        (imax - imin) * sizeof *dstptr);
          };
          tasks.push_back(std::move(task));

        } else if (!read_dst && N >= 1 && factors[0] == 1) {
          // Write without scaling

          auto task = [=]() {
#pragma omp simd
            for (std::ptrdiff_t i = imin; i < imax; ++i) {
              CCTK_REAL accum = srcptrs[0][i];
              for (std::size_t n = 1; n < N; ++n)
                accum += factors[n] * srcptrs[n][i];
              dstptr[i] = accum;
            }
          };
          tasks.push_back(std::move(task));

        } else if (!read_dst) {
          // Write

          auto task = [=]() {
#pragma omp simd
            for (std::ptrdiff_t i = imin; i < imax; ++i) {
              CCTK_REAL accum = 0;
              for (std::size_t n = 0; n < N; ++n)
                accum += factors[n] * srcptrs[n][i];
              dstptr[i] = accum;
            }
          };
          tasks.push_back(std::move(task));

        } else if (scale == 1) {
          // Update without scaling

          auto task = [=]() {
#pragma omp simd
            for (std::ptrdiff_t i = imin; i < imax; ++i) {
              CCTK_REAL accum = dstptr[i];
              for (std::size_t n = 0; n < N; ++n)
                accum += factors[n] * srcptrs[n][i];
              dstptr[i] = accum;
            }
          };
          tasks.push_back(std::move(task));

        } else {
          // Update

          auto task = [=]() {
#pragma omp simd
            for (std::ptrdiff_t i = imin; i < imax; ++i) {
              CCTK_REAL accum = scale * dstptr[i];
              for (std::size_t n = 0; n < N; ++n)
                accum += factors[n] * srcptrs[n][i];
              dstptr[i] = accum;
            }
          };
          tasks.push_back(std::move(task));
        }
      } // for imin

#else
      // GPU

      const CCTK_REAL scale1 = scale;
      assert(npoints < INT_MAX);
      const amrex::Box box(
          amrex::IntVect(0, 0, 0), amrex::IntVect(npoints - 1, 0, 0),
          amrex::IntVect(amrex::IndexType::CELL, amrex::IndexType::CELL,
                         amrex::IndexType::CELL));

      if (!read_dst) {

        amrex::launch(
            box, [=] CCTK_DEVICE(const amrex::Box &box) __attribute__((
                     __always_inline__, __flatten__)) {
              const int i = box.smallEnd()[0];
              // const int j = box.smallEnd()[1];
              // const int k = box.smallEnd()[2];
              CCTK_REAL accum = 0;
              // The ROCM 6.2 compiler can't handle
              // `std::array::operator[]`, so we avoid it via pointers:
              // for (std::size_t n = 0; n < N; ++n)
              //   accum += factors[n] * srcptrs[n][i];
              const CCTK_REAL *restrict const factors_ptr = factors.data();
              const CCTK_REAL *restrict const *restrict const srcptrs_ptr =
                  srcptrs.data();
              for (std::size_t n = 0; n < N; ++n)
                accum += factors_ptr[n] * srcptrs_ptr[n][i];
              dstptr[i] = accum;
            });

      } else {

        amrex::launch(
            box, [=] CCTK_DEVICE(const amrex::Box &box) __attribute__((
                     __always_inline__, __flatten__)) {
              const int i = box.smallEnd()[0];
              // const int j = box.smallEnd()[1];
              // const int k = box.smallEnd()[2];
              CCTK_REAL accum = scale1 * dstptr[i];
              // The ROCM 6.2 compiler can't handle
              // `std::array::operator[]`, so we avoid it via pointers:
              // for (std::size_t n = 0; n < N; ++n)
              //   accum += factors[n] * srcptrs[n][i];
              const CCTK_REAL *restrict const factors_ptr = factors.data();
              const CCTK_REAL *restrict const *restrict const srcptrs_ptr =
                  srcptrs.data();
              for (std::size_t n = 0; n < N; ++n)
                accum += factors_ptr[n] * srcptrs_ptr[n][i];
              dstptr[i] = accum;
            });
      }

#endif
    }
  }

#ifndef AMREX_USE_GPU
  // run all tasks
#pragma omp parallel for schedule(dynamic)
  for (std::size_t i = 0; i < tasks.size(); ++i)
    tasks[i]();
#else
  // wait for all tasks
  amrex::Gpu::synchronize();
  AMREX_GPU_ERROR_CHECK();
#endif
}

namespace detail {
template <std::size_t N>
void call_lincomb(const statecomp_t &dst, const CCTK_REAL scale,
                  const std::vector<CCTK_REAL> &factors,
                  const std::vector<const statecomp_t *> &srcs,
                  const std::vector<std::size_t> &indices,
                  const CarpetX::valid_t where) {
  assert(indices.size() == N);
  std::array<CCTK_REAL, N> factors1;
  std::array<const statecomp_t *, N> srcs1;
  for (std::size_t n = 0; n < N; ++n) {
    factors1[n] = factors.at(indices[n]);
    srcs1[n] = srcs.at(indices[n]);
  }
  statecomp_t::lincomb(dst, scale, factors1, srcs1, where);
}
} // namespace detail

void statecomp_t::lincomb(const statecomp_t &dst, const CCTK_REAL scale,
                          const std::vector<CCTK_REAL> &factors,
                          const std::vector<const statecomp_t *> &srcs,
                          const CarpetX::valid_t where) {
  const std::size_t N = factors.size();
  assert(srcs.size() == N);

  std::size_t NNZ = 0;
  for (std::size_t n = 0; n < N; ++n)
    NNZ += factors[n] != 0;
  std::vector<std::size_t> indices;
  indices.reserve(NNZ);
  for (std::size_t n = 0; n < N; ++n)
    if (factors[n] != 0)
      indices.push_back(n);
  assert(indices.size() == NNZ);

  switch (NNZ) {
  case 0:
    return detail::call_lincomb<0>(dst, scale, factors, srcs, indices, where);
  case 1:
    return detail::call_lincomb<1>(dst, scale, factors, srcs, indices, where);
  case 2:
    return detail::call_lincomb<2>(dst, scale, factors, srcs, indices, where);
  case 3:
    return detail::call_lincomb<3>(dst, scale, factors, srcs, indices, where);
  case 4:
    return detail::call_lincomb<4>(dst, scale, factors, srcs, indices, where);
  case 5:
    return detail::call_lincomb<5>(dst, scale, factors, srcs, indices, where);
  case 6:
    return detail::call_lincomb<6>(dst, scale, factors, srcs, indices, where);
  case 7:
    return detail::call_lincomb<7>(dst, scale, factors, srcs, indices, where);
  case 8:
    return detail::call_lincomb<8>(dst, scale, factors, srcs, indices, where);
  case 9:
    return detail::call_lincomb<9>(dst, scale, factors, srcs, indices, where);
  case 10:
    return detail::call_lincomb<10>(dst, scale, factors, srcs, indices, where);
  case 11:
    return detail::call_lincomb<11>(dst, scale, factors, srcs, indices, where);
  case 12:
    return detail::call_lincomb<12>(dst, scale, factors, srcs, indices, where);
  case 13:
    return detail::call_lincomb<13>(dst, scale, factors, srcs, indices, where);
  case 14:
    return detail::call_lincomb<14>(dst, scale, factors, srcs, indices, where);
  case 15:
    return detail::call_lincomb<15>(dst, scale, factors, srcs, indices, where);
  case 16:
    return detail::call_lincomb<16>(dst, scale, factors, srcs, indices, where);
  default:
    CCTK_VERROR("Unsupported vector length: %d", (int)NNZ);
  }
}

template void statecomp_t::lincomb<1>(const statecomp_t &dst, CCTK_REAL scale,
                                      const array<CCTK_REAL, 1> &factors,
                                      const array<const statecomp_t *, 1> &srcs,
                                      const valid_t where);
template void statecomp_t::lincomb<2>(const statecomp_t &dst, CCTK_REAL scale,
                                      const array<CCTK_REAL, 2> &factors,
                                      const array<const statecomp_t *, 2> &srcs,
                                      const valid_t where);
template void statecomp_t::lincomb<3>(const statecomp_t &dst, CCTK_REAL scale,
                                      const array<CCTK_REAL, 3> &factors,
                                      const array<const statecomp_t *, 3> &srcs,
                                      const valid_t where);
template void statecomp_t::lincomb<4>(const statecomp_t &dst, CCTK_REAL scale,
                                      const array<CCTK_REAL, 4> &factors,
                                      const array<const statecomp_t *, 4> &srcs,
                                      const valid_t where);
template void statecomp_t::lincomb<5>(const statecomp_t &dst, CCTK_REAL scale,
                                      const array<CCTK_REAL, 5> &factors,
                                      const array<const statecomp_t *, 5> &srcs,
                                      const valid_t where);
template void statecomp_t::lincomb<6>(const statecomp_t &dst, CCTK_REAL scale,
                                      const array<CCTK_REAL, 6> &factors,
                                      const array<const statecomp_t *, 6> &srcs,
                                      const valid_t where);
template void statecomp_t::lincomb<7>(const statecomp_t &dst, CCTK_REAL scale,
                                      const array<CCTK_REAL, 7> &factors,
                                      const array<const statecomp_t *, 7> &srcs,
                                      const valid_t where);
template void statecomp_t::lincomb<8>(const statecomp_t &dst, CCTK_REAL scale,
                                      const array<CCTK_REAL, 8> &factors,
                                      const array<const statecomp_t *, 8> &srcs,
                                      const valid_t where);
template void statecomp_t::lincomb<9>(const statecomp_t &dst, CCTK_REAL scale,
                                      const array<CCTK_REAL, 9> &factors,
                                      const array<const statecomp_t *, 9> &srcs,
                                      const valid_t where);
template void
statecomp_t::lincomb<10>(const statecomp_t &dst, CCTK_REAL scale,
                         const array<CCTK_REAL, 10> &factors,
                         const array<const statecomp_t *, 10> &srcs,
                         const valid_t where);
template void
statecomp_t::lincomb<11>(const statecomp_t &dst, CCTK_REAL scale,
                         const array<CCTK_REAL, 11> &factors,
                         const array<const statecomp_t *, 11> &srcs,
                         const valid_t where);
template void
statecomp_t::lincomb<12>(const statecomp_t &dst, CCTK_REAL scale,
                         const array<CCTK_REAL, 12> &factors,
                         const array<const statecomp_t *, 12> &srcs,
                         const valid_t where);
template void
statecomp_t::lincomb<13>(const statecomp_t &dst, CCTK_REAL scale,
                         const array<CCTK_REAL, 13> &factors,
                         const array<const statecomp_t *, 13> &srcs,
                         const valid_t where);
template void
statecomp_t::lincomb<14>(const statecomp_t &dst, CCTK_REAL scale,
                         const array<CCTK_REAL, 14> &factors,
                         const array<const statecomp_t *, 14> &srcs,
                         const valid_t where);
template void
statecomp_t::lincomb<15>(const statecomp_t &dst, CCTK_REAL scale,
                         const array<CCTK_REAL, 15> &factors,
                         const array<const statecomp_t *, 15> &srcs,
                         const valid_t where);
template void
statecomp_t::lincomb<16>(const statecomp_t &dst, CCTK_REAL scale,
                         const array<CCTK_REAL, 16> &factors,
                         const array<const statecomp_t *, 16> &srcs,
                         const valid_t where);

} // namespace ODESolvers
