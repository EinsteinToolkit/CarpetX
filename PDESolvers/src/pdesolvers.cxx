#include <iostream> //TODO

#include "pdesolvers.hxx"

#if defined _OPENMP
#include <omp.h>
#elif defined __HIPCC__
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
#else
static inline int omp_get_max_threads() { return 1; }
static inline int omp_get_thread_num() { return 0; }
#endif

#include <cassert>
#include <cstring>

namespace PDESolvers {

bool csr_t::invariant() const {
  assert(m >= 0);
  assert(n >= 0);
  assert(int(rowptrs.size()) == m + 1);
  assert(int(colvals.size()) == rowptrs.back());
  assert(int(nzvals.size()) == rowptrs.back());
  assert(rowptrs.at(0) == 0);
  for (int i = 0; i < m; ++i) {
    assert(rowptrs.at(i) <= rowptrs.at(i + 1));
    for (int jp = rowptrs.at(i); jp < rowptrs.at(i + 1); ++jp) {
      const int j = colvals.at(jp);
      assert(j >= 0 && j < n);
      // We could check that j are sorted
    }
  }
  return true;
}

void csr_t::insert_element(const int i, const int j, const CCTK_REAL v) {
  if (!(i >= 0 && i < m))
    std::cout << "m=" << m << " n=" << n << " i=" << i << " j=" << j
              << " v=" << v << "\n";
  assert(i >= 0 && i < m);
  assert(j >= 0 && j < n);
  assert(i >= int(rowptrs.size()) - 1);
  while (int(rowptrs.size()) <= i)
    rowptrs.push_back(colvals.size());
  colvals.push_back(j);
  nzvals.push_back(v);
}
void csr_t::finish_inserting() {
  assert(int(rowptrs.size()) <= m);
  while (int(rowptrs.size()) <= m)
    rowptrs.push_back(colvals.size());
}

csr_t::csr_t(
    const int m, const int n,
    const std::vector<
        std::vector<std::vector<std::tuple<int, int, CCTK_REAL> > > > &values)
    : m(m), n(n) {
  for (const auto &values1 : values) {
    for (const auto &values2 : values1) {
      for (const auto &ijv : values2) {
        const int i = std::get<0>(ijv);
        const int j = std::get<1>(ijv);
        const CCTK_REAL v = std::get<2>(ijv);
        insert_element(i, j, v);
      }
    }
  }
  finish_inserting();
  assert(invariant());
}

csr_t::csr_t(const int m, const int n,
             const Arith::spvect<std::tuple<int, int>, CCTK_REAL> &values)
    : m(m), n(n) {
  for (const auto &ijv : values) {
    const int i = std::get<0>(ijv.first);
    const int j = std::get<1>(ijv.first);
    const CCTK_REAL v = ijv.second;
    insert_element(i, j, v);
  }
  finish_inserting();
  assert(invariant());
}

std::size_t csr_t::size() const { return nzvals.size(); }

// void csr_t::count_nz(int ilocal_min, int ilocal_max, int &restrict nlocal,
//                      int &restrict ntotal) const {
//   for (int i = 0; i < m; ++i) {
//     const int ncols = rowptrs.at(i + 1) - rowptrs.at(i);
//     if (i >= ilocal_min && i < ilocal_max)
//       nlocal += ncols;
//     ntotal += ncols;
//   }
// }

////////////////////////////////////////////////////////////////////////////////

std::size_t jacobian_t::size() const { return entries.size(); }

// void jacobian_t::count_nz(int ilocal_min, int ilocal_max, int &restrict
// nlocal,
//                           int &restrict ntotal) const {
//   for (const auto &e : entries) {
//     const auto i = std::get<0>(e);
//     nlocal += i >= ilocal_min && i < ilocal_max;
//     ++ntotal;
//   }
// }

void jacobian_t::clear() { entries.clear(); }

void jacobian_t::count_matrix_entries(const csr_t &Jp, int ilocal_min,
                                      int ilocal_max, int &restrict nlocal,
                                      int &restrict ntotal) const {
  const auto count_point = [&](const int i, const int j) {
    nlocal += j >= ilocal_min && j < ilocal_max;
    ntotal += j >= 0;
  };
  for (const auto &e : entries) {
    const auto i = std::get<0>(e);
    const auto j = std::get<1>(e);
    assert(i >= 0);
    if (j < 0) {
      // ignore this point
    } else if (j >= 0 && j < prolongation_index_offset) {
      // regular point
      count_point(i, j);
    } else {
      // prolongated point
      const int row = j - prolongation_index_offset;
      assert(0 <= row && row < Jp.m);
      const int rowptr0 = Jp.rowptrs.at(row);
      const int rowptr1 = Jp.rowptrs.at(row + 1);
      for (int colindp = 0; colindp < rowptr1 - rowptr0; ++colindp)
        count_point(i, Jp.colvals.at(rowptr0 + colindp));
    }
  }
}

void jacobian_t::set_matrix_entries(const csr_t &Jp, Mat J) const {
  std::vector<CCTK_REAL> values;
  for (const auto &e : entries) {
    const auto i = std::get<0>(e);
    const auto j = std::get<1>(e);
    const auto v = std::get<2>(e);
    assert(i >= 0);
    if (j < 0) {
      // ignore this point
    } else if (j >= 0 && j < prolongation_index_offset) {
      // regular point
      MatSetValue(J, i, j, v, ADD_VALUES);
    } else {
      // prolongated point
      const int row = j - prolongation_index_offset;
      assert(0 <= row && row < Jp.m);
      const int rowptr0 = Jp.rowptrs.at(row);
      const int rowptr1 = Jp.rowptrs.at(row + 1);
      values.resize(rowptr1 - rowptr0);
      for (int rowptr = rowptr0; rowptr < rowptr1; ++rowptr)
        values.at(rowptr - rowptr0) = v * Jp.nzvals.at(rowptr);
      MatSetValues(J, 1, &i, rowptr1 - rowptr0, &Jp.colvals.at(rowptr0),
                   values.data(), ADD_VALUES);
    }
  }
}

void jacobian_t::set_matrix_entries(
    const csr_t &Jp,
    Arith::spvect<std::tuple<int, int>, CCTK_REAL> &Jsp) const {
  for (const auto &e : entries) {
    const auto i = std::get<0>(e);
    const auto j = std::get<1>(e);
    const auto v = std::get<2>(e);
    assert(i >= 0);
    if (j < 0) {
      // ignore this point
    } else if (j >= 0 && j < prolongation_index_offset) {
      // regular point
      Jsp.emplace_back(std::make_tuple(i, j), v);
    } else {
      // prolongated point
      const int row = j - prolongation_index_offset;
      assert(0 <= row && row < Jp.m);
      const int rowptr0 = Jp.rowptrs.at(row);
      const int rowptr1 = Jp.rowptrs.at(row + 1);
      for (int rowptr = rowptr0; rowptr < rowptr1; ++rowptr) {
        const int col = Jp.colvals.at(rowptr);
        const CCTK_REAL val = Jp.nzvals.at(rowptr);
        Jsp.emplace_back(std::make_tuple(i, col), v * val);
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

jacobians_t::jacobians_t() : jacobians(omp_get_max_threads()) {}

jacobian_t &jacobians_t::get_local() {
  return jacobians.at(omp_get_thread_num());
}

std::size_t jacobians_t::size() const {
  std::size_t sz = 0;
  for (auto &j : jacobians)
    sz += j.size();
  return sz;
}

// void jacobians_t::count_nz(int ilocal_min, int ilocal_max, int &restrict
// nlocal,
//                            int &restrict ntotal) const {
//   for (auto &j : jacobians)
//     j.count_nz(ilocal_min, ilocal_max, nlocal, ntotal);
// }

void jacobians_t::clear() {
  for (auto &j : jacobians)
    j.clear();
}

void jacobians_t::define_matrix(const csr_t &Jp, Mat J) const {
  PetscErrorCode ierr;
  ierr = MatZeroEntries(J);
  assert(!ierr);
  for (const auto &j : jacobians)
    j.set_matrix_entries(Jp, J);
  ierr = MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
  assert(!ierr);
  ierr = MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
  assert(!ierr);
  ierr = MatSetOption(J, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
  assert(!ierr);

#if 0
  PetscErrorCode ierr;
  int ilocal_min, ilocal_max;
  {
    MatInfo info;
    PetscErrorCode ierr;
    ierr = MatGetInfo(J, MAT_GLOBAL_SUM, &info);
    assert(!ierr);
    CCTK_VINFO("Jacobian info: nz_allocated=%g nz_used=%g nz_unneeded=%g "
               "memory=%g",
               double(info.nz_allocated), double(info.nz_used),
               double(info.nz_unneeded), double(info.memory));
  }
  // ierr = MatSetUp(J);
  // assert(!ierr);
  ierr = MatGetOwnershipRange(J, &ilocal_min, &ilocal_max);
  assert(!ierr);
  int nlocal = 0, ntotal = 0;
  for (const auto &j : jacobians)
    j.count_matrix_entries(Jp, ilocal_min, ilocal_max, nlocal, ntotal);
  const int dnz = nlocal;
  const int onz = ntotal - nlocal;
  ierr = MatMPIAIJSetPreallocation(J, dnz, NULL, onz, NULL);
  assert(!ierr);
  {
    MatInfo info;
    PetscErrorCode ierr;
    ierr = MatGetInfo(J, MAT_GLOBAL_SUM, &info);
    assert(!ierr);
    CCTK_VINFO("Jacobian info: nz_allocated=%g nz_used=%g nz_unneeded=%g "
               "memory=%g",
               double(info.nz_allocated), double(info.nz_used),
               double(info.nz_unneeded), double(info.memory));
  }
  // ierr = MatSetOption(J, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
  // assert(!ierr);
  ierr = MatZeroEntries(J);
  assert(!ierr);
  for (const auto &j : jacobians)
    j.set_matrix_entries(Jp, J);
  ierr = MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
  assert(!ierr);
  ierr = MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
  assert(!ierr);
#endif

#if 0
  PetscErrorCode ierr;
  Arith::spvect<std::tuple<int, int>, CCTK_REAL> Jsp;
  for (const auto &j : jacobians)
    j.set_matrix_entries(Jp, Jsp);
  Jsp.make_sorted();
  csr_t Jcsr(Jp.n, Jp.n, Jsp);
  Jsp = Arith::spvect<std::tuple<int, int>, CCTK_REAL>(); // release storage

  const char *type;
  ierr = MatGetType(J, &type);
  assert(!ierr);
  if (std::strcmp(type, MATSEQAIJ) == 0) {
    ierr = MatSeqAIJSetPreallocationCSR(
        J, Jcsr.rowptrs.data(), Jcsr.colvals.data(), Jcsr.nzvals.data());
    assert(!ierr);
  } else if (std::strcmp(type, MATMPIAIJ) == 0) {
    ierr = MatMPIAIJSetPreallocationCSR(
        J, Jcsr.rowptrs.data(), Jcsr.colvals.data(), Jcsr.nzvals.data());
    assert(!ierr);
  } else {
    assert(0);
  }
  // We don't need to zero since we insert instead of adding
  // ierr = MatZeroEntries(J);
  // assert(!ierr);
  for (int i = 0; i < Jcsr.m; ++i) {
    const int jp0 = Jcsr.rowptrs.at(i);
    const int jp1 = Jcsr.rowptrs.at(i + 1);
    for (int jp = jp0; jp < jp1; ++jp)
      assert(Jcsr.colvals.at(jp) >= 0);
    ierr = MatSetValues(J, 1, &i, jp1 - jp0, &Jcsr.colvals.at(jp0),
                        &Jcsr.nzvals.at(jp0), INSERT_VALUES);
    assert(!ierr);
  }
  Jcsr = csr_t(); // release storage

  ierr = MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
  assert(!ierr);
  ierr = MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
  assert(!ierr);
  ierr = MatSetOption(J, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
  assert(!ierr);
#endif
}

////////////////////////////////////////////////////////////////////////////////

std::optional<jacobians_t> jacobians;

} // namespace PDESolvers
