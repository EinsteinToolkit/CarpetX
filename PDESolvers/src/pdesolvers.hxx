#ifndef CARPETX_PDESOLVERS_PDESOLVERS_HXX
#define CARPETX_PDESOLVERS_PDESOLVERS_HXX

#include <spvect.hxx>

#include <cctk.h>

#include <petsc.h>

#include <cassert>
#include <optional>
#include <tuple>
#include <vector>

namespace PDESolvers {

struct csr_t {
  // Adapted from Julia SparseMatrixCSC
  int m;                         // rows
  int n;                         // columns
  std::vector<int> rowptrs;      // row i is in rowptrs[i] ... rowptrs[i+1]
  std::vector<int> colvals;      // column indices of stored values
  std::vector<CCTK_REAL> nzvals; // stored values

  csr_t(const csr_t &) = default;
  csr_t(csr_t &&) = default;
  csr_t &operator=(const csr_t &) = default;
  csr_t &operator=(csr_t &&) = default;

  csr_t() : m(0), n(0) {}

private:
  void insert_element(int i, int j, CCTK_REAL v);
  void finish_inserting();

public:
  csr_t(int m, int n,
        const std::vector<
            std::vector<std::vector<std::tuple<int, int, CCTK_REAL> > > >
            &values);
  csr_t(int m, int n,
        const Arith::spvect<std::tuple<int, int>, CCTK_REAL> &values);

  bool invariant() const;

  // void count_nz(int ilocal_min, int ilocal_max, int &restrict nlocal,
  //               int &restrict ntotal) const;
  std::size_t size() const;
};

////////////////////////////////////////////////////////////////////////////////

const int prolongation_index_offset = 0x40000000U;

class jacobian_t {
  std::vector<std::tuple<int, int, CCTK_REAL> > entries;

public:
  jacobian_t(const jacobian_t &) = default;
  jacobian_t(jacobian_t &&) = default;
  jacobian_t &operator=(const jacobian_t &) = default;
  jacobian_t &operator=(jacobian_t &&) = default;

  jacobian_t() = default;
  std::size_t size() const;
  // void count_nz(int ilocal_min, int ilocal_max, int &restrict nlocal,
  //               int &restrict ntotal) const;
  void clear();
  void add_value(int i, int j, CCTK_REAL v) { entries.emplace_back(i, j, v); }
  void count_matrix_entries(const csr_t &Jp, int ilocal_min, int ilocal_max,
                            int &restrict nlocal, int &restrict ntotal) const;
  void set_matrix_entries(const csr_t &Jp, Mat J) const;
  void
  set_matrix_entries(const csr_t &Jp,
                     Arith::spvect<std::tuple<int, int>, CCTK_REAL> &J) const;
};

////////////////////////////////////////////////////////////////////////////////

class jacobians_t {
  std::vector<jacobian_t> jacobians;

public:
  jacobians_t(const jacobians_t &) = delete;
  jacobians_t(jacobians_t &&) = default;
  jacobians_t &operator=(const jacobians_t &) = delete;
  jacobians_t &operator=(jacobians_t &&) = default;

  jacobians_t();
  jacobian_t &get_local();
  std::size_t size() const;
  // void count_nz(int ilocal_min, int ilocal_max, int &restrict nlocal,
  //               int &restrict ntotal) const;
  void clear();
  void define_matrix(const csr_t &Jp, Mat J) const;
};

////////////////////////////////////////////////////////////////////////////////

extern std::optional<jacobians_t> jacobians;

} // namespace PDESolvers

#endif // #ifndef CARPETX_PDESOLVERS_PDESOLVERS_HXX
