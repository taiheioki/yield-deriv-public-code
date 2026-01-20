#ifndef RCMC_RATE_CONSTANT_MATRIX_HPP
#define RCMC_RATE_CONSTANT_MATRIX_HPP

#include "network.hpp"
#include "numdef.hpp"
#include "utility.hpp"

namespace rcmc {
class RateConstantMatrix {
  private:
    SparseMatrix m_L;
    Vector m_pi;

  public:
    explicit RateConstantMatrix(const int n)
        : m_L(n, n), m_pi(Vector::Zero(n)) {}

    RateConstantMatrix(const SparseMatrix &L, const Vector &pi)
        : m_L(L), m_pi(pi) {}

    SparseMatrix &L() { return m_L; }

    const SparseMatrix &L() const { return m_L; }

    Real &L(const int i, const int j) { return L().coeffRef(i, j); }

    Real L(const int i, const int j) const { return L().coeff(i, j); }

    Vector &pi() { return m_pi; }

    const Vector &pi() const { return m_pi; }

    Real &pi(const int i) { return pi()(i); }

    Real pi(const int i) const { return pi()(i); }

    auto Pi() const { return pi().asDiagonal(); }

    auto K() const { return -L() * Pi().inverse(); }

    Real K(const int i, const int j) const { return -L(i, j) / pi(j); }

    int num_eq() const { return pi().size(); }

    int num_ts() const { return (nnz(L()) - num_eq()) / 2; }

    RateConstantMatrix
    permuted(const Eigen::PermutationMatrix<Eigen::Dynamic> &P) const {
        return RateConstantMatrix(P * L() * P.inverse(), P * pi());
    }

    UndirectedNetwork to_undirected_network() const;

    DirectedNetwork to_network() const;
};
} // namespace rcmc

#endif
