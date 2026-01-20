#ifndef CHOLESKY_HPP
#define CHOLESKY_HPP

#include <cassert>

#include "numdef.hpp"

namespace rcmc {

// Decomposition of a symmetric matrix A to B * d.diagonal() * B.transpose()
class CholeskyFactor {
  private:
    SparseMatrix m_B; // Lower uni-triangular
    Vector m_d;       // Vector of diagonals
    int num_S;

  public:
    CholeskyFactor(const Matrix &B, const Vector &d)
        : m_B(B.sparseView()), m_d(d), num_S(d.size()) {
        assert(B.rows() >= num_S && B.cols() == num_S);
    }

    template <class Derived>
    CholeskyFactor(const Eigen::MatrixBase<Derived> &L, const int num_S)
        : m_B(L), m_d(num_S), num_S(num_S) {
        const int n = L.rows();
        assert(L.cols() == n);
        assert(num_S <= n);

        Matrix M = L;

        for (int i = 0; i < num_S; ++i) {
            m_d(i) = M(i, i);
            M.col(i).tail(n - i) /= m_d(i);
            for (int j = i + 1; j < n; ++j) {
                if (M(i, j) != 0.0) {
                    M.col(j).tail(n - i - 1) -=
                        M(i, j) * M.col(i).tail(n - i - 1);
                    M(i, j) = 0.0;
                    M(j, j) = 0.0;
                    M(j, j) = -M.col(j).tail(n - i - 1).sum();
                }
            }
        }

        m_B = M.sparseView();
    }

    int size() const { return num_S; }

    const SparseMatrix &B() const { return m_B; }

    const Vector &d() const { return m_d; }

    auto D() const { return d().asDiagonal(); }

    Matrix product() const { return B() * D() * B().transpose(); }

    // Solve A.topLeftCorner(b.size(), b.size()) * x = b
    Vector solve(const Vector &b) const {
        assert(b.size() <= size());
        const auto B_SS = B().topLeftCorner(b.size(), b.size());
        const auto B_SS_sol = B_SS.triangularView<Eigen::UnitLower>();
        const auto D_S_inv = d().head(b.size()).asDiagonal().inverse();
        const auto B_SS_t_sol =
            B_SS.transpose().triangularView<Eigen::UnitUpper>();
        return B_SS_t_sol.solve(D_S_inv * B_SS_sol.solve(b));
    }
};

} // namespace rcmc

#endif
