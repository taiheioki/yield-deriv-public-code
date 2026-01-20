#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <algorithm>
#include <cassert>

#include "numdef.hpp"

namespace rcmc {

// Set each diagonals to the minus of the sum of off-diagonals in each column.
// Original diagonals are ignored.
inline void set_diagonals(Matrix &A) {
    assert(A.rows() == A.cols());
    const int n = A.rows();

    for (int i = 0; i < n; ++i) {
        A(i, i) = 0.0;
        A(i, i) = -A.col(i).sum();
    }
}

// Count the number of nonzero entries
inline int nnz(const Matrix &A) {
    return A.size() - std::count(A.data(), A.data() + A.size(), 0.0);
}

// Count the number of nonzero entries
inline int nnz(const SparseMatrix &A) { return A.nonZeros(); }

// e_i \in R^n
inline Vector standard_vector(const int n, const int i) {
    assert(0 <= i && i < n);
    Vector e_i = Vector::Zero(n);
    e_i(i) = 1.0;
    return e_i;
}

inline Real log_sum_exp(const Vector &x) {
    using boost::multiprecision::log;
    using std::log;

    const Real max_x = x.maxCoeff();
    return max_x + log((x.array() - max_x).exp().sum());
}

} // namespace rcmc

#endif
