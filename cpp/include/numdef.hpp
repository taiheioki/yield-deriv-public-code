#ifndef RCMC_NUMDEF_HPP
#define RCMC_NUMDEF_HPP

#include <boost/multiprecision/gmp.hpp>

// https://github.com/microsoft/vscode-cpptools/issues/7413#issuecomment-827172897
#if __INTELLISENSE__
#undef __ARM_NEON
#undef __ARM_NEON__
#endif

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace rcmc {

// using Real = double;
using Real =
    boost::multiprecision::number<boost::multiprecision::gmp_float<308>>;

using Matrix = Eigen::Matrix<rcmc::Real, Eigen::Dynamic, Eigen::Dynamic>;
using Vector = Eigen::Matrix<Real, Eigen::Dynamic, 1>;
using RowVector = Eigen::Matrix<Real, 1, Eigen::Dynamic>;

using SparseMatrix = Eigen::SparseMatrix<Real>;
using SparseVector = Eigen::SparseVector<Real>;
} // namespace rcmc

#endif
