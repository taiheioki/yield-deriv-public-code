#ifndef RCMC_RCMC_HPP
#define RCMC_RCMC_HPP

#include <optional>
#include <tuple>
#include <vector>

#include "cholesky.hpp"
#include "numdef.hpp"
#include "rate_constant_matrix.hpp"

namespace rcmc {

std::tuple<std::vector<int>, std::vector<Real>, CholeskyFactor>
compute_steady_eqs(const RateConstantMatrix &K,
                   const std::optional<Real> t_max = std::nullopt);

Eigen::PermutationMatrix<Eigen::Dynamic>
get_permutation_matrix(const std::vector<int> &steady_eqs, const int n);

Vector compute_population_perm(const RateConstantMatrix &K_perm,
                               const CholeskyFactor &L_SS_factor,
                               const Vector &initial_population_perm,
                               const int num_S);

Vector compute_population(const RateConstantMatrix &K,
                          const std::vector<int> &steady_eqs,
                          const CholeskyFactor &L_SS_factor,
                          const Vector &initial_population, const int num_S);

std::vector<Vector> compute_population_history(
    const RateConstantMatrix &K, const std::vector<int> &steady_eqs,
    const std::vector<Real> &qssa_times, const CholeskyFactor &L_SS_factor,
    const Vector &initial_population,
    const std::optional<Real> t_max = std::nullopt);

std::pair<Matrix, Vector> compute_population_derivative(
    const RateConstantMatrix &K, const std::vector<int> &steady_eqs,
    const CholeskyFactor &L_SS_factor, const Vector &initial_population,
    const int product_eq);

DirectedNetwork compute_flow(const RateConstantMatrix &K,
                             const std::vector<int> &steady_eqs,
                             const CholeskyFactor &L_SS_factor,
                             const Vector &x0, const Vector &x1);
} // namespace rcmc

#endif
