#ifndef RCMC_CHEMISTRY_HPP
#define RCMC_CHEMISTRY_HPP

#include "network.hpp"
#include "rate_constant_matrix.hpp"

namespace rcmc {

const static Real AVOGADRO = 6.02214076e23;   // mol^{-1}
const static Real BOLTZMANN = 1.380649e-26;   // kJ K^{-1}
const static Real GAS = AVOGADRO * BOLTZMANN; // kJ K^{-1} mol^{-1}
const static Real PLANCK = 6.62607015e-37;    // kJ s

RateConstantMatrix rate_constant_matrix(const UndirectedNetwork &energy_network,
                                        const Real temperature);

UndirectedNetwork rcm_to_energy_network(const RateConstantMatrix &K,
                                        const Real temperature);

UndirectedNetwork derivative_by_energy(const Matrix &DL, const Vector &Dpi,
                                       const RateConstantMatrix &K,
                                       const Real temperature);

RateConstantMatrix compute_oss_rcm(const SparseMatrix &input_K,
                                   const Vector &eqs_pi,
                                   const std::vector<std::vector<int>> &groups);

} // namespace rcmc

#endif
