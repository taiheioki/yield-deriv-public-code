#include <cmath>

#include "chemistry.hpp"
#include "utility.hpp"

namespace rcmc {

RateConstantMatrix rate_constant_matrix(const UndirectedNetwork &energy_network,
                                        const Real temperature) {
    using boost::multiprecision::exp;
    using boost::multiprecision::log;
    using std::exp;
    using std::log;

    const Real gamma = 1.0;

    const int n = energy_network.num_vertices();
    RateConstantMatrix K(n);

    Vector x(n);
    for (int i = 0; i < n; ++i) {
        x(i) = -energy_network.vertices[i] / (GAS * temperature);
    }
    const Real log_Z = log_sum_exp(x);

    for (int i = 0; i < n; ++i) {
        K.pi(i) =
            exp(-energy_network.vertices[i] / (GAS * temperature) - log_Z);
    }

    std::vector<Eigen::Triplet<Real>> triples;
    for (int i = 0; i < n; ++i) {
        for (const auto &[j, v] : energy_network.edges[i]) {
            const Real l_ij =
                -exp(log(gamma * BOLTZMANN * temperature / PLANCK) - log_Z -
                     v / (GAS * temperature));
            triples.emplace_back(i, j, l_ij);
        }
    }

    for (const auto &triple : triples) {
        const int i = triple.row();
        const int j = triple.col();
        K.L(i, j) = 0.0;
        K.L(j, i) = 0.0;
    }
    for (const auto &triple : triples) {
        const int i = triple.row();
        const int j = triple.col();
        const Real v = triple.value();
        K.L(j, i) += v;
    }

    for (int i = 0; i < n; ++i) {
        K.L(i, i) = -K.L().col(i).sum();
    }

    return K;
}

UndirectedNetwork rcm_to_energy_network(const RateConstantMatrix &K,
                                        const Real temperature) {
    using boost::multiprecision::log;
    using std::log;

    const int n = K.num_eq();
    UndirectedNetwork N(n);

    // The energy of EQ0 is set to the origin.
    const Real Z = Real(1.0) / K.pi(0);

    for (int i = 0; i < n; ++i) {
        N.vertices[i] = -GAS * temperature * log(Z * K.pi(i));
    }

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (K.L(i, j) != 0.0) {
                N.set_value(i, j,
                            -GAS * temperature *
                                log(-K.L(i, j) * Z * PLANCK /
                                    (BOLTZMANN * temperature)));
            }
        }
    }

    return N;
}

UndirectedNetwork derivative_by_energy(const Matrix &DL, const Vector &Dpi,
                                       const RateConstantMatrix &K,
                                       const Real temperature) {
    const int n = K.num_eq();
    UndirectedNetwork N(n);

    for (int i = 0; i < n; ++i) {
        N.vertices[i] = -Dpi(i) * K.pi(i) / (GAS * temperature);
    }

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (K.L(i, j) != 0.0) {
                N.set_value(i, j, -DL(i, j) * K.L(i, j) / (GAS * temperature));
            }
        }
    }

    return N;
}

RateConstantMatrix
compute_oss_rcm(const SparseMatrix &input_K, const Vector &eqs_pi,
                const std::vector<std::vector<int>> &groups) {
    using boost::multiprecision::abs;
    using boost::multiprecision::min;
    using std::abs;
    using std::min;

    const int num_oss = groups.size();

    RateConstantMatrix K(num_oss);

    for (int i = 0; i < num_oss; ++i) {
        K.pi(i) = Real(0.0);
        for (const int j : groups[i]) {
            K.pi(i) += eqs_pi(j);
        }
    }

    Real max_rel_err(0.0);

    for (int j = 0; j < input_K.outerSize(); ++j) {
        for (Eigen::SparseMatrix<Real>::InnerIterator itr(input_K, j); itr;
             ++itr) {
            const int i = itr.row();
            const Real Kij = itr.value();
            const Real Kji = input_K.coeff(j, i);

            if (Kij != 0.0 && Kji == 0.0) {
                K.L(i, j) = -Kij * K.pi(j);
            } else if (Kij == 0.0 && Kji != 0.0) {
                K.L(i, j) = -Kji * K.pi(i);
            } else {
                const Real rel_err = abs(Kij * K.pi(j) - Kji * K.pi(i)) /
                                     min(Kij * K.pi(j), Kji * K.pi(i));
                // if(rel_err > 0.01) {
                //     std::cerr << "Detailed balance is not satisfied."
                //               << std::endl;
                //     std::cerr << "K(" << i << ", " << j << ") = " << Kij
                //               << std::endl;
                //     std::cerr << "K(" << j << ", " << i << ") = " << Kji
                //               << std::endl;
                //     std::cerr << "pi(" << i << ") = " << K.pi(i) <<
                //     std::endl; std::cerr << "pi(" << j << ") = " << K.pi(j)
                //     << std::endl; std::cerr << "K(" << i << ", " << j << ") *
                //     pi(" << j
                //               << ") = " << Kij * K.pi(j) << std::endl;
                //     std::cerr << "K(" << j << ", " << i << ") * pi(" << i
                //               << ") = " << Kji * K.pi(i) << std::endl;
                //     std::cerr << "rel. err = " << rel_err << std::endl;

                //     std::exit(1);
                // }

                max_rel_err = std::max(max_rel_err, rel_err);
                K.L(i, j) = -(Kij * K.pi(j) + Kji * K.pi(i)) / 2.0;
            }

            const Real Lij = K.L(i, j);
            K.L(j, i) = Lij;
        }
    }

    std::cout << "Maximum relative error of detailed balance: " << max_rel_err
              << std::endl;

    for (int i = 0; i < num_oss; ++i) {
        K.L(i, i) = -K.L().col(i).sum();
    }

    return K;
}

} // namespace rcmc
