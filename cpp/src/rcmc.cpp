#include <cassert>
#include <iostream>
#include <unordered_set>
#include <utility>

#include "rcmc.hpp"
#include "utility.hpp"

namespace rcmc {

std::tuple<std::vector<int>, std::vector<Real>, CholeskyFactor>
compute_steady_eqs(const RateConstantMatrix &K,
                   const std::optional<Real> t_max) {
    const int n = K.num_eq();

    auto P = Eigen::PermutationMatrix<Eigen::Dynamic>(n);
    P.setIdentity();
    auto P_inv = P;

    // Lower uni-triangular. Cholesky factorization of L
    Matrix B = Matrix::Zero(n, n);

    Matrix L = K.L();   // Working matrix of L
    Vector pi = K.pi(); // Working vector of pi

    std::vector<int> steady_eqs;
    std::vector<Real> qssa_times;

    for (int i = 0; i < n - 1; ++i) {
        // Find the largest entry
        Real largest = 0.0;
        int pivot = -1;
        for (int j = i; j < n; ++j) {
            const Real value = L(j, j) / pi(j);
            if (largest < value) {
                largest = value;
                pivot = j;
            }
        }

        if (pivot == -1) {
            break;
        }
        const Real t = 1.0 / largest;

        if (t_max && t_max <= t) {
            break;
        }

        // Swap the i-th and pivot-th rows and columns
        B.row(i).head(i).swap(B.row(pivot).head(i));
        B.col(i).head(i).swap(B.col(pivot).head(i));
        L.row(i).tail(n - i).swap(L.row(pivot).tail(n - i));
        L.col(i).tail(n - i).swap(L.col(pivot).tail(n - i));
        std::swap(pi(i), pi(pivot));
        P.applyTranspositionOnTheLeft(i, pivot);
        P_inv.applyTranspositionOnTheRight(i, pivot);

        B.col(i).tail(n - i - 1) = L.col(i).tail(n - i - 1) / L(i, i);
        B(i, i) = 1.0;

        // Eliminate each column
        L.col(i).tail(n - i - 1).setZero();
        for (int j = i + 1; j < n; ++j) {
            if (L(i, j) != 0.0) {
                L.col(j).tail(n - i - 1) -= L(i, j) * B.col(i).tail(n - i - 1);
                L(i, j) = 0.0;
                L(j, j) = 0.0;
                L(j, j) = -L.col(j).tail(n - i - 1).sum();
            }
        }

        // Insert the epoch
        steady_eqs.push_back(P_inv.indices()[i]);
        qssa_times.push_back(t);
    }

    // By the following, B_SS * D_S * B_SS.inverse() is equal to L_SS
    const int num_S = steady_eqs.size();
    const Matrix B_S = B.leftCols(num_S);
    const Vector d_S = L.diagonal().head(num_S);

    return {steady_eqs, qssa_times, CholeskyFactor(B_S, d_S)};
}

Vector compute_population_perm(const RateConstantMatrix &K_perm,
                               const CholeskyFactor &L_SS_factor,
                               const Vector &initial_population_perm,
                               const int num_S) {
    const int n = K_perm.num_eq();
    const int num_T = n - num_S;

    const auto pi_S = K_perm.pi().head(num_S);
    const auto pi_T = K_perm.pi().tail(num_T);
    const auto Pi_S = pi_S.asDiagonal();
    const auto Pi_T = pi_T.asDiagonal();

    const auto L_ST = K_perm.L().topRightCorner(num_S, num_T);
    const auto L_TS = K_perm.L().bottomLeftCorner(num_T, num_S);

    const auto y_S = initial_population_perm.head(num_S);
    const auto y_T = initial_population_perm.tail(num_T);

    Vector z(n);
    auto z_S = z.head(num_S);
    auto z_T = z.tail(num_T);

    // cppcheck-suppress redundantInitialization
    z_T = (y_T - L_TS * L_SS_factor.solve(y_S))
              .cwiseQuotient(pi_T - L_TS * L_SS_factor.solve(pi_S));

    // cppcheck-suppress [unreadVariable, redundantInitialization]
    z_S = -1.0 * Pi_S * L_SS_factor.solve(L_ST * z_T);
    // cppcheck-suppress unreadVariable
    z_T = Pi_T * z_T;

    return z;
}

Vector compute_population_perm2(const RateConstantMatrix &K_perm,
                                const CholeskyFactor &L_SS_factor,
                                const Vector &initial_population_perm,
                                const int num_S) {
    const int n = K_perm.num_eq();
    const int num_T = n - num_S;

    const auto pi_S = K_perm.pi().head(num_S);
    const auto pi_T = K_perm.pi().tail(num_T);
    const auto Pi_S = pi_S.asDiagonal();
    const auto Pi_T = pi_T.asDiagonal();

    const auto L_ST = K_perm.L().topRightCorner(num_S, num_T);
    const auto L_TS = K_perm.L().bottomLeftCorner(num_T, num_S);

    const auto y_S = initial_population_perm.head(num_S);
    const auto y_T = initial_population_perm.tail(num_T);

    Vector z(n);
    auto z_S = z.head(num_S);
    auto z_T = z.tail(num_T);

    const Vector psi_y =
        y_S + Pi_S * L_SS_factor.solve(L_ST * Pi_T.inverse() * y_T);
    z_S = psi_y -
          Pi_S * L_SS_factor.solve(
                     L_ST *
                     (L_TS * L_SS_factor.solve(psi_y))
                         .cwiseQuotient(pi_T - L_TS * L_SS_factor.solve(pi_S)));
    z_T = L_TS * L_SS_factor.solve(z_S);

    return initial_population_perm - z;
}

// Construct a permutation matrix so that the first indices are
// steady EQs
Eigen::PermutationMatrix<Eigen::Dynamic>
get_permutation_matrix(const std::vector<int> &steady_eqs, const int n) {
    auto P = Eigen::PermutationMatrix<Eigen::Dynamic>(n);
    std::unordered_set<int> S;

    int cnt = 0;
    for (const int eq : steady_eqs) {
        P.indices()[eq] = cnt;
        S.insert(eq);
        ++cnt;
    }

    for (int i = 0; i < n; ++i) {
        if (S.find(i) == S.end()) {
            P.indices()[i] = cnt;
            ++cnt;
        }
    }

    return P;
}

Vector compute_population(const RateConstantMatrix &K,
                          const std::vector<int> &steady_eqs,
                          const CholeskyFactor &L_SS_factor,
                          const Vector &initial_population, const int num_S) {
    const auto P = get_permutation_matrix(steady_eqs, K.num_eq());
    const auto K_perm = K.permuted(P);
    const Vector initial_pop_perm = P * initial_population;

    return P.inverse() * compute_population_perm(K_perm, L_SS_factor,
                                                 initial_pop_perm, num_S);
}

std::vector<Vector> compute_population_history(
    const RateConstantMatrix &K, const std::vector<int> &steady_eqs,
    const std::vector<Real> &qssa_times, const CholeskyFactor &L_SS_factor,
    const Vector &initial_population, const std::optional<Real> t_max) {
    const auto P = get_permutation_matrix(steady_eqs, K.num_eq());
    const auto K_perm = K.permuted(P);
    const Vector initial_pop_perm = P * initial_population;

    std::vector<Vector> history;
    for (int num_S = 1; num_S <= int(steady_eqs.size()) &&
                        (!t_max || qssa_times[num_S - 1] <= t_max);
         ++num_S) {
        const auto z = compute_population_perm2(K_perm, L_SS_factor,
                                                initial_pop_perm, num_S);
        history.emplace_back(P.inverse() * z);
    }

    return history;
}

std::pair<Matrix, Vector> compute_population_derivative(
    const RateConstantMatrix &K, const std::vector<int> &steady_eqs,
    const CholeskyFactor &L_SS_factor, const Vector &initial_population,
    const int product_eq) {
    // Permute
    const int n = K.num_eq();
    const auto P = get_permutation_matrix(steady_eqs, n);
    const auto K_perm = K.permuted(P);
    const Vector initial_pop_perm = P * initial_population;
    const int product_eq_perm = P.indices()[product_eq];

    // Separate to parts
    const int num_S = steady_eqs.size();
    const int num_T = n - num_S;

    const auto K_matrix = K_perm.K();
    const auto K_ST = K_matrix.topRightCorner(num_S, num_T);
    const auto K_TS = K_matrix.bottomLeftCorner(num_T, num_S);

    const auto &pi = K_perm.pi();
    const auto pi_S = pi.head(num_S);
    const auto pi_T = pi.tail(num_T);
    const auto Pi_S = pi_S.asDiagonal();
    const auto Pi_S_inv = Pi_S.inverse();

    const auto y_S = initial_pop_perm.head(num_S);
    const auto y_T = initial_pop_perm.tail(num_T);

    Matrix DK_tr = Matrix::Zero(n, n);
    auto DK_SS_tr = DK_tr.topLeftCorner(num_S, num_S);
    auto DK_ST_tr = DK_tr.bottomLeftCorner(num_T, num_S);
    auto DK_TS_tr = DK_tr.topRightCorner(num_S, num_T);

    Vector Dpi = Vector::Zero(n);
    auto Dpi_S = Dpi.head(num_S);
    auto Dpi_T = Dpi.tail(num_T);

    // K_SS.inverse() * b
    const auto lsolve = [&](const auto &b) {
        return (Pi_S * L_SS_factor.solve(-b)).eval();
    };

    // b * K_SS.inverse()
    const auto rsolve = [&](const auto &b) {
        return L_SS_factor.solve(-1.0 * Pi_S * b.transpose())
            .transpose()
            .eval();
    };

    const RowVector sol_1_tr = -L_SS_factor.solve(pi_S).transpose();
    const auto sol_y_S = lsolve(y_S);

    if (product_eq_perm < num_S) { // target is in S
        const int p = product_eq_perm;

        Vector v_TT(num_T);
        for (int i = 0; i < num_T; ++i) {
            v_TT(i) =
                1.0 / (rcmc::Real(1.0) - (sol_1_tr * K_ST.col(i)).coeff(0, 0));
        }
        const auto V_TT = v_TT.asDiagonal();

        const Vector z_T = V_TT * (y_T - K_TS * sol_y_S);
        const auto sol_e_p_tr = rsolve(standard_vector(num_S, p).transpose());
        const RowVector a = sol_e_p_tr * K_ST * V_TT;
        const Vector b = a.asDiagonal() * z_T;
        const auto c = lsolve(K_ST * z_T);
        const auto d = lsolve(K_ST * b);

        // cppcheck-suppress redundantInitialization
        DK_SS_tr = c * sol_e_p_tr + d * sol_1_tr - sol_y_S * rsolve(a * K_TS);
        DK_ST_tr = -z_T * sol_e_p_tr - b * sol_1_tr;
        DK_TS_tr = sol_y_S * a;

        // cppcheck-suppress redundantInitialization
        Dpi_S = -1.0 * Pi_S_inv * d;
        Dpi_S(p) -= c(p) / pi_S(p);
        for (int j = 0; j < num_T; ++j) {
            Dpi_T(j) = -(DK_ST_tr.row(j) * K_ST.col(j)).coeff(0, 0) / pi_T(j);
        }
    } else { // target is in T
        const int p = product_eq_perm - num_S;
        const auto k_Sp = K_ST.col(p);
        const auto k_pS = K_TS.row(p);
        const Real v_pp =
            1.0 / (rcmc::Real(1.0) - (sol_1_tr * k_Sp).coeff(0, 0));
        const Real z_p = v_pp * (y_T(p) - (k_pS * sol_y_S).coeff(0, 0));
        const auto sol_K_Sp = lsolve(k_Sp);

        DK_ST_tr.row(p) = v_pp * z_p * sol_1_tr;
        DK_TS_tr.col(p) = -v_pp * sol_y_S;
        DK_SS_tr = -sol_K_Sp * DK_ST_tr.row(p) - DK_TS_tr.col(p) * rsolve(k_pS);

        Dpi_S = v_pp * z_p * Pi_S_inv * sol_K_Sp;
        Dpi_T(p) = -(DK_ST_tr.row(p) * k_Sp).coeff(0, 0) / pi_T(p);
    }

    Matrix DL(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            const Real DL_ij = (DK_tr(i, i) - DK_tr(i, j)) / pi(i) +
                               (DK_tr(j, j) - DK_tr(j, i)) / pi(j);
            DL(i, j) = DL_ij;
            DL(j, i) = DL_ij;
        }
    }

    return {P.inverse() * DL * P, P.inverse() * Dpi};
}

DirectedNetwork compute_flow(const RateConstantMatrix &K,
                             const std::vector<int> &steady_eqs,
                             const CholeskyFactor &L_SS_factor,
                             const Vector &x0, const Vector &x1) {
    const auto n = K.num_eq();
    const auto P = get_permutation_matrix(steady_eqs, n);
    const auto K_perm = K.permuted(P);
    const Vector x_perm = P * (x0 - x1);

    Vector y_perm(n);

    const Vector pi_all = K_perm.pi();
    const Vector pi_head = pi_all.head(n - 1).eval();
    Vector tmp = L_SS_factor.solve(x_perm.head(n - 1)).eval();
    y_perm.head(n - 1) = (pi_head.array() * tmp.array()).matrix();
    y_perm(n - 1) = 0.0;

    const Vector y = P.inverse() * y_perm;

    DirectedNetwork flow(n);
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (K.L(i, j) != 0.0) {
                flow.set_value(i, j, K.K(j, i) * y(i) - K.K(i, j) * y(j));
            }
        }
    }

    return flow;
}

} // namespace rcmc
