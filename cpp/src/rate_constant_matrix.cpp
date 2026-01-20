#include "rate_constant_matrix.hpp"

namespace rcmc {

UndirectedNetwork RateConstantMatrix::to_undirected_network() const {
    const int n = num_eq();
    UndirectedNetwork N(n);

    for (int i = 0; i < n; ++i) {
        N.vertices[i] = pi(i);
    }

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (L(j, i) != 0.0) {
                N.set_value(i, j, -L(j, i));
            }
        }
    }

    return N;
}

DirectedNetwork RateConstantMatrix::to_network() const {
    const int n = num_eq();
    DirectedNetwork N(n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j && L(j, i) != 0.0) {
                N.set_value(i, j, -L(j, i) / pi(i));
            }
        }
    }

    return N;
}
} // namespace rcmc
