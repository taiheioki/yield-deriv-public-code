#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

#include "rcmc.hpp"

namespace {

void require(const bool condition, const std::string &message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void check_population_history(const int center, const int initial) {
    // The first contraction of this three-state star used to produce -1/6
    // at the initially unoccupied leaf in the population history.
    rcmc::Matrix L = rcmc::Matrix::Zero(3, 3);
    for (int leaf = 0; leaf < 3; ++leaf) {
        if (leaf != center) {
            L(center, leaf) = L(leaf, center) = -1;
            L(leaf, leaf) = 1;
        }
    }
    L(center, center) = 2;
    const rcmc::RateConstantMatrix K(L.sparseView(), rcmc::Vector::Ones(3));
    const auto [steady_eqs, times, factor] = rcmc::compute_steady_eqs(K);
    const auto y = rcmc::standard_vector(3, initial);
    const auto history =
        rcmc::compute_population_history(K, steady_eqs, times, factor, y);
    const auto context = "center EQ" + std::to_string(center) + ": ";
    const rcmc::Real tolerance = 1e-12;

    require(history.size() == 2, context + "expected two contractions");
    require(steady_eqs.front() == center,
            context + "expected the center to be contracted first");

    rcmc::Vector expected = rcmc::Vector::Zero(3);
    expected(center) = rcmc::Real(1) / 3;
    expected(initial) = rcmc::Real(2) / 3;
    require((history.front() - expected).cwiseAbs().maxCoeff() < tolerance,
            context + "first population differs from the known Type A result");

    for (std::size_t step = 0; step < history.size(); ++step) {
        const auto &population = history[step];
        const auto single = rcmc::compute_population(
            K, steady_eqs, factor, y, static_cast<int>(step + 1));
        require(population.minCoeff() >= 0,
                context + "population has a negative component");
        const rcmc::Real sum_error = population.sum() - 1;
        require(-tolerance < sum_error && sum_error < tolerance,
                context + "population does not sum to one");
        require((population - single).cwiseAbs().maxCoeff() < tolerance,
                context + "history differs from a single population calculation");
    }
}

} // namespace

int main() {
    try {
        check_population_history(0, 1);
        // Relabel the same network to exercise a nonidentity permutation.
        check_population_history(2, 0);
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
