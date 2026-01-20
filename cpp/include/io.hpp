#ifndef RCMC_IO_HPP
#define RCMC_IO_HPP

#include <filesystem>
#include <vector>

#include "network.hpp"
#include "numdef.hpp"

namespace rcmc {

const std::filesystem::path DATA_DIR = "data";
const std::filesystem::path RESULT_DIR = "result";

UndirectedNetwork load_network(const std::filesystem::path &path,
                               Real temperature);

SparseMatrix load_OSs_rcm(const std::filesystem::path &path);

void save_network(const UndirectedNetwork &N,
                  const std::filesystem::path &path);

void save_network(const DirectedNetwork &N, const std::filesystem::path &path);

void save_two_networks(const UndirectedNetwork &N1, const UndirectedNetwork &N2,
                       const std::filesystem::path &path);

void save_steady_eqs(const std::vector<int> &S,
                     const std::filesystem::path &path);

void save_population_history(const std::vector<Real> &times,
                             const std::vector<Vector> &populations,
                             const std::filesystem::path &path);

void save_improve_history(
    const std::vector<std::string> &improve_history_header,
    const std::vector<std::vector<double>> &improve_history,
    const std::filesystem::path &path);
} // namespace rcmc

#endif
