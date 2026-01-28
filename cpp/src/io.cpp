#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <optional>
#include <regex>
#include <tuple>

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>

#include "chemistry.hpp"
#include "io.hpp"

namespace rcmc {
namespace fs = std::filesystem;

// Convert std::string to T.
template <class T> std::optional<T> from_string(const std::string &str) {
    try {
        return boost::lexical_cast<T>(str);
    } catch (...) {
        return std::nullopt;
    }
}

// Read one line from `in` and check whether the line matches to `pattern`. If
// it is, every `i`th submatch is converted to type of the `i`th template
// argument and is packed as a tuple. Returns std::nullopt when some operations
// fail.
template <class... Args>
std::optional<std::tuple<Args...>>
get_line_and_parse(std::istream &in, const std::regex &pattern) {
    std::string line;
    if (!std::getline(in, line)) {
        return std::nullopt;
    }

    std::smatch match;
    if (!std::regex_match(line, match, pattern)) {
        std::cout << "does not match" << std::endl;
        return std::nullopt;
    }

    if (match.size() < sizeof...(Args)) {
        return std::nullopt;
    }

    int i = 0;
    const std::tuple tuple = {
        [&]() { return from_string<Args>(match[++i].str()); }()...};

    if (std::apply([](auto... args) { return (args && ...); }, tuple)) {
        return std::apply([](auto... args) { return std::tuple{*args...}; },
                          tuple);
    } else {
        return std::nullopt;
    }
}

std::optional<UndirectedNetwork> load_MinPATH(const fs::path &path) {
    std::ifstream in(path);
    if (!in) {
        std::cerr << "Could not open " << path << std::endl;
        std::exit(1);
    }

    // List of EQ
    int num_eqs;
    const std::regex re_list_of_eqs(R"(List of EQs \((\d+)\):.*)");
    if (const auto tuple = get_line_and_parse<int>(in, re_list_of_eqs); tuple) {
        num_eqs = std::get<0>(*tuple);
    } else {
        return std::nullopt;
    }

    UndirectedNetwork G(num_eqs);

    const std::regex re_eq(R"( EQ +(\d+) \( *(-?\d*\.\d*)\).*)");
    for (int k = 0; k < num_eqs; ++k) {
        const auto tuple = get_line_and_parse<int, double>(in, re_eq);
        if (!tuple) {
            return std::nullopt;
        }
        const auto [i, v] = *tuple;
        if (i < 0 || num_eqs <= i) {
            return std::nullopt;
        }
        G.vertices[i] = v;
    }

    if (!get_line_and_parse(in, std::regex(""))) {
        return std::nullopt;
    }

    // List of TSs
    int num_tss;
    const std::regex re_list_of_tss(R"(List of TSs \((\d+)\):.*)");
    if (const auto tuple = get_line_and_parse<int>(in, re_list_of_tss); tuple) {
        num_tss = std::get<0>(*tuple);
    } else {
        return std::nullopt;
    }

    const std::regex re_ts(
        R"( TS +\d+: +(-?\d+) - +(-?\d+) \( *(-?\d*\.\d*)\).*)");
    for (int k = 0; k < num_tss; ++k) {
        const auto tuple = get_line_and_parse<int, int, double>(in, re_ts);
        if (!tuple) {
            std::cerr << "Failed to parse " << path << std::endl;
            std::exit(1);
        }

        const auto [i, j, v] = *tuple;
        if (0 <= i && i < num_eqs && 0 <= j && j < num_eqs && i != j) {
            G.set_value(i, j, v);
        }
    }

    if (!get_line_and_parse(in, std::regex(""))) {
        return std::nullopt;
    }

    // List of PTs
    int num_pts;
    const std::regex re_list_of_pts(R"(List of PTs \((\d+)\):.*)");
    if (const auto tuple = get_line_and_parse<int>(in, re_list_of_pts); tuple) {
        num_pts = std::get<0>(*tuple);
    } else {
        return std::nullopt;
    }

    const std::regex re_pt(
        R"( PT +\d+: +(-?\d+) - +(-?\d+) \( *(-?\d*\.\d*)\).*)");
    for (int k = 0; k < num_pts; ++k) {
        const auto tuple = get_line_and_parse<int, int, double>(in, re_pt);
        if (!tuple) {
            return std::nullopt;
        }

        const auto [i, j, v] = *tuple;
        if (0 <= i && i < num_eqs && 0 <= j && j < num_eqs && i != j) {
            G.set_value(i, j, v);
        }
    }

    return G;
}

std::optional<UndirectedNetwork> load_kINP(const fs::path &path,
                                           const Real temperature) {
    std::ifstream in(path);
    if (!in) {
        std::cerr << "Could not open " << path << std::endl;
        std::exit(1);
    }

    // List of EQ
    int num_eqs;
    const std::regex re_list_of_eqs(R"(List of EQs \((\d+)\):.*)");
    if (const auto tuple = get_line_and_parse<int>(in, re_list_of_eqs); tuple) {
        num_eqs = std::get<0>(*tuple);
    } else {
        return std::nullopt;
    }

    Vector eqs_pi(num_eqs);

    for (int k = 0; k < num_eqs; ++k) {
        const std::regex e(R"(.*\[([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?).*)");
        if (const auto tuple = get_line_and_parse<double>(in, e); tuple) {
            eqs_pi(k) = std::get<0>(*tuple);
        } else {
            return std::nullopt;
        }
    }

    // OS groups
    std::string line;
    do {
        if (!std::getline(in, line)) {
            return std::nullopt;
        }
    } while (!boost::algorithm::contains(line, "input EQs were grouped into"));

    const std::regex re_num_oss(
        R"(\d+ input EQs were grouped into (\d+) original states \(OSs\) by GPA)");
    std::smatch match;
    if (!std::regex_match(line, match, re_num_oss)) {
        return std::nullopt;
    }

    int num_oss;
    if (const auto result = from_string<int>(match[1].str()); result) {
        num_oss = *result;
    } else {
        return std::nullopt;
    }

    std::vector<std::vector<int>> groups(num_oss);

    std::vector<std::string> split;

    int k = -1;
    for (;;) {
        if (!std::getline(in, line)) {
            return std::nullopt;
        }
        boost::algorithm::trim(line);
        boost::algorithm::split(split, line, boost::algorithm::is_any_of(" "),
                                boost::algorithm::token_compress_on);

        if (split.empty() || split[0].empty()) {
            continue;
        } else if (split[0] == "Connectivity") {
            if (k < num_oss - 1) {
                return std::nullopt;
            } else {
                break;
            }
        } else if (split[0] == "OS") {
            ++k;
            for (std::size_t i = 3; i < split.size(); ++i) {
                if (const auto result = from_string<int>(split[i]); result) {
                    groups[k].push_back(*result);
                } else {
                    return std::nullopt;
                }
            }
        } else {
            for (const auto &s : split) {
                if (const auto result = from_string<int>(s); result) {
                    groups[k].push_back(*result);
                } else {
                    return std::nullopt;
                }
            }
        }
    }

    if (std::any_of(groups.begin(), groups.end(),
                    [](const auto &g) { return g.empty(); })) {
        std::cout << "There is an empty group." << std::endl;
        return std::nullopt;
    }

    // Rate constant matrix
    do {
        if (!std::getline(in, line)) {
            return std::nullopt;
        }
    } while (!boost::starts_with(line, "Rate constant matrix of OSs"));

    std::vector<Eigen::Triplet<double>> triples;

    for (int j = 0; j < num_oss; ++j) {
        int i = 0;
        for (int b = 0; b < (num_oss + 4) / 5; ++b) {
            if (!std::getline(in, line)) {
                return std::nullopt;
            }

            boost::algorithm::split(split, line.substr(14),
                                    boost::is_any_of(" "),
                                    boost::algorithm::token_compress_on);

            for (const auto &s : split) {
                if (s == "") {
                    continue;
                }
                double value = from_string<double>(s).value();
                if (value != 0.0) {
                    triples.emplace_back(i, j, value);
                }
                ++i;
            }
        }
    }

    SparseMatrix input_K(num_oss, num_oss);
    input_K.setFromTriplets(triples.begin(), triples.end());

    const auto K = compute_oss_rcm(input_K, eqs_pi, groups);
    const auto N = rcm_to_energy_network(K, temperature);

    return N;
}

UndirectedNetwork load_network(const fs::path &path, const Real temperature) {
    const std::string suffix = path.stem();
    if (boost::algorithm::ends_with(suffix, "MinPATH")) {
        if (const auto G = load_MinPATH(path); G) {
            return *G;
        } else {
            std::cerr << "Failed to parse " << path << std::endl;
            std::exit(1);
        }
    } else if (boost::algorithm::ends_with(suffix, "kINP")) {
        if (const auto G = load_kINP(path, temperature); G) {
            return *G;
        } else {
            std::cerr << "Failed to parse " << path << std::endl;
            std::exit(1);
        }
    } else {
        std::cerr << "Unknown file type: \"" << path << "\"" << std::endl;
        std::cerr << "File name must end with \"MinPATH\" or \"kINP\"."
                  << std::endl;
        std::exit(1);
    }
}

SparseMatrix load_OSs_rcm(const fs::path &path) {
    std::ifstream in(path);
    if (!in) {
        std::cerr << "Could not open " << path << std::endl;
        std::exit(1);
    }

    std::string line;

    do {
        std::getline(in, line);
    } while (!boost::starts_with(line, "Rate"));

    int n;
    std::sscanf(line.c_str(), "Rate constant matrix of OSs (%d by %*d)", &n);

    std::vector<Eigen::Triplet<double>> triples;
    std::vector<std::string> buf;

    for (int j = 0; j < n; ++j) {
        int i = 0;
        for (int b = 0; b < (n + 4) / 5; ++b) {
            std::getline(in, line);
            line += ' ';
            std::size_t next, last = 14;
            while ((next = line.find(' ', last)) != std::string::npos) {
                double value;
                std::sscanf(line.data() + last, "%lf", &value);
                if (value != 0.0) {
                    triples.emplace_back(i, j, value);
                }
                last = next + 1;
                ++i;
            }
        }
    }

    SparseMatrix K(n, n);
    K.setFromTriplets(triples.begin(), triples.end());
    // set_diagonals(K);

    return K;
}

void save_network(const UndirectedNetwork &N,
                  const std::filesystem::path &path) {
    fs::create_directory(path.parent_path());

    std::ofstream fout(path);
    if (!fout) {
        std::cerr << "File open error: " << path << std::endl;
        std::exit(EXIT_FAILURE);
    }
    fout << std::setprecision(17);

    const int n = N.num_vertices();
    fout << n << '\n';

    for (const Real &v : N.vertices) {
        fout << v << '\n';
    }

    fout << N.num_edges() << '\n';

    for (int i = 0; i < n; ++i) {
        for (const auto &[j, v] : N.edges[i]) {
            if (i < j) {
                fout << i << ' ' << j << ' ' << v << '\n';
            }
        }
    }
}

void save_network(const DirectedNetwork &N, const std::filesystem::path &path) {
    fs::create_directory(path.parent_path());

    std::ofstream fout(path);
    if (!fout) {
        std::cerr << "File open error: " << path << std::endl;
        std::exit(EXIT_FAILURE);
    }
    fout << std::setprecision(17);

    const int n = N.num_vertices();
    fout << n << '\n';

    for (const Real &v : N.vertices) {
        fout << v << '\n';
    }

    fout << N.num_edges() << '\n';

    for (int i = 0; i < n; ++i) {
        for (const auto &[j, v] : N.edges[i]) {
            fout << i << ' ' << j << ' ' << v << '\n';
        }
    }
}

void save_two_networks(const UndirectedNetwork &N1, const UndirectedNetwork &N2,
                       const std::filesystem::path &path) {
    fs::create_directory(path.parent_path());

    std::ofstream fout(path);
    if (!fout) {
        std::cerr << "File open error: " << path << std::endl;
        std::exit(EXIT_FAILURE);
    }
    fout << std::setprecision(17);

    const int n = N1.num_vertices();
    fout << n << '\n';

    for (int i = 0; i < n; ++i) {
        fout << N1.vertices[i] << ' ' << N2.vertices[i] << '\n';
    }

    fout << N2.num_edges() << '\n';

    for (int i = 0; i < n; ++i) {
        for (const auto &[j, v] : N1.edges[i]) {
            if (i < j) {
                fout << i << ' ' << j << ' ' << v << ' ' << N2.edges[i].at(j)
                     << '\n';
            }
        }
    }
}

void save_steady_eqs(const std::vector<int> &S, const fs::path &path) {
    fs::create_directory(path.parent_path());

    std::ofstream fout(path);
    if (!fout) {
        std::cerr << "File open error: " << path << std::endl;
        std::exit(EXIT_FAILURE);
    }

    for (int s : S) {
        fout << s << '\n';
    }
}

void save_population_history(const std::vector<Real> &times,
                             const std::vector<Vector> &populations,
                             const fs::path &path) {
    const int T(times.size());
    assert(int(populations.size()) == T);

    const int n(populations[0].size());
    assert(std::all_of(populations.begin(), populations.end(),
                       [n](const auto &p) { return int(p.size()) == n; }));

    fs::create_directory(path.parent_path());

    std::ofstream fout(path);
    if (!fout) {
        std::cerr << "File open error: " << path << std::endl;
        std::exit(EXIT_FAILURE);
    }
    fout << std::setprecision(17);

    fout << "time";
    for (int i = 0; i < n; ++i) {
        fout << "," << i;
    }
    fout << '\n';

    for (int t = 0; t < T; ++t) {
        fout << times[t];
        for (int i = 0; i < n; ++i) {
            fout << "," << populations[t](i);
        }
        fout << '\n';
    }
}

void save_improve_history(
    const std::vector<std::string> &improve_history_header,
    const std::vector<std::vector<double>> &improve_history,
    const fs::path &path) {
    fs::create_directory(path.parent_path());

    std::ofstream fout(path);
    if (!fout) {
        std::cerr << "File open error: " << path << std::endl;
        std::exit(EXIT_FAILURE);
    }
    fout << std::setprecision(17);

    const int n(improve_history[0].size());
    const int T(improve_history.size());

    for (int i = 0; i < n; ++i) {
        fout << improve_history_header[i] << ",";
    }
    fout << '\n';

    for (int t = 0; t < T; ++t) {
        for (int i = 0; i < n; ++i) {
            fout << improve_history[t][i] << ",";
        }
        fout << '\n';
    }
}

} // namespace rcmc
