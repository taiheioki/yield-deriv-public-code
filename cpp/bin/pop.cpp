#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>

#include "chemistry.hpp"
#include "io.hpp"
#include "rcmc.hpp"

void run(const std::string &name, const int initial,
         const std::optional<double> t_max, const double temperature,
         const bool last_only, const std::filesystem::path &output_path) {
    // Output formatting
    // std::cout.imbue(std::locale(""));
    std::cout << std::setprecision(17);

    // Input data
    const auto input_path = rcmc::DATA_DIR / name;
    const auto energy_network = rcmc::load_network(input_path, temperature);
    const auto K = rcmc::rate_constant_matrix(energy_network, temperature);

    std::cout << energy_network.get_value(20, 2) << std::endl;
    std::cout << energy_network.get_value(12, 2) << std::endl;
    std::cout << energy_network.get_value(20, 12) << std::endl;
    std::cout << "Loaded " << input_path << "." << std::endl;

    // Print data information
    const int n = K.num_eq();
    std::cout << "#EQ = " << n << std::endl;
    std::cout << "#TS = " << K.num_ts() << std::endl;
    std::cout << std::endl;

    // Check the input
    if (initial >= n) {
        std::cerr << "Error: the initial EQ is larger than or equal to #EQ."
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }
    std::cout << "The initial is EQ" << initial << std::endl;

    // Compute steady variables
    std::cout << "Computing steady EQs until " << t_max.value_or(INFINITY)
              << " sec..." << std::endl;

    using clock = std::chrono::steady_clock;
    auto t0 = clock::now();
    auto [steady_eqs, qssa_times, L_SS_factor] =
        rcmc::compute_steady_eqs(K, t_max);
    const auto num_S = steady_eqs.size();
    const auto nnz_B = rcmc::nnz(L_SS_factor.B());
    const auto num_B = num_S * (num_S + 1) / 2;
    auto t1 = clock::now();
    auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
    std::cout << "elapsed = " << dt.count() << " ms\n";

    std::cout << "#Steady EQs = " << num_S << std::endl;
    std::cout << "#Nonzeros of B_QQ = " << nnz_B << " / " << num_B << " ("
              << std::setprecision(3) << double(nnz_B) / num_B * 100.0 << "%)"
              << std::endl;
    std::cout << std::endl;

    // Compute population
    std::cout << "Computing populations..." << std::endl;
    boost::timer::cpu_timer timer2;

    const auto initial_pop = rcmc::standard_vector(n, initial);
    std::vector<rcmc::Vector> pop_history;

    if (last_only) {
        pop_history = {rcmc::compute_population(K, steady_eqs, L_SS_factor,
                                                initial_pop, num_S)};
        qssa_times = {qssa_times.back()};
    } else {
        pop_history = rcmc::compute_population_history(
            K, steady_eqs, qssa_times, L_SS_factor, initial_pop, t_max);

        // Insert the initial population
        qssa_times.insert(qssa_times.begin(), 0.0);
        pop_history.insert(pop_history.begin(), initial_pop);
    }

    std::cout << timer2.format();
    std::cout << "Last population sum = " << pop_history.rbegin()->sum()
              << std::endl;
    std::cout << std::endl;

    // Save the population history
    rcmc::save_population_history(qssa_times, pop_history, output_path);
    std::cout << "Saved population history to " << output_path << "."
              << std::endl;
}

int main(const int argc, const char *const *const argv) {
    namespace po = boost::program_options;

    std::string name;
    std::string output_path;
    int initial;
    double temperature;
    bool last_only;

    po::options_description description("Options");

    // clang-format off
    description.add_options()
        ("name", po::value(&name), "Input file name")
        ("output,o", po::value(&output_path), "Output file path")
        ("initial,i", po::value(&initial)->default_value(0), "The initial EQ")
        ("time", po::value<double>(), "The reaction time [sec]")
        ("temperature", po::value(&temperature)->default_value(300.0), "The temperature [K]")
        ("last-only,l", po::value(&last_only), "Compute last population only")
        ("help,h", "Print this help message")
    ;
    // clang-format on

    po::positional_options_description positional_description;
    positional_description.add("name", -1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv)
                  .options(description)
                  .positional(positional_description)
                  .run(),
              vm);
    po::notify(vm);

    if (vm.count("help") || !vm.count("name")) {
        std::cout << description << std::endl;
        return EXIT_SUCCESS;
    }

    std::optional<double> time = std::nullopt;
    if (vm.count("time")) {
        time = vm["time"].as<double>();
    }

    if (output_path == "") {
        output_path = rcmc::RESULT_DIR / (name + "_pop.txt");
    }

    run(name, initial, time, temperature, last_only, output_path);
}