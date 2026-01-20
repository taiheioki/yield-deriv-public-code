#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <optional>
#include <string>

#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>

#include "chemistry.hpp"
#include "io.hpp"
#include "numdef.hpp"
#include "rcmc.hpp"

void run(const std::string& name,
         const int initial,
         const std::optional<int> product_opt,
         const double time,
         const double temperature)
{
    std::cout << "Temperature = " << temperature << " K" << std::endl;

    // Input data
    const auto input_path     = rcmc::DATA_DIR / name;
    const auto energy_network = rcmc::load_network(input_path, temperature);
    const auto K = rcmc::rate_constant_matrix(energy_network, temperature);
    std::cout << "Loaded " << input_path << "." << std::endl;

    // Print data information
    const int n = K.num_eq();
    std::cout << "#EQ = " << n << std::endl;
    std::cout << "#TS = " << K.num_ts() << std::endl;
    std::cout << std::endl;

    // Check the input
    if(initial >= n) {
        std::cerr << "Error: the initial EQ is larger than or equal to #EQ."
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }

    std::cout << "The initial is EQ" << initial << std::endl;

    if(product_opt && *product_opt >= n) {
        std::cerr << "Error: the product EQ is larger than or equal to #EQ."
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Compute steady variables
    const rcmc::Real t_max = time;
    std::cout << "Computing steady EQs until " << time << " sec" << std::endl;
    boost::timer::cpu_timer timer1;
    auto [steady_eqs, qssa_times, L_SS_factor] =
        rcmc::compute_steady_eqs(K, t_max);
    std::cout << timer1.format();
    const int num_S = steady_eqs.size();
    std::cout << "#Steady EQs = " << num_S << std::endl;
    std::cout << std::endl;

    const auto initial_pop = rcmc::standard_vector(n, initial);

    // Determine the product
    std::cout << "Computing the population..." << std::endl;
    boost::timer::cpu_timer timer;
    const auto population = rcmc::compute_population(
        K, steady_eqs, L_SS_factor, initial_pop, num_S);

    std::cout << timer.format();

    int product;
    if(product_opt) {
        product = *product_opt;
        std::cout << "The product is EQ" << product;
    }
    else {
        population.maxCoeff(&product);
        std::cout << "The product is set to EQ" << product;
    }

    std::cout << " (population = " << std::setprecision(17)
              << population(product) << ")";
    const bool is_S = std::find(steady_eqs.begin(), steady_eqs.end(), product)
                      != steady_eqs.end();
    std::cout << ", which is " << (is_S ? "steady" : "transient") << ".\n"
              << std::endl;

    // Compute derivative
    std::cout << "Computing the derivative of product's population..."
              << std::endl;
    boost::timer::cpu_timer timer2;
    const auto [DL, Dpi] = rcmc::compute_population_derivative(
        K, steady_eqs, L_SS_factor, initial_pop, product);
    const auto derivative_network =
        rcmc::derivative_by_energy(DL, Dpi, K, 300.0);
    std::cout << timer2.format() << std::endl;

    // Save the derivative
    const auto output_path =
        rcmc::RESULT_DIR / (input_path.string() + "_deriv.txt");
    rcmc::save_two_networks(energy_network, derivative_network, output_path);
    std::cout << "Saved the derivative to " << output_path << "." << std::endl;
}

int main(const int argc, const char* const* const argv)
{
    namespace po = boost::program_options;

    std::string name;
    int initial;
    int product;
    double time;
    double temperature;

    po::options_description description("Options");

    // clang-format off
    description.add_options()
        ("name", po::value(&name), "Input file name")
        ("initial,i", po::value(&initial)->default_value(0), "The initial EQ")
        ("product,p", po::value(&product), "The product EQ")
        ("time", po::value(&time)->default_value(86400.0), "The reaction time [sec]")
        ("temperature", po::value(&temperature)->default_value(300.0), "The temperature [K]")
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

    if(vm.count("help") || !vm.count("name")) {
        std::cout << description << std::endl;
        return EXIT_SUCCESS;
    }

    const std::optional product_opt =
        vm.count("product") ? std::optional(product) : std::nullopt;
    run(name, initial, product_opt, time, temperature);
}
