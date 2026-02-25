#include <iostream>
#include <chrono>
#include <string>
#include <fstream>
#include <cstdlib>   // std::exit
#include <stdexcept> // std::invalid_argument, std::out_of_range

#include "argparse/argparse.hpp"
#include "ggbp/graph.hpp"
#include "ggbp/ggbp.hpp"
#include "utils/utils.hpp"
#include "utils/statistics.hpp"
#include "utils/fix_c++_print.hpp"
#include "utils/logger.hpp"
#include "utils/contract.hpp"

using std::vector;
using std::string;
using argparse::ArgumentParser;
using ggbp::stats::ElapsedTimer;

void
configure_command_line(ArgumentParser& program) {
    program.add_argument("file_path")
           .help("File path of the input graph")
           .required();

    program.add_argument("r")
           .help("Define the r constant for the GGBP")
           .required();

    program.add_argument("--verbosity")
           .help("Set verbosity level");

    program.add_argument("-v", "--verbose")
           .help("Enable verbosity level 1")
           .flag();

}

int main(int argc, char* argv[]) {

    #if 0
    ArgumentParser config("Greedy");
    configure_command_line(config);

    try {
        config.parse_args(argc, argv);
    } catch (const std::exception& err) {
        std::cerr << "Error parsing command line: " << err.what() << std::endl;
        std::cerr << config;
        return EXIT_FAILURE;
    }

    string graph_fpath = config.get<string>("file_path");
    string r_str = config.get<string>("r");
    int r = 0;
    try {
        r = std::stoi(r_str);
    } catch (...) {
        std::cerr << "Invalid r value: " << r_str << std::endl;
        return EXIT_FAILURE;
    }

    int verbosity = 0;
    if (config.present("--verbosity")) {
        string v_str = config.get<string>("--verbosity");
        try {
            verbosity = std::stoi(v_str);
        } catch (...) {
            std::cerr << "Invalid verbosity value: " << v_str << std::endl;
            return EXIT_FAILURE;
        }
    }

    if (config.get<bool>("-v")) {
        verbosity = std::max(verbosity, 1);
    }
    logger.set_level(verbosity);

    logger.info(1) << "reading from file" << std::endl;
    Graph G = Graph::from_mtx_file(graph_fpath);
    
    ElapsedTimer greedy_timer("[timer]\t[greedy]");
    greedy_timer.start();
    vector<Vertex> burning_sequence = ggbp::heuristics::greedy(G, r);
    greedy_timer.stop();


    #endif

    
    Graph G("IhCGGC@?G");
    vector<Vertex> burning_sequence = ggbp::heuristics::greedy(G, 2);
    std::cout << "burning_number: " << burning_sequence.size() << std::endl;
    std::cout << "burning_sequence: " <<  burning_sequence << std::endl;
    //greedy_timer.print();

    return 0;
}
