#include <iostream>
#include <chrono>
#include <string>
#include <fstream>
#include "argparse/argparse.hpp"
#include "ggbp/ggbp.hpp"
#include "ggbp/graph.hpp"

using std::vector;
using std::cout, std::endl;

int main(int argc, char* argv[]) {
    using std::chrono::high_resolution_clock;
    using std::chrono::duration;
    using std::chrono::seconds;

    argparse::ArgumentParser program("RunBRKGA");

    std::string graph_file_path;
    program.add_argument("file_path")
        .store_into(graph_file_path);
    
    int r;
    program.add_argument("r")
        .help("Define the r constant for the GGBP")
        .required()
        .store_into(r);
    
    int T;
    program.add_argument("-T", "--iter")
        .help("Define the number of iterations of GRASP")
        .store_into(T)
        .default_value(100);

    program.add_argument("-V", "--verbose")
        .help("Increase output verbosity")
          .flag();

    program.add_argument("-d", "--debug")
        .help("Choose to run only some functions")
        .default_value("normal")
        .required();

        
    try {
        program.parse_args(argc, argv);
    }
    catch (const std::exception& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        std::exit(1);
    }

    auto t1 = high_resolution_clock::now();

    Graph graph(graph_file_path);

    std::string::iterator ptr, end_ptr;
    for (ptr = graph_file_path.end(); *ptr != '/'; --ptr) {
        if (*ptr == '.')
            end_ptr = ptr;
    }
    std::string graph_name(ptr+1, end_ptr);

    int verbose = program["--verbose"] == true ? 1 : 0;

    if (verbose)
        std::cout << graph_name << ": " << graph.order() << " " << graph.size() << std::endl;

    vector<int> burning_sequence;

    ggbp::grasp::grasp(graph, r, T, burning_sequence, verbose);

    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;

    std::cout << graph_name << "; "
              << graph.order() << "; "<< graph.size() << "; "
              << burning_sequence.size() << "; "
              << ms_double.count()/1000 << "; [";
    for (auto i = burning_sequence.begin(); i < burning_sequence.end(); ++i)
        std::cout << *i+1 << ","; // adding 1 to match the original vertex index 
    std::cout << "]" << std::endl;

    std::ofstream benchmark_csv("./benchmarks/runs/benchmark-r" + std::to_string(r) + ".csv", std::ofstream::app);
    if (!benchmark_csv)
        std::cout << "Erro ao abrir o arquivo .csv" << std::endl;
    else {	
        benchmark_csv << graph_name << "; "
                << graph.order() << "; "<< graph.size() << "; "
                << burning_sequence.size() << "; "
                << ms_double.count()/1000 << "; [";
        for (auto i = burning_sequence.begin(); i < burning_sequence.end(); ++i)
            benchmark_csv << *i+1 << ","; // adding 1 to match the original vertex index 
        benchmark_csv << "]" << std::endl;
    }
    return 0;
}
