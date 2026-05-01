#include <iostream>
#include <chrono>
#include <random>
#include <string>
#include <fstream>
#include "argparse/argparse.hpp"
#include "ggbp/ggbp.hpp"
#include "ggbp/graph.hpp"

#include "utils/utils.hpp"
#include <print>

using std::vector;
using std::cout, std::endl;

int main(int argc, char* argv[]) {
    using std::chrono::high_resolution_clock;
    using std::chrono::duration;
    using std::chrono::seconds;

    argparse::ArgumentParser program("GGBP-GRASP");

    std::string format;
    program.add_argument("format")
        .store_into(format);

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
        .default_value(100)
        .store_into(T);

    program.add_argument("-V", "--verbose")
        .help("Increase output verbosity")
        .flag();

    program.add_argument("-S", "--seed")
        .help("Define seed for the RNG")
        .scan<'u', unsigned>();

    program.add_argument("-o", "--out-file")
        .help("Choose the output file");

        
    try {
        program.parse_args(argc, argv);
    }
    catch (const std::exception& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        std::exit(1);
    }

    auto t1 = high_resolution_clock::now();

    Graph graph;

    if (format == "mtx") {
        graph = Graph::from_mtx_file(graph_file_path);
    } else if (format == "g6") {
        std::ifstream g6_file(graph_file_path);
        std::string g6_string {std::istreambuf_iterator<char>(g6_file), std::istreambuf_iterator<char>() };
        graph = Graph(g6_string);
    } else {
        cout << "Invalid graph format" << endl;
    }
    
    std::string::iterator ptr, end_ptr;
    for (ptr = graph_file_path.end(); *ptr != '\\' && *ptr != '/'; --ptr) {
        if (*ptr == '.')
            end_ptr = ptr;
    }
    std::string graph_name(ptr+1, end_ptr);

    int verbose = program["--verbose"] == true ? 1 : 0;

    unsigned seed;

    if (program.present("-S")) {
        seed = program.get<unsigned>("-S");
    } else {
        seed = std::chrono::system_clock::now().time_since_epoch().count();
    }

    static std::default_random_engine e(seed);

    vector<int> burning_sequence;
    vector<int> iters;
    int first_best;
    
    ggbp::grasp::grasp(graph, r, T, e, burning_sequence, iters, first_best, verbose);

    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;

    std::cout << graph_name << "; "
              << graph.order() << "; "<< graph.size() << "; "
              << burning_sequence.size() << "; "
              << ms_double.count()/1000 << "; [";
    for (auto i = burning_sequence.begin(); i < burning_sequence.end(); ++i)
        std::cout << *i+1 << ","; // adding 1 to match the original vertex index 
    std::cout << "]" << std::endl;

    if (program.present("-o")) {
        std::ofstream out_file(program.get("-o"));
        
        if (!out_file) {
            cout << "Arquivo de saída inválido" << endl;
            
            return 1;
        }

        out_file << graph_name << "\n"
                 << graph.order() << " " << graph.size() << "\n\n"
                 << "r=" << r << endl
                 << seed << endl
                 << "Iters: " << T << "\n\n";
        
        out_file << burning_sequence.size() << endl;

        for (auto u : burning_sequence)
            out_file << u << ", ";
        out_file << "\n\n";

        out_file << ms_double.count()/1000 << "\n\n";

        out_file << first_best << endl;
        for (auto i : iters)
            out_file << i << ", ";
        out_file << "\n\n";
    }
    return 0;
}
