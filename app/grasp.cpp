#include <chrono>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>

#include <argparse/argparse.hpp>
#include <nlohmann/json.hpp>
#include "ggbp/ggbp.hpp"
#include "ggbp/graph.hpp"
#include "utils/utils.hpp"

using std::vector;
using std::cout, std::endl;
namespace fs = std::filesystem;
using json = nlohmann::json;

using std::chrono::high_resolution_clock;
using std::chrono::duration;
using std::chrono::seconds;

int main(int argc, char* argv[]) {
    argparse::ArgumentParser program("GGBP-GRASP");

    std::string graph_file_path;
    program.add_argument("file_path")
        .store_into(graph_file_path);
    
    int r;
    program.add_argument("r")
        .help("Define the r constant for the GGBP")
        .required()
        .store_into(r);
    
    std::string format;
    program.add_argument("--format")
        .help("Set the format for the instance. Choices: {mtx, g6}")
        .default_value(std::string{"mtx"})
        .choices("mtx", "g6")
        .store_into(format);

    int T;
    program.add_argument("-T", "--iter")
        .help("Define the number of iterations of GRASP")
        .default_value(100)
        .store_into(T);

    bool verbose;
    program.add_argument("-V", "--verbose")
        .help("Increase output verbosity")
        .flag()
        .store_into(verbose);

    program.add_argument("-S", "--seed")
        .help("Define seed for the RNG")
        .scan<'u', unsigned>();

    program.add_argument("--json")
        .help("Choose the output json file");

        
    try {
        program.parse_args(argc, argv);
    }
    catch (const std::exception& err) {
        std::cerr << err.what() << endl;
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
    
    std::string graph_name = fs::path(graph_file_path).stem().string();

    unsigned seed;

    if (program.is_used("-S")) {
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

    cout << "graph: " << graph_name << ", n = " << graph.order() << ", m = " << graph.size() << "\n"
         << "runtime = " << ms_double.count() / 1000 << endl
         << "b(G) = " << burning_sequence.size() << "\n";

    cout << "burning Sequence = ";
    for (auto u : burning_sequence)
        cout << u + 1 << " ";
    cout << endl;

    if (program.present("--json")) {
        std::ofstream out_file(program.get("--json"));
        
        if (!out_file) {
            cout << "Arquivo de saída inválido" << endl;
            
            return 1;
        }
        
        json j;

        j["instance"] = graph_name;
        j["order"] = graph.order();
        j["size"] = graph.size();
        j["r"] = r;
        j["burning_number"] = burning_sequence.size();
        j["burning_sequence"] = burning_sequence;
        j["execution_time"] = ms_double.count()/1000;
        j["iters"] = iters;
        j["first_best"] = first_best;
        j["algorithm"] = "grasp";
        j["params"] = {
            {"iterations", T},
            {"seed", seed}
        };

        out_file << std::setw(4) << j << endl;
    }
    
    return 0;
}
