#include <iostream>
#include <string>
#include <fstream>
#include <filesystem>
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
    argparse::ArgumentParser program("GGBP-Unburned");

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

    program.add_argument("-V", "--verbose")
        .help("Increase output verbosity")
        .flag();

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

    auto burning_sequence = ggbp::heuristics::unburned_neighbors(graph, r);

    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;

    cout << graph_name << "; "
              << graph.order() << "; "<< graph.size() << "; "
              << burning_sequence.size() << "; "
              << ms_double.count()/1000 << "; [";

    for (auto u : burning_sequence)
        cout << u + 1 << ","; // adding 1 to match the original vertex index 
    cout << "]" << endl;

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
        j["algorithm"] = "unburned";
        j["burning_number"] = burning_sequence.size();
        j["burning_sequence"] = burning_sequence;
        j["execution_time"] = ms_double.count()/1000;

        out_file << std::setw(4) << j << endl;
    }
    
    return 0;
}
