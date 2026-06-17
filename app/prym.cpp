#include <iostream>
#include <chrono>
#include <string>
#include <fstream>
#include <filesystem>
#include <argparse/argparse.hpp>
#include <nlohmann/json.hpp>
#include "ggbp/graph.hpp"
#include "ggbp/ggbp.hpp"

using std::vector;
using std::cout, std::endl;
namespace fs = std::filesystem;
using json = nlohmann::json;

int main(int argc, char* argv[]) {
    using std::chrono::high_resolution_clock;
    using std::chrono::duration;
    using std::chrono::seconds;

    argparse::ArgumentParser program("GBP-PRYM");
        
    std::string graph_file_path;
    program.add_argument("file_path")
        .store_into(graph_file_path);

    std::string format;
    program.add_argument("--format")
        .help("Set the format for the instance. Choices: {mtx, g6}")
        .default_value(std::string{"mtx"})
        .choices("mtx", "g6")
        .store_into(format);

    program.add_argument("-V", "--verbose")
        .help("Increase output verbosity")
        .flag();

    program.add_argument("-d", "--debug")
        .help("Choose to run only some functions")
        .default_value("normal")
        .required();
        
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

    std::string graph_name = fs::path(graph_file_path).stem();

    int grb_verbose = 0;
    if (program["--verbose"] == true) {
        std::cout << graph_name << ": " << graph.order() << " " << graph.size() << std::endl;
        grb_verbose = 1;
    }
    
    vector<int> burning_sequence;
    ggbp::ilp::prym::run_prym(graph, burning_sequence);

    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;

    std::cout << graph_name << "; "
              << graph.order() << "; "<< graph.size() << "; "
              << burning_sequence.size() << "; "
              << ms_double.count()/1000 << "; [";
    for (auto i = burning_sequence.begin(); i < burning_sequence.end(); ++i)
        std::cout << *i+1 << ","; // adding 1 to match the original vertex index 
    std::cout << "]" << std::endl;

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
        j["burning_number"] = burning_sequence.size();
        j["burning_sequence"] = burning_sequence;

        out_file << std::setw(4) << j << endl;
    }

    return 0;
}
