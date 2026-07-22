#include <chrono>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <string>

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

    ggbp::ilp::GurobiParams params;
       program.add_argument("-V", "--verbose")
        .help("Increase output verbosity")
        .flag()
        .store_into(params.verbose);
        
    program.add_argument("-t", "--timeout")
        .help("Define the timeout for the LP solver")
        .default_value(1800.0)
        .store_into(params.time_limit);

    program.add_argument("-M", "--mem-lim")
        .help("Define the memory limit for the LP solver")
        .default_value(14.0)
        .store_into(params.mem_limit);

    program.add_argument("-m", "--soft-mem-lim")
        .help("Define the soft memory limit for the LP solver")
        .default_value(12.0)
        .store_into(params.soft_mem_limit);
        
    program.add_argument("--json")
        .help("Choose the output json file");

        
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

    std::string graph_name = fs::path(graph_file_path).stem().string();
    
    vector<int> burning_sequence;
    ggbp::ilp::GurobiErrorInfo err_info;
    ggbp::ilp::prym::run_prym(graph, burning_sequence, params, err_info);

    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;

    cout << "graph: " << graph_name << ", n = " << graph.order() << ", m = " << graph.size() << "\n"
         << "runtime = " << ms_double.count() / 1000 << endl;

    if (!err_info.error) {     
        cout << "b(G) = " << burning_sequence.size() << "\n";

        cout << "burning Sequence = ";
        for (auto u : burning_sequence)
            cout << u + 1 << " ";
        cout << endl;
    } else {
        cout << "Failed to compute burning sequence. Error: " << err_info.msg << endl;
    }

    
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
        j["execution_time"] = ms_double.count()/1000;
        j["algorithm"] = "prym";

        j["completed"] = !err_info.error;

        if (err_info.error) {
            j["error_message"] = err_info.msg;
            j["upper_bound"] = err_info.ObjVal;
            j["lower_bound"] = err_info.ObjBound;
        }

        out_file << std::setw(4) << j << endl;
    }

    return 0;
}
