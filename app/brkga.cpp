#include <chrono>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <format>
#include <string>

#include <argparse/argparse.hpp>
#include <nlohmann/json.hpp>
#include "ggbp/ggbp.hpp"
#include "ggbp/graph.hpp"

using std::vector;
using std::cout, std::endl;
namespace fs = std::filesystem;
using json = nlohmann::json;

using std::chrono::high_resolution_clock;
using std::chrono::duration;
using std::chrono::seconds;

int main(int argc, char* argv[]) {
    argparse::ArgumentParser program("GGBP-BRKGA");

    program.add_argument("file_path");
    
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

    program.add_argument("-V", "--verbose")
        .help("Increase output verbosity")
        .flag();

    ggbp::brkga::BRKGAParams params;
    program.add_argument("-p", "--pop")
        .help("Define the size of each population")
        .default_value(100U)
        .required()
        .store_into(params.p);
    
    program.add_argument("-G", "--gens")
        .help("Define the number of generations")
        .default_value(100U)
        .required()
        .store_into(params.MAX_GENS);
    
    program.add_argument("-e", "--elite")
        .help("Define the fraction of population to be the elite-set")
        .default_value(0.20)
        .required()
        .store_into(params.pe);
    
    program.add_argument("-m", "--mutant")
        .help("Define the fraction of population to be replaced by mutants")
        .default_value(0.10)
        .required()
        .store_into(params.pm);

    program.add_argument("-r", "--rhoe")
        .help("Define the probability that offspring inherit an allele from elite parent")
        .default_value(0.70)
        .required()
        .store_into(params.rhoe);

    program.add_argument("-K", "--n_pops")
        .help("Define the number of independent populations")
        .default_value(2U)
        .required()
        .store_into(params.K);
    
    program.add_argument("-T", "--MAXT")
        .help("Define the number of threads for parallel decoding")
        .default_value(4U)
        .required()
        .store_into(params.MAXT);
    
    program.add_argument("-S", "--seed")
        .help("Define seed for the RNG");

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

    std::string graph_file_path = program.get("file_path");

    Graph graph;

    if (format == "mtx") {
        graph = Graph::from_mtx_file(graph_file_path);
    } else if (format == "g6") {
        std::ifstream g6_file(graph_file_path);
        std::string g6_string {std::istreambuf_iterator<char>(g6_file), std::istreambuf_iterator<char>() };
        graph = Graph(g6_string);
    } else {
        cout << "Invalid graph format" << endl;
        return 1;
    }
    
    std::string graph_name = fs::path(graph_file_path).stem().string();

    if (program.present("-S")) {
        params.seed = program.get<unsigned>("-S");
    } else {
        params.seed = std::chrono::system_clock::now().time_since_epoch().count();
    }

    vector<int> burning_sequence, iters;
    int first_best;

    ggbp::brkga::brkga(graph, r, params, burning_sequence, iters, first_best);


    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;

    cout << std::format("graph: {}, n = {}, m = {}\n", graph_name, graph.order(), graph.size())
         << std::format("runtime = {}", ms_double.count() / 1000) << endl
         << std::format("b_{}(G) = {}\n", r, burning_sequence.size());

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
        j["burning_number"] = burning_sequence.size();
        j["burning_sequence"] = burning_sequence;
        j["execution_time"] = ms_double.count()/1000;
        j["iters"] = iters;
        j["first_best"] = first_best;
        j["algorithm"] = "brkga";
        j["params"] = {
            {"population", params.p},
            {"generations", params.MAX_GENS},
            {"elite_fraction", params.pe},
            {"mutant_fraction", params.pm},
            {"alele_probability", params.rhoe},
            {"number_populations", params.K},
            {"seed", params.seed}
        };

        out_file << std::setw(4) << j << endl;
    }

    return 0;
}
