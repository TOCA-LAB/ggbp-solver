#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <string>

#include <argparse/argparse.hpp>
#include <nlohmann/json.hpp>
#include "ggbp/graph.hpp"
#include "ggbp/gurobi_status.hpp"
#include "ggbp/ggbp.hpp"

using std::vector;
using std::cout, std::endl;
namespace fs = std::filesystem;
using json = nlohmann::json;

using std::chrono::high_resolution_clock;
using std::chrono::duration;
using std::chrono::seconds;

int main(int argc, char* argv[]) {
    argparse::ArgumentParser program("GGBP-Naive");
    program.add_description("Use the Naive ILP to determine b_r(G)");

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

    ggbp::ilp::naive::LPParams params;
       program.add_argument("-V", "--verbose")
        .help("Increase output verbosity")
        .flag()
        .store_into(params.verbose);
        
    program.add_argument("-t", "--timeout")
        .help("Define the timeout (in seconds) for the LP solver")
        .default_value(7200.0)
        .store_into(params.time_limit);

    program.add_argument("-M", "--mem-lim")
        .help("Define the memory limit (in GB) for the LP solver")
        .default_value(14.0)
        .store_into(params.mem_limit);

    program.add_argument("-m", "--soft-mem-lim")
        .help("Define the soft memory (in GB) limit for the LP solver")
        .default_value(12.0)
        .store_into(params.soft_mem_limit);

    program.add_argument("-lr", "--linear-relax")
        .help("Run the linear relaxation of the LP")
        .flag()
        .store_into(params.linear_relax);

    program.add_argument("-opt", "--optimization")
        .help("Run the optmization LP")
        .flag()
        .store_into(params.opt);

    program.add_argument("-rc", "--opt-const-radius-cut")
        .help("Add the Radius Cut optimization constraint")
        .flag()
        .store_into(params.radius_cut);

    program.add_argument("-ac", "--opt-const-active-neighborhood")
        .help("Add the Active Neighborhood optimization constraint")
        .flag()
        .store_into(params.active_neighborhood);

    program.add_argument("-lb", "--lower-bound")
        .help("Define the lower bound for this instance")
        .scan<'i', int>();

    program.add_argument("-ub", "--upper-bound")
        .help("Define the upper bound for this instance")
        .scan<'i', int>();
    
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
        return 1;
    }

    std::string graph_name = fs::path(graph_file_path).stem().string();

    if (params.verbose) {
        cout << graph_name << ": " << graph.order() << " " << graph.size() << endl;
    }
    
    vector<int> burning_sequence;
    ggbp::ilp::GurobiStatus solve_status;
    
    int UB, LB;
    std::optional<double> lr;

    if (program.is_used("-lb")) {
        LB = program.get<int>("-lb");
    } else {
        LB = ggbp::ilp::get_lower_bound(graph, r, params.verbose);
    }

    if (program.is_used("-ub")) {
        UB = program.get<int>("-ub");
    } else {
        UB = ggbp::ilp::get_upper_bound(graph, r, params.verbose);
    }

    if (params.opt) {
        lr = ggbp::ilp::naive::solve_ilp(graph, r, UB, params, burning_sequence, solve_status);
    
        if (solve_status.error) {
            UB = solve_status.ObjVal;
            LB = solve_status.ObjBound;
        } else {
            UB = burning_sequence.size();

            if (lr.has_value()) {
                LB = lr.value();
                LB += (lr.value() - LB) > 0.001;
            } else {
                LB = burning_sequence.size();
            }
        }
    } else {
        ggbp::ilp::naive::find_burning_sequence(graph, r, LB, UB, params, burning_sequence, solve_status);
    }
    
    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;

    cout << "graph: " << graph_name << ", n = " << graph.order() << ", m = " << graph.size() << "\n"
         << "runtime = " << ms_double.count() / 1000 << endl;

    if (!solve_status.error) {     
        cout << "b(G) = " << burning_sequence.size() << "\n";

        cout << "burning Sequence = ";
        for (auto u : burning_sequence)
            cout << u + 1 << " ";
        cout << endl;
    } else {
        cout << "Failed to compute burning sequence. Error: " << solve_status.msg << endl;
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
        
        if (params.linear_relax)
            j["linear_relaxation"] = lr.value_or(UB);
        
        j["upper_bound"] = UB;
        j["lower_bound"] = LB;

        std::string alg = "dec";

        if (params.opt) {
            if (params.linear_relax)
                alg = "relax";
            else
                alg = "opt";
        }
        j["algorithm"] = "naive-" + alg;

        j["completed"] = !solve_status.error;

        if (solve_status.error)
            j["error_message"] = solve_status.msg;

        out_file << std::setw(4) << j << endl;
    }

    return 0;
}
