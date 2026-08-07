#include <chrono>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <optional>
#include <string>

#include <argparse/argparse.hpp>
#include <nlohmann/json.hpp>
#include "ggbp/graph.hpp"
#include "ggbp/ggbp.hpp"
#include "ggbp/gurobi_status.hpp"

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

    if (params.linear_relax && !params.opt) {
        std::cerr
            << "--linear-relax requires --optimization because an LP "
               "relaxation cannot certify feasibility in the binary search"
            << endl;
        return 2;
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

    if (params.verbose) {
        cout << graph_name << ": " << graph.order() << " " << graph.size() << endl;
    }
    
    vector<int> burning_sequence;
    ggbp::ilp::GurobiErrorInfo err_info;
    std::optional<double> relaxation_bound;
    
    int UB, LB;
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
        relaxation_bound = ggbp::ilp::naive::solve_ilp(
            graph, r, UB, params, burning_sequence, err_info);
    
        if (err_info.error) {
            if (err_info.solution_count > 0) {
                UB = err_info.ObjVal;
                LB = err_info.ObjBound;
            }
        } else if (params.linear_relax && !relaxation_bound) {
            err_info.error = true;
            err_info.msg = "solver returned no LP relaxation value";
        } else if (!params.linear_relax) {
            UB = burning_sequence.size();
            LB = burning_sequence.size();
        }
    } else {
        ggbp::ilp::naive::find_burning_sequence(graph, r, LB, UB, params, burning_sequence, err_info);
    }
    
    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;

    cout << "graph: " << graph_name << ", n = " << graph.order() << ", m = " << graph.size() << "\n"
         << "runtime = " << ms_double.count() / 1000 << endl;

    if (err_info.error) {
        cout << "Failed to compute burning sequence. Error: " << err_info.msg << endl;
    } else if (relaxation_bound) {
        cout << "LP relaxation lower bound = " << *relaxation_bound << "\n";
    } else {
        cout << "b(G) = " << burning_sequence.size() << "\n";

        cout << "burning Sequence = ";
        for (auto u : burning_sequence)
            cout << u + 1 << " ";
        cout << endl;
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
        j["execution_time"] = ms_double.count()/1000;
        j["algorithm"] = "naive";
        j["gurobi_status"] = err_info.native_status;
        j["gurobi_status_name"] =
            ggbp::ilp::gurobi_status_name(err_info.native_status);
        j["gurobi_solution_count"] = err_info.solution_count;
        j["gurobi_runtime"] = err_info.runtime_seconds;
        j["gurobi_mip_gap"] = err_info.mip_gap
            ? json(*err_info.mip_gap)
            : json(nullptr);

        if (params.linear_relax) {
            j["solution_kind"] = "lp_relaxation_bound";
            j["upper_bound"] = UB;
            if (relaxation_bound && !err_info.error) {
                j["relaxation_bound"] = *relaxation_bound;
                j["lower_bound"] = *relaxation_bound;
            }
        } else {
            j["solution_kind"] = "integer_sequence";
            j["burning_number"] = burning_sequence.size();
            j["burning_sequence"] = burning_sequence;
            j["upper_bound"] = UB;
            j["lower_bound"] = LB;
        }

        j["completed"] = !err_info.error;

        if (err_info.error)
            j["error_message"] = err_info.msg;

        out_file << std::setw(4) << j << endl;
    }

    return err_info.error ? 1 : 0;
}
