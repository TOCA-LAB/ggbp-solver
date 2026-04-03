#include <iostream>
#include <chrono>
#include <string>
#include <fstream>
#include "argparse/argparse.hpp"
#include "ggbp/graph.hpp"
#include "ggbp/ggbp.hpp"

using std::vector;
using std::cout, std::endl;

int main(int argc, char* argv[]) {
    using std::chrono::high_resolution_clock;
    using std::chrono::duration;
    using std::chrono::seconds;

    argparse::ArgumentParser program("GGBP-Naive");

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

    program.add_argument("-o", "--out-file")
        .help("Choose the output file");

    ggbp::ilp::naive::LPParams params;
       program.add_argument("-V", "--verbose")
        .help("Increase output verbosity")
        .flag()
        .store_into(params.verbose);
        
    program.add_argument("-t", "--timeout")
        .help("Define the timeout for the LP solver")
        .default_value(7200.0)
        .store_into(params.time_limit);

    program.add_argument("-M", "--mem-lim")
        .help("Define the memory limit for the LP solver")
        .default_value(14.0)
        .store_into(params.mem_limit);

    program.add_argument("-m", "--soft-mem-lim")
        .help("Define the soft memory limit for the LP solver")
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
        .help("Define the timeout for the LP solver")
        .flag()
        .store_into(params.radius_cut);

    program.add_argument("-ac", "--opt-const-active-neighborhood")
        .help("Define the timeout for the LP solver")
        .flag()
        .store_into(params.active_neighborhood);
    
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

    std::string::iterator ptr, end_ptr;
    for (ptr = graph_file_path.end(); (*ptr != '/') && (*ptr != '\\'); --ptr) {
        if (*ptr == '.')
            end_ptr = ptr;
    }
    std::string graph_name(ptr+1, end_ptr);

    if (params.verbose) {
        cout << graph_name << ": " << graph.order() << " " << graph.size() << endl;
    }
    
    vector<int> burning_sequence;
    ggbp::ilp::GurobiErrorInfo err_info;
    err_info.msg = "";
    
    if (params.opt) {
        int B = ggbp::ilp::naive::get_upper_bound(graph, r, params.verbose);
        ggbp::ilp::naive::solve_ilp(graph, r, B, burning_sequence, err_info, params);
    } else {
        ggbp::ilp::naive::find_burning_sequence(graph, r, burning_sequence, params);
    }


    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;

    cout << graph_name << "; "
         << graph.order() << "; "<< graph.size() << "; ";

    if (err_info.msg == "") {
        cout << burning_sequence.size() << "; " << burning_sequence.size() << "; " 
             << burning_sequence.size() << "; ";

        for (auto u : burning_sequence)
            cout << u << ", ";
        cout << endl;

    } else {
        cout << err_info.ObjBound << "; " << err_info.ObjVal << "; --; ";
    }

    cout << ms_double.count()/1000 << " " + err_info.msg << endl;

    if (program.present("-o")) {
        std::ofstream out_file(program.get("-o"));
        
        if (!out_file)
            cout << "Arquivo de saída inválido" << endl;
        
        out_file << graph_name << "; "
                 << graph.order() << "; "<< graph.size() << "; ";

        if (err_info.msg == "") {
            out_file << burning_sequence.size() << "; ";
        } else {
            out_file << err_info.ObjBound << "; " << err_info.ObjVal << "; --; ";
        }

        out_file << ms_double.count()/1000 << " " + err_info.msg << endl ;
    }

    // for (auto i = burning_sequence.begin(); i < burning_sequence.end(); ++i)
    // 	benchmark_csv << *i+1 << ","; // adding 1 to match the original vertex index 
    // benchmark_csv << "]" << endl;

    return 0;
}
