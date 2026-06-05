#include "argparse/argparse.hpp"

#include <iostream>
#include <string>

#include "ggbp/ggbp.hpp"
#include "ggbp/graph.hpp"
#include "utils/combinations.hpp"


using namespace std;

int main(int argc, char *argv[]) {
    argparse::ArgumentParser program("GGBP-Cover");

    string graph_file_path;
    program.add_argument("file_path")
        .store_into(graph_file_path);
    
    int r;
    program.add_argument("r")
        .help("Define the r constant for the GGBP")
        .required()
        .store_into(r);

    program.add_argument("-o", "--out-file")
        .help("Choose the output file");

    ggbp::ilp::GurobiParams params;
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

    program.add_argument("-lb", "--lower-bound")
        .help("Define the lower bound for this instance")
        .scan<'i', int>();

    program.add_argument("-ub", "--upper-bound")
        .help("Define the upper bound for this instance")
        .scan<'i', int>();
    
    try {
        program.parse_args(argc, argv);
    }
    catch (const std::exception& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        std::exit(1);
    }
    
    Graph g = Graph::from_mtx_file(graph_file_path);

    int UB, LB;
    if (program.is_used("-lb")) {
        LB = program.get<int>("-lb");
    } else {
        LB = ggbp::ilp::get_lower_bound(g, r, params.verbose);
    }

    if (program.is_used("-ub")) {
        UB = program.get<int>("-ub");
    } else {
        UB = ggbp::ilp::get_upper_bound(g, r, params.verbose);
    }

    vector<Vertex> burning_sequence;
    
    ggbp::ilp::cover::find_burning_sequence(g, r, UB, LB, burning_sequence, params);

    cout << "Burning Sequence: ";
    for (auto u : burning_sequence)
        cout << u + 1 << " ";
    cout << endl;

    
    return 0;
}
