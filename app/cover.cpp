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

    bool verbose;
    program.add_argument("-V", "--verbose")
        .help("Increase output verbosity")
          .flag()
        .default_value(false)
        .store_into(verbose);
        
    try {
        program.parse_args(argc, argv);
    }
    catch (const std::exception& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        std::exit(1);
    }
    
    Graph g = Graph::from_mtx_file(graph_file_path);

    vector<Vertex> burning_sequence;
    
    ggbp::ilp::cover::ggbp_covering(g, r, burning_sequence, verbose);

    cout << "Burning Sequence: ";
    for (auto u : burning_sequence)
        cout << u + 1 << " ";
    cout << endl;

    
    return 0;
}
