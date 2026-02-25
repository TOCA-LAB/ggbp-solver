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

    // vector<int> teste = g.vertices();
    // auto combs = ggbp::utils::Combinations(teste, r);
    // unsigned long long i = 0;
    // for (auto s : combs) {
    //     cout << i << ": ";
    //     for (auto u : s)
    //         cout << u << " ";
    //     cout << endl;

    //     vector<int> tree = generate_covering_tree(g, r, s);

    //     for (int i = 0; i < g.order(); i++) {
    //         cout << i << ": " << tree[i] << endl;
    //     }
    //     cout << endl;
    //     i++;
    // }

    vector<int> burning_sequence;
    for (int i = 7; i <= g.order(); i++) {
        ggbp::ilp::solve_covering_ilp(g, r, i, burning_sequence, verbose);
        // ggbp_covering(g, r, burning_sequence, verbose);

        cout << i << " Burning Sequence: ";
        for (auto u : burning_sequence)
            cout << u + 1 << " ";
        cout << endl;
    }


    // Vertex u;
    // burning_sequence = {1, 3, 4};
    // cout << is_burning_sequence(g, r, burning_sequence, u) << endl;

    return 0;
}
