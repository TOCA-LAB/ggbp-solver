#include <iostream>
#include <chrono>
#include <string>
#include <fstream>
#include "argparse/argparse.hpp"
#include "ggbp/ggbp.hpp"
#include "ggbp/graph.hpp"
#include "brkga-lib/MTRand.h"
#include "brkga-lib/BRKGA.h"

using std::vector;
using std::cout, std::endl;

int main(int argc, char* argv[]) {
    using std::chrono::high_resolution_clock;
    using std::chrono::duration;
    using std::chrono::seconds;

    argparse::ArgumentParser program("GGBP-BRKGA");

    std::string format;
    program.add_argument("format")
        .store_into(format);

    program.add_argument("file_path");
    
    int r;
    program.add_argument("r")
        .help("Define the r constant for the GGBP")
        .required()
        .store_into(r);

    program.add_argument("-V", "--verbose")
        .help("Increase output verbosity")
        .flag();

    program.add_argument("-d", "--debug")
        .help("Choose to run only some functions")
        .default_value("normal")
        .required();

    program.add_argument("-p", "--pop")
        .help("Define the size of each population")
        .default_value(100U)
        .required()
        .scan<'u', unsigned>();
    
    program.add_argument("-G", "--gens")
        .help("Define the number of generations")
        .default_value(100U)
        .required()
        .scan<'u', unsigned>();
    
    program.add_argument("-e", "--elite")
        .help("Define the fraction of population to be the elite-set")
        .default_value(0.20)
        .required()
        .scan<'f', double>();
    
    program.add_argument("-m", "--mutant")
        .help("Define the fraction of population to be replaced by mutants")
        .default_value(0.10)
        .required()
        .scan<'f', double>();

    program.add_argument("-r", "--rhoe")
        .help("Define the probability that offspring inherit an allele from elite parent")
        .default_value(0.70)
        .required()
        .scan<'f', double>();

    program.add_argument("-K", "--n_pops")
        .help("Define the number of independent populations")
        .default_value(2U)
        .required()
        .scan<'u', unsigned>();
    
    program.add_argument("-T", "--MAXT")
        .help("Define the number of threads for parallel decoding")
        .default_value(4U)
        .required()
        .scan<'u', unsigned>();
    
    program.add_argument("-S", "--seed")
        .help("Define seed for the RNG");

    program.add_argument("-o", "--out-file")
        .help("Choose the output file");
        
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
    }
    
    std::string::iterator ptr, end_ptr;
    for (ptr = graph_file_path.end(); *ptr != '/' && *ptr != '\\'; --ptr) {
        if (*ptr == '.')
            end_ptr = ptr;
    }
    std::string graph_name(ptr+1, end_ptr);

    const unsigned n = graph.order();                   // size of chromosomes
    const unsigned p = program.get<unsigned>("-p");     // size of population
    const double pe = program.get<double>("-e");        // fraction of population to be the elite-set
    const double pm = program.get<double>("-m");        // fraction of population to be replaced by mutants
    const double rhoe = program.get<double>("-r");      // probability that offspring inherit an allele from elite parent
    const unsigned K = program.get<unsigned>("-K");     // number of independent populations
    const unsigned MAXT = program.get<unsigned>("-T");  // number of threads for parallel decoding

    ggbp::brkga::Decoder decoder(graph, r);	// initialize the decoder

    // seed to the random number generator
    const long unsigned rngSeed = std::chrono::system_clock::now().time_since_epoch().count();
    
    // initialize the random number generator
    MTRand rng(rngSeed);
    
    // initialize the BRKGA-based heuristic
    BRKGA< ggbp::brkga::Decoder, MTRand > algorithm(n, p, pe, pm, rhoe, decoder, rng, K, MAXT);
    
    unsigned generation = 0;                                // current generation
    const unsigned X_INTVL = 10;                            // exchange best individuals at every 100 generations
    const unsigned X_NUMBER = 2;                            // exchange top 2 best
    const unsigned MAX_GENS = program.get<unsigned>("-G");  // run for 1000 gens

    vector<int> iters;
    iters.reserve(MAX_GENS);

    unsigned best = graph.order();
    unsigned first_best = 0;
    
    do {
        algorithm.evolve();	// evolve the population for one generation
        int fitness = algorithm.getBestFitness();
        
        if (program["--verbose"] == true) {
            cout << "Generation " << generation << " best solution found has objective value = "
                 << fitness << endl;
        }

        iters.push_back(fitness);

        if (fitness < best) {
            best = fitness;
            first_best = generation;
        }

        if((++generation) % X_INTVL == 0) {
            algorithm.exchangeElite(X_NUMBER);	// exchange top individuals
        }
    } while (generation < MAX_GENS && duration<double, std::milli>(high_resolution_clock::now() - t1).count()/1000 < 7200);
    
    vector<int> burning_sequence(decoder.getBurningSequence(algorithm.getBestChromosome()));

    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;

    cout << algorithm.getBestFitness() << endl
         << ms_double.count()/1000 << endl;

    if (program.present("-o")) {
        std::ofstream out_file(program.get("-o"));
        
        if (!out_file) {
            cout << "Arquivo de saída inválido" << endl;
            
            return 1;
        }

        out_file << graph_name << "\n"
                 << graph.order() << " " << graph.size() << "\n\n"
                 << "r=" << r << endl
                 << rngSeed << endl
                 << "p=" << p << endl
                 << "gens=" << MAX_GENS << endl
                 << "pe=" << pe << endl
                 << "pm=" << pm << endl
                 << "rhoe=" << rhoe << endl
                 << "K=" << K << endl
                 << "MAXT=" << MAXT << "\n\n";
        
        out_file << burning_sequence.size() << endl;

        for (auto u : burning_sequence)
            out_file << u << ", ";
        out_file << "\n\n";

        out_file << ms_double.count()/1000 << "\n\n";

        out_file << first_best << endl;
        for (auto i : iters)
            out_file << i << ", ";
        out_file << "\n\n";
    }

    return 0;
}
