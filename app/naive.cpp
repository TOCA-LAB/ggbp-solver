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

	argparse::ArgumentParser program("Naive");

	std::string graph_file_path;
	program.add_argument("file_path")
		.store_into(graph_file_path);
	
	int r;
	program.add_argument("r")
		.help("Define the r constant for the GGBP")
		.required()
		.store_into(r);
        
    int verbose = 0;
	program.add_argument("-V", "--verbose")
		.help("Increase output verbosity")
  		.flag();

    program.add_argument("-o", "--out-file")
        .help("Choose the output file");

    int timeout;
    program.add_argument("-t", "--timeout")
        .help("Define the timeout for the LP solver")
        .default_value(7200)
        .store_into(timeout);
	
	try {
		program.parse_args(argc, argv);
	}
	catch (const std::exception& err) {
		std::cerr << err.what() << endl;
		std::cerr << program;
		std::exit(1);
	}
	
	auto t1 = high_resolution_clock::now();

	Graph graph = Graph::from_mtx_file(graph_file_path);

	std::string::iterator ptr, end_ptr;
	for (ptr = graph_file_path.end(); *ptr != '/'; --ptr) {
		if (*ptr == '.')
			end_ptr = ptr;
	}
	std::string graph_name(ptr+1, end_ptr);

	if (program["--verbose"] == true) {
		cout << graph_name << ": " << graph.order() << " " << graph.size() << endl;
		verbose = 1;
	}
	vector<int> burning_sequence;
	ggbp::ilp::GurobiErrorInfo err_info;
	err_info.msg = "";
	
	ggbp::ilp::naive::solve_ilp(graph, r, burning_sequence, err_info, verbose);

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
            out_file << burning_sequence.size() << "; " << burning_sequence.size() << "; " 
                     << burning_sequence.size() << "; ";
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
