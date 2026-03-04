#include <cmath>
#include "gurobi_c++.h"
#include "utils/contract.hpp"
#include "brkga-lib/MTRand.h"
#include "brkga-lib/BRKGA.h"
#include "ggbp/graph.hpp"
#include "ggbp/decoder.hpp"

#include "ggbp/ggbp.hpp"

namespace ggbp::ilp::naive {
using namespace std;

int get_upper_bound(Graph const &graph, int const r, int verbose) {
    const unsigned n = graph.order();
    const unsigned p = 60;
    const unsigned MAX_GENS = 25;
    const double pe = 0.20;
    const double pm = 0.10;
    const double rhoe = 0.70;	
    const unsigned K = 2;
    const unsigned MAXT = 2;

    Decoder decoder(graph, r);

    const long unsigned rngSeed = 0;
    MTRand rng(rngSeed);

	BRKGA< Decoder, MTRand > algorithm(n, p, pe, pm, rhoe, decoder, rng, K, MAXT);
	unsigned generation = 0;
	const unsigned X_INTVL = 10;
	const unsigned X_NUMBER = 2;
	do {
		algorithm.evolve();	// evolve the population for one generation
		if (verbose) {
			std::cout << "Generation " << generation << " best solution found has objective value = "
					<< algorithm.getBestFitness() << std::endl;
		}

		if((++generation) % X_INTVL == 0) {
			algorithm.exchangeElite(X_NUMBER);	// exchange top individuals
		}
	} while (generation < MAX_GENS);

    return algorithm.getBestFitness();
}

void solve_ilp(Graph const &graph, int const r, vector<int> &burning_sequence, GurobiErrorInfo &err_info, bool verbose) {
    int const n_vertices = graph.order();
    int const B = get_upper_bound(graph, r, false);
    

    try {        
        GRBEnv env(true);
        env.set(GRB_IntParam_OutputFlag, verbose);
        env.start();
        GRBModel model(env);

        model.set(GRB_DoubleParam_TimeLimit, 7200.0);
        model.set(GRB_DoubleParam_SoftMemLimit, 13.0);
        model.set(GRB_DoubleParam_MemLimit, 13.5);

        // Defining the variables
        vector<vector<GRBVar>> vars_s(n_vertices, vector<GRBVar>(B));
        vector<vector<GRBVar>> vars_q(n_vertices, vector<GRBVar>(B));
        vector<GRBVar> vars_z(B);

        for (int i = 0; i < n_vertices; i++) {
            for (int j = 0; j < B; j++) {
                vars_s[i][j] = model.addVar(0, 1, 0, GRB_BINARY, "s_" + to_string(i) + "_" + to_string(j));
                vars_q[i][j] = model.addVar(0, 1, 0, GRB_BINARY, "q_" + to_string(i) + "_" + to_string(j));
            }
        }

        for (int i = 0; i < B; i++)
            vars_z[i] = model.addVar(0, 1, 0, GRB_BINARY, "z_" + to_string(i)); 
        
        // Adding the constraints

        // One ignitor vertex for each used turn
        for (int i = 0; i < B; i++) {
            GRBLinExpr expr = 0;
            for (int j = 0; j < n_vertices; j++)
                expr += vars_s[j][i];
            model.addConstr(expr == vars_z[i], "MIT_" + to_string(i)).set(GRB_IntAttr_Lazy, 2);
        }

        // A vertex can only be ignitor once
        for (int i = 0; i < n_vertices; i++) {
            GRBLinExpr expr = 0;
            for (int j = 0; j < B; j++)
                expr += vars_s[i][j];
            model.addConstr(expr <= 1, "I_" + to_string(i)).set(GRB_IntAttr_Lazy, 2);
        }

        // An ignitor also burns
        for (int i = 0; i < n_vertices; i++)
            for (int j = 0; j < B; j++)
                model.addConstr(vars_s[i][j] <= vars_q[i][j], "QI_" + to_string(i) + "_" + to_string(j));

        // A vertex burns only if there are r burning neighbours or if it is a source 
        for (int i = 0; i < n_vertices; i++)
            for (int j = 0; j < B-1; j++) {
                GRBLinExpr expr = 0, e_expr, ne_expr;
                for (auto neighbor : graph.neighbors(i))
                    expr += vars_q[neighbor][j];
                e_expr = (expr + 1 - r) * 1.0/graph.degree(i);
                ne_expr = expr * 1.0/r;
                model.addConstr(e_expr <= vars_q[i][j+1], 
                "ESP_"  + to_string(i) + "_" + to_string(j));
                model.addConstr(vars_q[i][j+1] <= (ne_expr + vars_s[i][j+1] + vars_q[i][j]), 
                "NESP_"  + to_string(i) + "_" + to_string(j));  
            }

        for (int u = 0; u < n_vertices; u++) {
            if (graph.degree(u) < r) {
                GRBLinExpr expr = 0;
                // cout << "Vertex: " << u << ", d(u) = " << graph.degree(u) << endl; 
                for (int i = 0; i < B; i++)
                    expr += vars_s[u][i];
                model.addConstr(expr == 1, "S_" + to_string(u));
            }
        }
        
        // A burned vertex continues burning throughout the entire process
        for (int i = 0; i < n_vertices; i++)
            for (int j = 0; j < B-1; j++)
                model.addConstr(vars_q[i][j] <= vars_q[i][j+1], "Q_" + to_string(i) + "_" + to_string(j));

        // On the first turn, all vertices are not burned, except the first source
        for (int i = 0; i < n_vertices; i++)
            model.addConstr(vars_q[i][0] == vars_s[i][0], "T0_" + to_string(i));

        // The first turn always happens
        model.addConstr(vars_z[0] == 1);

        // If there's a vertex not burning, a new turn is needed
        for (int i = 0; i < n_vertices; i++)
            for (int j = 0; j < B-1; j++)
                model.addConstr((1 - vars_q[i][j]) <= vars_z[j+1], "NT_" + to_string(i) + "_" + to_string(j));

        // All vertices must be burned at the end
        GRBLinExpr expr = 0;
        for (int i = 0; i < n_vertices; i++)
            expr += vars_q[i][B-1];
        model.addConstr(expr == n_vertices, "T_K");

        // Minimize the number of turns
        GRBLinExpr obj = 0;
        for (int i = 0; i < B; i++)
            obj += vars_z[i];
        model.setObjective(obj, GRB_MINIMIZE);

        for (int i = 0; i < n_vertices; i++)
            for (int j = 0; j < B-1; j++) {
                // vars_s[i][j].set(GRB_IntAttr_BranchPriority, B-j);
                // vars_z[j].set(GRB_IntAttr_BranchPriority, B-j);
            } 

        model.optimize();

        int opt_status = model.get(GRB_IntAttr_Status);

        if (opt_status == GRB_OPTIMAL) {
			for(int i = 0; i < B; i++) {
            	for(int j = 0; j < n_vertices; j++) {
					if(vars_s[j][i].get(GRB_DoubleAttr_X) > 0.5)
						burning_sequence.push_back(j);
                }
            }
            // printVector(burning_sequence, "Burning sequence");
        } else if (opt_status == GRB_TIME_LIMIT) {
            if (model.get(GRB_DoubleAttr_ObjVal) < B)
                err_info.ObjVal = model.get(GRB_DoubleAttr_ObjVal);
            else
                err_info.ObjVal = B;
            
            err_info.ObjBound = model.get(GRB_DoubleAttr_ObjBound);
            err_info.msg = "//TIME LIMIT";
        } else if (opt_status == GRB_MEM_LIMIT) {
            if (model.get(GRB_DoubleAttr_ObjVal) < B)
                err_info.ObjVal = model.get(GRB_DoubleAttr_ObjVal);
            else
                err_info.ObjVal = B;

            err_info.ObjBound = model.get(GRB_DoubleAttr_ObjBound);
            err_info.msg = "//OUT OF MEMORY (soft)";
        }
    } catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;

        if (e.getErrorCode() == GRB_ERROR_OUT_OF_MEMORY) {
            err_info.ObjVal = B;
            err_info.ObjBound = 0;
            err_info.msg = "//OUT OF MEMORY";
        }
	} catch (...) {
		cout << "Exception during optimization" << endl;
	}
}

void solve_decision_ilp(Graph const &graph, int const r, int const B, std::vector<int> &burning_sequence, GurobiErrorInfo &err_info, bool verbose) {
    // TODO
    return;
}

void find_burning_sequence(Graph const &g, int const r, std::vector<int> &burning_sequence, bool verbose) {
    PRECONDITION(r >= 1);
    PRECONDITION(G.order() >= r);

    int LB = 0;
    int UB = get_upper_bound(g, r, false);
    vector<int> temp_burning_sequence;

	ggbp::ilp::GurobiErrorInfo err_info;

    int iteration = 0;
    while(LB <= UB) {
        int B = (UB + LB)/2;

        if (1) {
            cout << "Iteration #" << iteration << endl;
            cout << "Current lower bound = " << LB << endl;
            cout << "Current upper bound = " << UB << endl;
            cout << "Trying to find a burning sequence of length " << B << "." << endl;
        }

    	err_info.msg = "";

        solve_decision_ilp(g, r, B, temp_burning_sequence, err_info, verbose);

        if(int(temp_burning_sequence.size()) == 0) {
            if (1) {
                cout << "There is no burning sequence of length " << B << "." << endl << endl;
            }
            LB = B+1;
        } else {
            if (1) {
                cout << "A burning sequence of length " << B << " was found." << endl << endl;
            }
            UB = B;
            swap(burning_sequence, temp_burning_sequence);

            if (UB == LB)
                break;
        }

        iteration++;
    }

    cout << "Optimal burning sequence has length " << int(burning_sequence.size()) << "." << endl << endl;
}


}