#include <cmath>
#include "gurobi_c++.h"
#include "brkga-lib/MTRand.h"
#include "brkga-lib/BRKGA.h"
#include "ggbp/graph.hpp"
#include "ggbp/ggbp.hpp"
#include "ggbp/decoder.hpp"

namespace ggbp::ilp {
using namespace std;

int get_upper_bound(Graph &graph, int const r, int verbose) {
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

void solve_ilp(Graph &graph, int const r, vector<int> &burning_sequence, GurobiErrorInfo &err_info, int verbose) {
    int const n_vertices = graph.order();
    int const K = get_upper_bound(graph, r, verbose);

    try {        
        GRBEnv env(true);
        env.set(GRB_IntParam_OutputFlag, verbose);
        env.start();
        GRBModel model(env);

        model.set(GRB_DoubleParam_TimeLimit, 7200.0);
        model.set(GRB_DoubleParam_SoftMemLimit, 13.0);
        model.set(GRB_DoubleParam_MemLimit, 13.5);

        // Defining the variables
        vector<vector<GRBVar>> vars_s(n_vertices, vector<GRBVar>(K));
        vector<vector<GRBVar>> vars_q(n_vertices, vector<GRBVar>(K));
        vector<GRBVar> vars_z(K);

        for (int i = 0; i < n_vertices; i++) {
            for (int j = 0; j < K; j++) {
                vars_s[i][j] = model.addVar(0, 1, 0, GRB_BINARY, "s_" + to_string(i) + "_" + to_string(j));
                vars_q[i][j] = model.addVar(0, 1, 0, GRB_BINARY, "q_" + to_string(i) + "_" + to_string(j));
            }
        }

        for (int i = 0; i < K; i++)
            vars_z[i] = model.addVar(0, 1, 0, GRB_BINARY, "z_" + to_string(i)); 
        
        // Adding the constraints

        // One ignitor vertex for each used turn
        for (int i = 0; i < K; i++) {
            GRBLinExpr expr = 0;
            for (int j = 0; j < n_vertices; j++)
                expr += vars_s[j][i];
            model.addConstr(expr == vars_z[i], "MIT_" + to_string(i));
        }

        // A vertex can only be ignitor once
        for (int i = 0; i < n_vertices; i++) {
            GRBLinExpr expr = 0;
            for (int j = 0; j < K; j++)
                expr += vars_s[i][j];
            model.addConstr(expr <= 1, "I_" + to_string(i));
        }

        // An ignitor also burns
        for (int i = 0; i < n_vertices; i++)
            for (int j = 0; j < K; j++)
                model.addConstr(vars_s[i][j] <= vars_q[i][j], "QI_" + to_string(i) + "_" + to_string(j));

        // A vertex burns only if there are r burning neighbours or if it is a source 
        for (int i = 0; i < n_vertices; i++)
            for (int j = 0; j < K-1; j++) {
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
                for (int i = 0; i < K; i++)
                    expr += vars_s[u][i];
                model.addConstr(expr == 1, "S_" + to_string(u));
            }
        }
        
        // A burned vertex continues burning throughout the entire process
        for (int i = 0; i < n_vertices; i++)
            for (int j = 0; j < K-1; j++)
                model.addConstr(vars_q[i][j] <= vars_q[i][j+1], "Q_" + to_string(i) + "_" + to_string(j));

        // On the first turn, all vertices are not burned, except the first source
        for (int i = 0; i < n_vertices; i++)
            model.addConstr(vars_q[i][0] == vars_s[i][0], "T0_" + to_string(i));

        // The first turn always happens
        model.addConstr(vars_z[0] == 1);

        // If there's a vertex not burning, a new turn is needed
        for (int i = 0; i < n_vertices; i++)
            for (int j = 0; j < K-1; j++)
                model.addConstr((1 - vars_q[i][j]) <= vars_z[j+1], "NT_" + to_string(i) + "_" + to_string(j));

        // All vertices must be burned at the end
        GRBLinExpr expr = 0;
        for (int i = 0; i < n_vertices; i++)
            expr += vars_q[i][K-1];
        model.addConstr(expr == n_vertices, "T_K");

        // Minimize the number of turns
        GRBLinExpr obj = 0;
        for (int i = 0; i < K; i++)
            obj += vars_z[i];
        model.setObjective(obj, GRB_MINIMIZE);

        for (int i = 0; i < n_vertices; i++)
            for (int j = 0; j < K-1; j++) {
                // vars_s[i][j].set(GRB_IntAttr_BranchPriority, K-j);
                // vars_z[j].set(GRB_IntAttr_BranchPriority, K-j);
            } 

        model.optimize();

        int opt_status = model.get(GRB_IntAttr_Status);

        if (opt_status == GRB_OPTIMAL) {
			for(int i = 0; i < K; i++) {
            	for(int j = 0; j < n_vertices; j++) {
					if(vars_s[j][i].get(GRB_DoubleAttr_X) > 0.5)
						burning_sequence.push_back(j);
                }
            }
            // printVector(burning_sequence, "Burning sequence");
        } else if (opt_status == GRB_TIME_LIMIT) {
            if (model.get(GRB_DoubleAttr_ObjVal) < K)
                err_info.ObjVal = model.get(GRB_DoubleAttr_ObjVal);
            else
                err_info.ObjVal = K;
            
            err_info.ObjBound = model.get(GRB_DoubleAttr_ObjBound);
            err_info.msg = "//TIME LIMIT";
        } else if (opt_status == GRB_MEM_LIMIT) {
            if (model.get(GRB_DoubleAttr_ObjVal) < K)
                err_info.ObjVal = model.get(GRB_DoubleAttr_ObjVal);
            else
                err_info.ObjVal = K;

            err_info.ObjBound = model.get(GRB_DoubleAttr_ObjBound);
            err_info.msg = "//OUT OF MEMORY (soft)";
        }
    } catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;

        if (e.getErrorCode() == GRB_ERROR_OUT_OF_MEMORY) {
            err_info.ObjVal = K;
            err_info.ObjBound = 0;
            err_info.msg = "//OUT OF MEMORY";
        }
	} catch (...) {
		cout << "Exception during optimization" << endl;
	}
}
}