#include <cmath>
#include "gurobi_c++.h"
#include "utils/contract.hpp"
#include "brkga-lib/MTRand.h"
#include "brkga-lib/BRKGA.h"
#include "ggbp/graph.hpp"
#include "ggbp/ggbp.hpp"
#include "utils/utils.hpp"

#include "ggbp/ggbp.hpp"

namespace ggbp::ilp::naive {
using namespace std;

int get_upper_bound(Graph const &g, int const r, int verbose) {
    const unsigned n = g.order();
    const unsigned p = 60;
    const unsigned MAX_GENS = 25;
    const double pe = 0.20;
    const double pm = 0.10;
    const double rhoe = 0.70;	
    const unsigned K = 2;
    const unsigned MAXT = 2;

    ggbp::brkga::Decoder decoder(g, r);

    const long unsigned rngSeed = 0;
    MTRand rng(rngSeed);

	BRKGA< ggbp::brkga::Decoder, MTRand > algorithm(n, p, pe, pm, rhoe, decoder, rng, K, MAXT);
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

void add_active_neighborhood_logic(GRBModel &model, Graph const &g, vector<vector<GRBVar>> x, vector<vector<GRBVar>> w, int const r, int const B) {
    auto num_vertices = g.order();
    for (int u = 0; u < num_vertices; u++) {
        if (g.degree(u) < r)
            break;

        for (int t = 2; t < B; t++) {
            GRBLinExpr expr = 0;

            for (auto v : g.neighbors(u)) {
                expr += x[v][t-1] - x[v][t-2];
            }

            model.addConstr(w[u][t] <= expr, "Ngbr_" + to_string(u) + "_" + to_string(t));
        }
    }
}

void solve_ilp(Graph const &g, int const r, vector<int> &burning_sequence, GurobiErrorInfo &err_info, bool verbose) {
    int const n_vertices = g.order();
    int const B = get_upper_bound(g, r, false) + 1;
    cout << B << endl;
    
    bool linear_relax = false;

    char vtype = GRB_BINARY;

    if (linear_relax) {
        vtype = GRB_CONTINUOUS;
    }

    bool opt = true;


    try {        
        GRBEnv env(true);
        env.set(GRB_IntParam_OutputFlag, verbose);
        env.start();
        GRBModel model(env);

        model.set(GRB_DoubleParam_TimeLimit, 7200.0);
        model.set(GRB_DoubleParam_SoftMemLimit, 13.0);
        model.set(GRB_DoubleParam_MemLimit, 13.5);

        // --- Variáveis ---
        // x[v, t]: 1 se o vértice v está queimado no tempo t
        vector<vector<GRBVar>> x(n_vertices, vector<GRBVar>(B));

        // y[v, t]: 1 se o vértice v é inflamado manualmente no tempo t
        vector<vector<GRBVar>> y(n_vertices, vector<GRBVar>(B));
        
        // w[v, t]: 1 se o vértice v é inflamado espontaneamente no tempo t
        vector<vector<GRBVar>> w(n_vertices, vector<GRBVar>(B));

        for (int u = 0; u < n_vertices; u++) {
            x[u][0] = model.addVar(0, 1, 0, vtype, "x_" + to_string(u) + "_0");

            for (int t = 1; t < B; t++) {
                x[u][t] = model.addVar(0, 1, 0, vtype, "x_" + to_string(u) + "_" + to_string(t));
                y[u][t] = model.addVar(0, 1, 0, vtype, "y_" + to_string(u) + "_" + to_string(t));
                w[u][t] = model.addVar(0, 1, 0, vtype, "w_" + to_string(u) + "_" + to_string(t));
            }
        }

        // --- Restrições ---
        // 0. Estado inicial: ninguém queimado
        for (int u = 0; u < n_vertices; u++) {
            model.addConstr(x[u][0] == 0, "init");
        }

        GRBLinExpr expr;
        for (int t = 1; t < B; t++) {
            // 1. Exatamente uma ignição manual por passo
            if (!opt) {
                cout << t << endl;
                expr = 0;
                for (int u = 0; u < n_vertices; u++) {
                    expr += y[u][t];
                }
                model.addConstr(expr == 1, "ignitor_" + to_string(t));
            }

            for (int u = 0; u < n_vertices; u++) {
                // 2. Persistência: uma vez queimado, sempre queimado
                model.addConstr(x[u][t] >= x[u][t-1], "stay_" + to_string(u) + "_" + to_string(t));
            
                //  3. Lógica de Difusão (Threshold r)
                if (g.degree(u) >= r) {
                    // w[v,t] só pode ser 1 se r vizinhos estavam queimados em t-1
                    expr = 0;
                    for (int v : g.neighbors(u))
                        expr += x[v][t-1];
                    model.addConstr(r * w[u][t] <= expr, "diff_" + to_string(u) + "_" + to_string(t));
                } else {
                    model.addConstr(w[u][t] == 0);
                }

                // 4. Transição: x[v,t] é 1 se já estava queimado, foi ignitado ou difundido
                model.addConstr(x[u][t] <= x[u][t-1] + y[u][t] + w[u][t], "trans_" + to_string(u) + "_" + to_string(t));
            }
        }

        // 5. Condição de Finalização: TODOS os vértices devem estar queimados no tempo B
        for (int u = 0; u < n_vertices; u++) {
            model.addConstr(x[u][B-1] == 1, "burned_" + to_string(u));
        }

        add_active_neighborhood_logic(model, g, x, w, r, B);


        // vector<Vertex> sequence = {2, 33, 5, 6, 11};

        // for (int t = 1; t < B; t++) {
        //     model.addConstr(y[sequence[t-1]][t] == 1);
        // }

        // Se escolher versão de otimização
        if (opt) {
            // Variável de otimização
            // done[t]: 1 se todos os vértices estão queimados no turno t
            vector<GRBVar> done(B);

            for (int t = 0; t < B; t++) {
                done[t] = model.addVar(0, 1, 0, vtype, "done_" + to_string(t));
            }

            for (int t = 1; t < B; t++) {
                expr = 0;
                for (int u = 0; u < n_vertices; u++) {
                    expr += y[u][t];
                }
                model.addConstr(expr == 1 - done[t-1], "ignition_limit_done_" + to_string(t));

                model.addConstr(done[t] >= done[t-1], "done_persistence" + to_string(t));

                for (int u = 0; u < n_vertices; u++) {
                    model.addConstr(done[t] <= x[u][t], "check_done_" + to_string(u) + "_" + to_string(t));
                }
            }
            // Minimiza a quantidade de turnos até todos vértices serem queimados
            GRBLinExpr obj = 0;
            for (int t = 1; t < B; t++) {
                obj += 1 - done[t];
            }
            model.setObjective(obj + 1, GRB_MINIMIZE);
        }


        model.optimize();

        int opt_status = model.get(GRB_IntAttr_Status);

        burning_sequence.clear();
        if (opt_status == GRB_OPTIMAL) {
			for(int t = 1; t < B; t++) {
            	for(int u = 0; u < n_vertices; u++) {
                    if (y[u][t].get(GRB_DoubleAttr_X) > 0.5) {
                        burning_sequence.push_back(u);
                    }
                }
            }

            Vertex u;
            cout << ggbp::utils::is_burning_sequence(g, r, burning_sequence, u) << endl;
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

void solve_decision_ilp(Graph const &g, int const r, int const B, std::vector<int> &burning_sequence, GurobiErrorInfo &err_info, bool verbose) {
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