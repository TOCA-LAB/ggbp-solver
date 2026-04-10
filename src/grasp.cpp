#include <vector>
#include <iostream>
#include <chrono>
#include <random>

#include "ggbp/ggbp.hpp"


namespace ggbp::grasp {

using namespace std;

int calculateCost(int v, Graph const &graph, int r, const vector<int> &sequence, const vector<int> &burning, const vector<int> &burning_neighbors) {
    if (burning[v])
        return 0;

    int cost = -1, degree = graph.degree(v), n_vertices = graph.order();

    if (degree < r)
        cost -= 10 * n_vertices * graph.degree(v);

    for (auto neighbor : graph.neighbors(v))
        if (!burning[neighbor] && burning_neighbors[neighbor] < r)
            cost -= n_vertices/(graph.degree(neighbor)*(r - burning_neighbors[neighbor]));
       
    cost -= degree;

    return cost;
}

vector<int> construct(Graph const &graph, int r) {
    auto t1 = chrono::high_resolution_clock::now();

    int n_vertices = graph.order(), n_burning = 0;

    vector<int> burning_sequence, cost(n_vertices), burning(n_vertices), burning_neighbors(n_vertices), to_burn;

    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    static default_random_engine e(seed);
    static uniform_real_distribution<float> f(0,1);
    float a = f(e);

    auto t2 = chrono::high_resolution_clock::now();
    std::chrono::_V2::system_clock::time_point t31, t32;

    while (n_burning < n_vertices) {
        t31 = chrono::high_resolution_clock::now();

        to_burn.clear();
        vector<int> in_to_burn(n_vertices, 0);

        cost[0] = calculateCost(0, graph, r, burning_sequence, burning, burning_neighbors); 

        int g_min = cost[0], g_max = cost[0];
        for (int i = 1; i < n_vertices; i++) {
            cost[i] = calculateCost(i, graph, r, burning_sequence, burning, burning_neighbors);

            if (cost[i] < g_min)
                g_min = cost[i];
            if (cost[i] > g_max)
                g_max = cost[i];
        }
        
        int threshold = g_min + a * (g_max - g_min);

        vector<int> RCL;
        RCL.reserve(n_vertices/10);

        for (int i = 0; i < n_vertices; i++)
            if (cost[i] <= threshold)
                RCL.push_back(i);

        RCL.shrink_to_fit();

        uniform_int_distribution<> u(0, RCL.size()-1);
        int current;
        do {
            current = RCL[u(e)];
        } while (burning[current]);
    
        burning_sequence.push_back(current);
        burning[current] = 1;
        n_burning++;

        for (auto neighbor : graph.neighbors(current)) {
            burning_neighbors[neighbor]++;

            if (!burning[neighbor] && !in_to_burn[neighbor] && burning_neighbors[neighbor] >= r) {
                to_burn.push_back(neighbor);
                in_to_burn[neighbor] = 1;
            }
        }

        for (auto v : to_burn) {
            burning[v] = 1;
            n_burning++;
        }

        if (n_burning == n_vertices)
            burning_sequence.push_back(graph.neighbors(current)[0]);

        t32 = chrono::high_resolution_clock::now();
    }
    auto t4 = chrono::high_resolution_clock::now();
    chrono::duration<double, std::milli> pri = t2 - t1;
    chrono::duration<double, std::milli> seg = t32 - t31;
    chrono::duration<double, std::milli> ter = t4 - t1;
    cout << pri.count()/1000 << ", " << seg.count()/1000 << ", " << ter.count()/1000 << endl;
    return burning_sequence;
}

void localSearch(vector<int>& burning_sequence, Graph const &graph, int r) {
    int n_vertices = graph.order(), n_burning = 0;
    auto B = burning_sequence.size();
    vector<int> temp_burning_sequence, burning(n_vertices), burning_neighbors(n_vertices), to_burn;
    temp_burning_sequence.reserve(B);

    // static default_random_engine e;
    
    int i = 0;
    while (i < burning_sequence.size() && n_burning < n_vertices) {
        to_burn.clear();
        to_burn.reserve(B);
        vector<int> in_to_burn(n_vertices, 0);

        int current = burning_sequence[i];
        for (auto v : graph.neighbors(burning_sequence[i])) {
            if (!burning[v]) {
                current = v;
                break;
            }
        }
        if (current = burning_sequence[i]) {
            i++;
            continue;
        }

        temp_burning_sequence.push_back(current);
        burning[current] = 1;
        n_burning++;

        // cout << current << endl;
        for (auto neighbor : graph.neighbors(current)) {
            burning_neighbors[neighbor]++;
            if (!burning[neighbor] && !in_to_burn[neighbor] && burning_neighbors[neighbor] >= r) {
                to_burn.push_back(neighbor);
                in_to_burn[neighbor] = 1;
            }
        } 

        for (auto v : to_burn) {
            burning[v] = 1;
            n_burning++;
        }

        if (n_burning == n_vertices)
            temp_burning_sequence.push_back(graph.neighbors(current)[0]);
        
        i++;
        
    }
    if (i < burning_sequence.size())
        swap(burning_sequence, temp_burning_sequence);
}

void grasp(Graph const &graph, int r, int max_it, vector<int>& burning_sequence, bool verbose) {
    vector<int> temp_burning_sequence;
    
    for (int i = 0; i < max_it; i++) {
        cout << "Iteration " << i+1 << endl;
        temp_burning_sequence = construct(graph, r);
        localSearch(temp_burning_sequence, graph, r);

        if (burning_sequence.size() == 0)
            swap(temp_burning_sequence, burning_sequence);
        else if (temp_burning_sequence.size() < burning_sequence.size())
            swap(temp_burning_sequence, burning_sequence);
    }
}

}