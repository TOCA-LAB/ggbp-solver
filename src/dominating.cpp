void dominating(Graph &graph, int r, vector<int> &burning_sequence, bool verbose) {
    vector<int> dominating_set;

    try {
        min_dominating_set(graph, dominating_set, verbose);
    } catch (...) {
        throw;
    }

    int n_vertices = graph.size();

    int n_burning = 0;
    vector<int> burning(n_vertices);
    vector<int> num_burning_neighbors(n_vertices);

    forward_list<int> to_burn;
    forward_list<int> next_to_burn;
    vector<int> in_to_burn(n_vertices, 0);

    vector<int>::iterator curr_vertex = dominating_set.begin();
    while (n_burning < n_vertices && curr_vertex != dominating_set.end()) {
        for (auto v : to_burn) {
            burning[v] = 1;
            n_burning++;
        }
        
        for (auto v : to_burn) {
            for (auto u : graph.neighborhood(v)) {
                num_burning_neighbors[u]++;
                
                if (num_burning_neighbors[u] >= r && !burning[u] && !in_to_burn[u]) {
                    next_to_burn.push_front(u);
                    in_to_burn[u] = 1;
                }
            }
        }

        while (curr_vertex != dominating_set.end() && (burning[*curr_vertex] || in_to_burn[*curr_vertex])) {
            curr_vertex++;
        }
    
        if (curr_vertex == dominating_set.end())
            break;
                
        burning_sequence.push_back(*curr_vertex);
        burning[*curr_vertex] = 1;
        n_burning++;
        
        for (auto u : graph.neighborhood(*curr_vertex)) {
            num_burning_neighbors[u]++;
            
            if (num_burning_neighbors[u] >= r && !burning[u] && !in_to_burn[u]) {
                next_to_burn.push_front(u);
                in_to_burn[u] = 1;
            }
        }

        curr_vertex++;

        in_to_burn = vector<int>(n_vertices, 0);
        to_burn = next_to_burn;
        next_to_burn.clear();
    }

    to_burn.clear();
    next_to_burn.clear();

    while (n_burning < n_vertices) {

        for (auto v : to_burn) {
            burning[v] = 1;
            n_burning++;
        }
        
        for (auto v : to_burn) {
            for (auto u : graph.neighborhood(v)) {
                num_burning_neighbors[u]++;
                
                if (num_burning_neighbors[u] >= r && !burning[u] && !in_to_burn[u]) {
                    next_to_burn.push_front(u);
                    in_to_burn[u] = 1;
                }
            }
        }
        
        int max_vertex = 0;
        int max_neighbors = -1;

        for (int u = 0; u < n_vertices; u++) {
            if (!burning[u] && !in_to_burn[u]) {
                int new_spontaneous = 0;
                for (auto v : graph.neighborhood(u)) {
                    if (num_burning_neighbors[v] == r-1 && !burning[v] && !in_to_burn[v]) {
                        new_spontaneous++;
                    }
                }

                if (new_spontaneous > max_neighbors) {
                    max_vertex = u;
                    max_neighbors = new_spontaneous;
                }
            }
        }

        burning_sequence.push_back(max_vertex);
        burning[max_vertex] = 1;
        n_burning++;

        for (auto u : graph.neighborhood(max_vertex)) {
            num_burning_neighbors[u]++;
            
            if (num_burning_neighbors[u] >= r && !burning[u] && !in_to_burn[u]) {
                next_to_burn.push_front(u);
                in_to_burn[u] = 1;
            }
        }

        in_to_burn = vector<int>(n_vertices, 0);
        to_burn = next_to_burn;
        next_to_burn.clear();
    }
}

