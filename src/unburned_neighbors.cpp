void unburned_neighbors(Graph &graph, int r, vector<int> &burning_sequence) {
    int n_vertices = graph.size();
    
    int n_burning = 0;
    vector<int> burning(n_vertices, 0);

    vector<int> num_burning_neighbors(n_vertices);

    burning_sequence.clear();

    forward_list<int> to_burn;

    while (n_burning < n_vertices) {
        forward_list<int> next_to_burn;
        vector<int> in_to_burn(n_vertices, 0);

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
        int max_neighbors = 0;


        for (int u = 0; u < n_vertices; u++) {
            if (graph.degree(u) - num_burning_neighbors[u] > max_neighbors && !burning[u] && !in_to_burn[u]) {
                max_vertex = u;
                max_neighbors = graph.degree(u) - num_burning_neighbors[u];
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
        
        to_burn = next_to_burn;
        next_to_burn.clear();
    }
}

