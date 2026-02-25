
void burn_communities_un(Graph &graph, int r, string graph_bin_path, string louvain_path,
                         vector<int> &burning_sequence) {
    graph.get_communities(graph_bin_path, louvain_path);
    
    int n_vertices = graph.size();
    int n_burning = 0;
    vector<int> burning(n_vertices, 0);
    vector<int> num_burning_neighbors(n_vertices);

    forward_list<int> next_to_burn;
    vector<int> in_to_burn(n_vertices, 0);
    forward_list<int> to_burn;

    int num_communities = graph.communities.size();
    int current_community = 0;
    int community_iters = 0;
    vector<bool> burned_community(num_communities, false); 
    
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

        int max_vertex = -1, max_neighbors;
        while (max_vertex == -1 && n_burning < n_vertices) {
            max_vertex = -1;
            max_neighbors = -1;
            
            for (int u : graph.communities[current_community]) {
                if (graph.degree(u) - num_burning_neighbors[u] > max_neighbors && !burning[u] && !in_to_burn[u]) {
                    max_vertex = u;
                    max_neighbors = graph.degree(u) - num_burning_neighbors[u];
                }

            }

            if (max_vertex == -1) {
                int num_changes = 0;
                burned_community[current_community] = true;
                while (burned_community[current_community]) {
                    current_community = (current_community + 1) % num_communities;
                    num_changes++;

                    if (num_changes >= num_communities) {
                        max_vertex = 1;
                        break;
                    }
                }
                community_iters = 0;
            }
        }

        if (max_vertex == -1)
            max_vertex = 1;

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

        community_iters++;
        if (community_iters > r || burned_community[current_community]) {
            int num_changes = 0;
            while (burned_community[current_community]) {
                current_community = (current_community + 1) % num_communities;
                num_changes++;

                if (num_changes >= num_communities)
                    break;
            }
            community_iters = 0;
        }
    }
}
