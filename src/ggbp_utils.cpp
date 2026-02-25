
bool is_burning_sequence(const Graph& G, int r, vector<Vertex> burning_sequence) {
    

    vector<int> burned_round(G.order(), -1);
    MinPQ pq;
    map<Vertex, MinPQ::handle_type> idx_pq;

    for (auto [i, u] : enumerate(burning_sequence)) {
        idx_pq[u] = pq.push({static_cast<int>(i), u});
    }

    vector<int> nburned_neighbors(G.order(), 0);
    while (!pq.empty()) {
        auto [round, u] = pq.top();
        pq.pop();
        if (burned_round[u] != -1) 
            continue;
        burned_round[u] = round;
        for (Vertex w : G.neighbors(u)) {
            if (burned_round[w] == -1) {
                nburned_neighbors[w] += 1;
                if (nburned_neighbors[w] == r) {
                    pq.push({round + 1, w});
                }
            }
        }
    }

    for (Vertex u : G.vertices())
        if (nburned_neighbors[u] ==  -1)
            return false;
    return true;
}
