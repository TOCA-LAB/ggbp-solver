#include <boost/dynamic_bitset.hpp>
#include <iostream>
#include <fstream>
#include <numeric>
#include <queue>
#include <vector>

#include "ggbp/graph.hpp"


void Graph::add_edge(Vertex u, Vertex v) {
    check_vertex(u);
    check_vertex(v);

    // unordered_set prevents duplicates automatically
    adj_[u].push_back(v);
    adj_[v].push_back(u);

    ++nedges_;
}

std::vector<Vertex> Graph::vertices(void) const {
  std::vector<Vertex> vert(order(), 0);
  std::iota(vert.begin(), vert.end(), 0);
  return vert;
}

bool Graph::remove_edge(Vertex u, Vertex v) {
    check_vertex(u);
    check_vertex(v);

    auto uv = std::find(adj_[u].begin(), adj_[u].end(), v);
    auto vu = std::find(adj_[v].begin(), adj_[v].end(), u);

    if (uv != adj_[u].end() && vu != adj_[v].end()) {
        adj_[u].erase(uv);
        adj_[v].erase(vu);

        --nedges_;
        
        return true;
    }

    return false;
}

inline std::ostream& operator<<(std::ostream& os, const Graph& G) {
    os << G.order() << " " << G.size() << "\n";
    for (Vertex u = 0; u < G.order(); ++u) {
        for (auto v : G.adj_[u]) {
            if (u < v)
                os << u << " " << v << "\n";
        }
    }
    return os;
}

Graph Graph::from_file(const std::string& path) {

    std::ifstream in(path);
    if (!in) throw std::runtime_error("unable to open file: " + path);

    size_type n;
    size_type m;
    if (!(in >> n >> m))
        throw std::runtime_error("invalid header in graph file");

    Graph G(n);

    Vertex u, v;
    for (Vertex i = 0; i < m; ++i) {
        in >> u >> v;
        G.add_edge(u, v);
    }

    return G;
}

Graph Graph::from_mtx_file(const std::string& path) {

    std::ifstream in(path);
    if (!in) throw std::runtime_error("unable to open file: " + path);

    std::string line;

    // skip header's comments
    std::getline(in, line);
    while (line[0] == '%')
        std::getline(in, line);

    std::istringstream istream(line);
    size_type nvertices, _, nedges;
    istream >> nvertices >> _ >> nedges;

    Graph G(nvertices);

    Vertex u, v;
    while (in >> u >> v)
        G.add_edge(u - 1, v - 1);

    return G;
}

boost::dynamic_bitset<> string6_to_bitset(const std::string& s) {

    boost::dynamic_bitset<> bs(6 * s.size());

    int ne = s.size();
    for (int i = 0; i < ne; i++) {
        unsigned char c = s[ne - i - 1] - 63;
        for (int j = 0; j < 6; j++) {
            bs[6 * i + j] = c & 1;
            c = c >> 1;
        }
    }

    return bs;
}

void graph6_parse(const std::string& s, std::string& order, std::string& adjstr) {

    if (s[0] != '~') {
        order = s.substr(0, 1);
        adjstr = s.substr(1, s.size() - 1);
    } else if (s[1] == '~') {
        order = s.substr(2, 6);
        adjstr = s.substr(8, s.size() - 8);
    } else {
        order = s.substr(1, 3);
        adjstr = s.substr(4, s.size() - 4);
    }
}

Graph::Graph(const std::string& graph6_str) {

    std::string order_str;
    std::string adj_str;

    graph6_parse(graph6_str, order_str, adj_str);

    long order = string6_to_bitset(order_str).to_ulong();
    nvertices_ = order;
    nedges_ = 0;

    adj_.resize(order);
    distance_matrix.resize(order);
    auto adj_bs = string6_to_bitset(adj_str);

    int idx = adj_bs.size() - 1;
    for (int j = 1; j < order; j++) {
        for (int i = 0; i < j; i++) {
            if (adj_bs[idx] == 1) {
                add_edge(i, j);
            }
            idx -= 1;
        }
    }

}

void Graph::calculate_distances(int origin) {
    size_type n_vertices = order();
    distance_matrix[origin].resize(n_vertices);
    for (int i = 0; i < n_vertices; i++)
        distance_matrix[origin][i] = n_vertices;

    std::queue<int> bfs_queue;
    std::vector<int> visited(n_vertices, 0);
    bfs_queue.push(origin);
    visited[origin] = 1;
    int u, curr_size, distance = 0;
    
    while (!bfs_queue.empty()) {
        curr_size = bfs_queue.size();
        while (curr_size > 0) {
            u = bfs_queue.front();
            bfs_queue.pop();
            
            distance_matrix[origin][u] = distance;

            for (int v : neighbors(u)) {
                if (!visited[v]) {
                    visited[v] = 1;
                    bfs_queue.push(v);
                }
            }
            curr_size--;
        }
        distance++;
    }
}

int Graph::get_distance(int v1, int v2) {
    if (distance_matrix[v1].size() != 0)
        return distance_matrix[v1][v2];
    
    if (distance_matrix[v2].size() != 0)
        return distance_matrix[v2][v1];
    
    calculate_distances(v1);
    return distance_matrix[v1][v2];        
}