#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <unordered_set>
#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <ostream>
#include <ranges>

using Vertex = int;

class Graph {
public:
    //using Vertex = int;
    using AdjList = std::vector<std::vector<Vertex>>;
    using size_type = size_t; 

    explicit Graph(size_type nvertices) 
        : nvertices_(nvertices),
          nedges_(0),
          adj_(nvertices),
          distance_matrix(nvertices)
    {}

    Graph(const std::string& graph6_str);

    static Graph from_file(const std::string& fpath);
    static Graph from_mtx_file(const std::string& path);

    Graph() = default;
    Graph(const Graph&) = default;
    Graph(Graph&&) noexcept = default;
    Graph& operator=(const Graph&) = default;
    Graph& operator=(Graph&&) noexcept = default; 
    ~Graph() = default;

    size_type order() const noexcept { return nvertices_; } 
    size_type size()  const noexcept { return nedges_; }

    // preciso melhorar a implementacao desse método
    std::vector<Vertex> vertices(void) const;

    void add_edge(Vertex u, Vertex v);

    bool has_edge(Vertex u, Vertex v) const {
        check_vertex(u);
        check_vertex(v);

        return std::find(adj_[u].begin(), adj_[u].end(), v) != adj_[u].end();
    }

    size_type degree(Vertex v) const {
        check_vertex(v);

        return adj_[v].size();
    }

    const std::vector<Vertex>& neighbors(Vertex v) const {
        check_vertex(v);

        return adj_[v];
    }

    bool remove_edge(Vertex u, Vertex v);

    void calculate_distances(int origin);

    int get_distance(int v1, int v2);

    void clear() noexcept {
        for (auto& s : adj_) s.clear();
        nedges_ = 0;
    }
    
    std::vector<std::vector<int>> distance_matrix;


    friend std::ostream& operator<<(std::ostream& os, const Graph& G);

private:
    size_type nvertices_;
    size_type nedges_;
    AdjList adj_;

    void check_vertex(Vertex v) const {
        if (v > nvertices_)
            throw std::out_of_range("Vertex index out of range");
    }
};


#endif // GRAPH_H

