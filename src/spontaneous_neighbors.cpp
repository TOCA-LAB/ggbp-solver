#include <chrono>
#include <boost/heap/d_ary_heap.hpp>

#include <vector>
#include <algorithm>
#include <ranges>
#include <map>
#include "ggbp/graph.hpp"
#include "utils/combinations.hpp"
#include "utils/contract.hpp"
#include "utils/utils.hpp"
#include "utils/fix_c++_print.hpp"
#include "utils/logger.hpp"

namespace ggbp::heuristics {

using std::vector;
using ggbp::utils::Combinations;
using std::unordered_set;
using std::map;
using ggbp::utils::enumerate;

constexpr int INF_INT = std::numeric_limits<int>::max();
constexpr int NO_ROUND = -1;

// supondo que Vertex já esteja definido
struct Item {
    int priority;
    Vertex vertex;
};

struct MaxPriority {
    bool operator()(Item const& a, Item const& b) const {
        return a.priority < b.priority;
    }
};

struct MinPriority {
    bool operator()(Item const& a, Item const& b) const {
        return a.priority > b.priority;
    }
};

using MaxPQ = boost::heap::d_ary_heap<Item,
                                      boost::heap::arity<4>,          
                                      boost::heap::mutable_<true>,   
                                      boost::heap::compare<MaxPriority>>;

using MinPQ = boost::heap::d_ary_heap<Item,
                                      boost::heap::arity<4>,         
                                      boost::heap::mutable_<true>,    
                                      boost::heap::compare<MinPriority>>;


int 
value_starting_ignitors(const Graph& G, int r, const vector<Vertex>& ignitors) {
    PRECONDITION(r >= 1);
    PRECONDITION(ignitors.size() >= 1);
    PRECONDITION(all_in_interval(ignitors, 0, G.order() - 1));
    
    Vertex u = ignitors[0];
    vector<Vertex> common_neighbors = G.neighbors(u);

    for (size_t i = 1; i < ignitors.size(); i++) {
      u = ignitors[i];
      const vector<Vertex>& nei = G.neighbors(u);
      for (auto it = common_neighbors.begin(); it != common_neighbors.end(); ) 
        if (std::find(nei.begin(), nei.end(), *it) != nei.end())
          it = common_neighbors.erase(it);
        else
          ++it;
    }

    int count = common_neighbors.size();
    POSTCONDITION(count >= 0 && count <= G.order());
    return count;
}

vector<Vertex> 
best_starting_ignitors(const Graph& G, int r) {
    PRECONDITION(r >= 1);
    PRECONDITION(G.order() >= r);

    vector<Vertex> vertices = G.vertices();
    Combinations candidate_sets(vertices, r);

    int best_score = -1;
    vector<Vertex> best_set;
    best_set.reserve(r);

    for (auto candidate : candidate_sets)  {
        int score = value_starting_ignitors(G, r, candidate);
        if (score > best_score) {
            best_score = score;
            best_set = candidate;
        } 
    }

    logger.info(2) << "best_starting_ignitors: " << best_set;

    POSTCONDITION(static_cast<int>(best_set.size()) == r);
    return best_set;
}

Vertex 
choose_next_ignitor(MaxPQ& candidates, const vector<Vertex>& burning_order) {
    PRECONDITION(!candidates.empty());
    PRECONDITION(all_in_interval(burning_order, NO_ROUND, burning_order.size() - 1));

    Item it;
    do {
        it = candidates.top();
        candidates.pop();
        if (burning_order[it.vertex] == NO_ROUND)
            return it.vertex;
    } while (!candidates.empty());

    POSTCONDITION(it.vertex >= 0 && it.vertex < static_cast<int>(burning_order.size()));
    return it.vertex;
}

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
        if (burned_round[u] == -1)
            return false;
    return true;
}

void 
update_at_risk_neighbor_counts(const Graph& G, 
                               const int r, 
                               const vector<int> &burning_round,
                               vector<int> &at_risk_neighbor_count, 
                               Vertex u,
                               int delta) noexcept {
  PRECONDITION(r >= 1);
  PRECONDITION(G.order() == static_cast<int>(burning_round.size()));
  PRECONDITION(G.order() == static_cast<int>(at_risk_neighbor_count.size()));
  PRECONDITION(u >= 0 && u < G.order());
  PRECONDITION(delta == 1 || delta == -1);

  for (Vertex x : G.neighbors(u))
    if (burning_round[x] == NO_ROUND)
      at_risk_neighbor_count[x] += delta;
}


std::pair<MinPQ, std::vector<MinPQ::handle_type>>
setup_burning_queue(const Graph& G,
                   const std::vector<Vertex>& S,
                   int r) {
    PRECONDITION(r >= 1);
    PRECONDITION(S.size() == r);

    const std::size_t n = G.order();

    vector<bool> in_S(n, false);
    for (Vertex v : S) in_S[v] = true;

    MinPQ burning_queue;
    std::vector<MinPQ::handle_type> handles(n);

    for (Vertex w : G.vertices())
        handles[w] = (!in_S[w]) ? burning_queue.push({INF_INT, w})
                                : burning_queue.push({r, w});

    return {std::move(burning_queue), std::move(handles)};
}


std::pair<MaxPQ, std::vector<MaxPQ::handle_type>>
setup_candidate_ignitors_queues(const Graph& G,
                               const std::vector<Vertex>& S,
                               int r,
                               Vertex t) { 
    PRECONDITION(r >= 1);
    PRECONDITION(S.size() == r);
    const std::size_t n = G.order();

    vector<bool> in_S(n, false);
    for (Vertex v : S) in_S[v] = true;

    MaxPQ candidate_ignitors;
    vector<MaxPQ::handle_type> candidate_ignitors_handles(n);

    for (Vertex w : G.vertices()) {
        if (!in_S[w])
          candidate_ignitors_handles[w] = candidate_ignitors.push({0, w});
        else if (w == t)
          candidate_ignitors_handles[w] = candidate_ignitors.push({r, w});
    }

    return {std::move(candidate_ignitors), std::move(candidate_ignitors_handles)};
}

std::pair<std::vector<Vertex>, std::vector<int>>
initialize_burning_state(const Graph& G,
                         const vector<Vertex>& S)
{
    PRECONDITION(!S.empty());
    PRECONDITION(S.size() <= G.order());
    PRECONDITION(all_elements_in(S, G.vertices()));

    const std::size_t n = G.order();

    vector<Vertex> burning_sequence;
    burning_sequence.reserve(n);

    vector<int> burning_round(n, NO_ROUND);



    for (size_t i = 0; i < S.size() - 1; ++i) {
        Vertex u = S[i];
        burning_sequence.push_back(u);
        burning_round[u] = i;
    }

    POSTCONDITION(burning_sequence.size() == S.size() - 1);
    return {std::move(burning_sequence), std::move(burning_round)};
}

vector<Vertex> spontaneous_neighbors(const Graph &G, const int r) {
  PRECONDITION(r >= 1);

  logger.info(1) << "starting greedy...";

  vector<Vertex> S = best_starting_ignitors(G, r);

  Vertex t = S.back();
  logger.info(2) << "t: " << t << std::endl;

  int round = r - 1;

  auto [burning_sequence, burning_round] = initialize_burning_state(G, S);
  auto [burning_queue, burning_queue_handle] = setup_burning_queue(G, S, r);
  auto [candidate_ignitors, candidate_ignitors_handles] = setup_candidate_ignitors_queues(G, S, r, t);

  vector<int> burned_neighbors_count(G.order(), 0);
  vector<int> at_risk_neighbor_count(G.order(), 0);
  while (!burning_queue.empty()) {
    auto [v_round, v] = burning_queue.top();

    logger.info(2) << "[New iter] round: " << round << " v: " << v
                   << " v_round: " << v_round;

    if (v_round > round) {
      Vertex u = choose_next_ignitor(candidate_ignitors, burning_round);
      logger.info(2) << "ignitor: " << u;

      burning_sequence.push_back(u);
      burning_round[u] = round + 1;
      burning_queue.update(burning_queue_handle[u], {round + 1, u});
      ++round;

      if (burned_neighbors_count[u] == r - 1)
        update_at_risk_neighbor_counts(G, r, burning_round, at_risk_neighbor_count, u, -1);
    } else {
      auto [p, u] = burning_queue.top();
      burning_queue.pop();
      for (auto w : G.neighbors(u)) {
        if (burning_round[w] == NO_ROUND) {
          burned_neighbors_count[w] += 1;
          if (burned_neighbors_count[w] == r) {
            logger.info(3) << "Spontaneous: " << w;
            burning_queue.update(burning_queue_handle[w], {round + 1, w});
            burning_round[w] = round + 1;
            update_at_risk_neighbor_counts(G, r, burning_round, at_risk_neighbor_count, w, -1);
          } else if (burned_neighbors_count[w] == r - 1)
                update_at_risk_neighbor_counts(G, r, burning_round, at_risk_neighbor_count, w, 1);
        }
      }
    }
  }

  POSTCONDITION(is_burning_sequence(G, r, burning_sequence));
  return burning_sequence;
}
}

