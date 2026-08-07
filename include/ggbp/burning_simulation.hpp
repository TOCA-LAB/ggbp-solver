#ifndef GGBP_BURNING_SIMULATION_HPP
#define GGBP_BURNING_SIMULATION_HPP

#include <cstddef>
#include <optional>
#include <span>
#include <vector>

#include "ggbp/graph.hpp"

namespace ggbp::simulation {

struct BurningResult {
    bool feasible = false;
    std::size_t simulated_turns = 0;
    std::optional<std::size_t> completion_turns;
    std::optional<Vertex> first_unburned;
    std::vector<int> burned_round;
    std::vector<std::size_t> burned_per_turn;
};

class BurningState {
public:
    BurningState(const Graph& graph, int threshold);
    BurningState(Graph&& graph, int threshold) = delete;

    void burn(Vertex source);

    [[nodiscard]] bool complete() const noexcept;
    [[nodiscard]] bool is_burned(Vertex vertex) const;
    [[nodiscard]] bool is_in_frontier(Vertex vertex) const;
    [[nodiscard]] unsigned burned_neighbor_count(Vertex vertex) const;

    [[nodiscard]] const Graph& graph() const noexcept;
    [[nodiscard]] int threshold() const noexcept;
    [[nodiscard]] std::size_t turn_count() const noexcept;
    [[nodiscard]] std::size_t burned_count() const noexcept;
    [[nodiscard]] std::optional<std::size_t> completion_turns() const noexcept;
    [[nodiscard]] std::optional<Vertex> first_unburned() const noexcept;

    [[nodiscard]] std::span<const Vertex> frontier() const noexcept;
    [[nodiscard]] const std::vector<bool>& burned_vertices() const noexcept;
    [[nodiscard]] const std::vector<unsigned>& burned_neighbor_counts() const noexcept;
    [[nodiscard]] const std::vector<int>& burned_rounds() const noexcept;
    [[nodiscard]] const std::vector<std::size_t>& burned_per_turn() const noexcept;

    [[nodiscard]] BurningResult result() const;

private:
    const Graph* graph_;
    int threshold_;
    std::vector<bool> burned_;
    std::vector<unsigned> burned_neighbor_counts_;
    std::vector<Vertex> frontier_;
    std::vector<bool> in_frontier_;
    std::vector<int> burned_rounds_;
    std::vector<std::size_t> burned_per_turn_;
    std::size_t turn_count_ = 0;
    std::size_t burned_count_ = 0;
    std::optional<std::size_t> completion_turns_;

    void validate_vertex(Vertex vertex) const;
    void rebuild_frontier();
};

[[nodiscard]] BurningResult simulate_burning(
    const Graph& graph,
    int threshold,
    std::span<const Vertex> burning_sequence);

} // namespace ggbp::simulation

#endif
