#include <algorithm>
#include <limits>
#include <stdexcept>

#include "ggbp/burning_simulation.hpp"

namespace ggbp::simulation {

BurningState::BurningState(const Graph& graph, const int threshold)
    : graph_(&graph),
      threshold_(threshold),
      burned_(graph.order(), false),
      burned_neighbor_counts_(graph.order(), 0),
      in_frontier_(graph.order(), false),
      burned_rounds_(graph.order(), -1) {
    if (threshold < 1) {
        throw std::invalid_argument("burning threshold must be positive");
    }
    if (graph.order() == 0) {
        completion_turns_ = 0;
    }
}

void BurningState::validate_vertex(const Vertex vertex) const {
    if (vertex < 0 ||
        static_cast<Graph::size_type>(vertex) >= graph_->order()) {
        throw std::out_of_range("burning source is outside the graph");
    }
}

void BurningState::burn(const Vertex source) {
    validate_vertex(source);

    if (turn_count_ > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error("burning round does not fit in an int");
    }

    const int current_round = static_cast<int>(turn_count_);
    const bool source_burns = !burned_[source];
    std::size_t burned_this_turn = 0;

    if (source_burns) {
        burned_[source] = true;
        burned_rounds_[source] = current_round;
        ++burned_count_;
        ++burned_this_turn;
    }

    for (const Vertex vertex : frontier_) {
        in_frontier_[vertex] = false;
        if (!burned_[vertex]) {
            burned_[vertex] = true;
            burned_rounds_[vertex] = current_round;
            ++burned_count_;
            ++burned_this_turn;
        }
    }

    next_frontier_.clear();
    const auto propagate_from = [this](const Vertex vertex) {
        for (const Vertex neighbor : graph_->neighbors(vertex)) {
            if (burned_[neighbor]) {
                continue;
            }

            auto& count = burned_neighbor_counts_[neighbor];
            ++count;
            if (count == static_cast<unsigned>(threshold_)) {
                in_frontier_[neighbor] = true;
                next_frontier_.push_back(neighbor);
            }
        }
    };

    if (source_burns) {
        propagate_from(source);
    }

    for (const Vertex vertex : frontier_) {
        if (vertex != source) {
            propagate_from(vertex);
        }
    }

    std::sort(next_frontier_.begin(), next_frontier_.end());
    frontier_.swap(next_frontier_);
    next_frontier_.clear();

    burned_per_turn_.push_back(burned_this_turn);
    ++turn_count_;

    if (complete() && !completion_turns_) {
        completion_turns_ = turn_count_;
    }
}

bool BurningState::complete() const noexcept {
    return burned_count_ == graph_->order();
}

bool BurningState::is_burned(const Vertex vertex) const {
    validate_vertex(vertex);
    return burned_[vertex];
}

bool BurningState::is_in_frontier(const Vertex vertex) const {
    validate_vertex(vertex);
    return in_frontier_[vertex];
}

unsigned BurningState::burned_neighbor_count(const Vertex vertex) const {
    validate_vertex(vertex);
    return burned_neighbor_counts_[vertex];
}

const Graph& BurningState::graph() const noexcept {
    return *graph_;
}

int BurningState::threshold() const noexcept {
    return threshold_;
}

std::size_t BurningState::turn_count() const noexcept {
    return turn_count_;
}

std::size_t BurningState::burned_count() const noexcept {
    return burned_count_;
}

std::optional<std::size_t> BurningState::completion_turns() const noexcept {
    return completion_turns_;
}

std::optional<Vertex> BurningState::first_unburned() const noexcept {
    const auto position = std::find(burned_.begin(), burned_.end(), false);

    if (position == burned_.end()) {
        return std::nullopt;
    }

    return static_cast<Vertex>(std::distance(burned_.begin(), position));
}

std::span<const Vertex> BurningState::frontier() const noexcept {
    return frontier_;
}

const std::vector<bool>& BurningState::burned_vertices() const noexcept {
    return burned_;
}

const std::vector<unsigned>&
BurningState::burned_neighbor_counts() const noexcept {
    return burned_neighbor_counts_;
}

const std::vector<int>& BurningState::burned_rounds() const noexcept {
    return burned_rounds_;
}

const std::vector<std::size_t>& BurningState::burned_per_turn() const noexcept {
    return burned_per_turn_;
}

BurningResult BurningState::result() const {
    return BurningResult{
        .feasible = complete(),
        .simulated_turns = turn_count_,
        .completion_turns = completion_turns_,
        .first_unburned = first_unburned(),
        .burned_round = burned_rounds_,
        .burned_per_turn = burned_per_turn_,
    };
}

BurningResult simulate_burning(
    const Graph& graph,
    const int threshold,
    const std::span<const Vertex> burning_sequence) {
    BurningState state(graph, threshold);

    for (const Vertex source : burning_sequence) {
        state.burn(source);
    }
    
    return state.result();
}

}