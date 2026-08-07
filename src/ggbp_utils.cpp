#include <limits>
#include <stdexcept>
#include <vector>

#include "ggbp/burning_simulation.hpp"
#include "utils/utils.hpp"

namespace ggbp::utils {

bool is_burning_sequence(
    const Graph& graph,
    const int threshold,
    const std::vector<int>& burning_sequence,
    Vertex& first_unburned) {
    const auto result = simulation::simulate_burning(
        graph, threshold, burning_sequence);
    first_unburned = result.first_unburned.value_or(-1);
    return result.feasible;
}

int burning_turns(
    const Graph& graph,
    const int threshold,
    const std::vector<int>& burning_sequence) {
    const auto result = simulation::simulate_burning(
        graph, threshold, burning_sequence);
    if (!result.completion_turns) {
        return -1;
    }
    if (*result.completion_turns >
        static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error("burning turn count does not fit in an int");
    }
    return static_cast<int>(*result.completion_turns);
}

bool burning_history(
    const Graph& graph,
    const int threshold,
    const std::vector<int>& burning_sequence,
    std::vector<int>& history) {
    const auto result = simulation::simulate_burning(
        graph, threshold, burning_sequence);

    history.clear();
    history.reserve(result.burned_per_turn.size());
    for (const std::size_t count : result.burned_per_turn) {
        if (count > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            throw std::overflow_error(
                "burned vertex count does not fit in an int");
        }
        history.push_back(static_cast<int>(count));
    }
    return result.feasible;
}

} // namespace ggbp::utils
