#ifndef UTIL_H
#define UTIL_H

#include <cstddef>
#include <iostream>
#include <iterator>
#include <utility>
#include <vector>
#include <functional>

#include "ggbp/graph.hpp"


#if NDEBUG
    #define DEBUG(msg)  (void)0
#else
    #define DEBUG(msg) do { std::cerr << "[DEBUG] " << __FILE__ << ":" << __LINE__ << ": "  << msg << std::endl; } while(0)
#endif


namespace ggbp::utils {

// -----------------------------------------------------------------------------
//                            enumerate(Iterable)
// -----------------------------------------------------------------------------

template <typename Iterable>
class EnumerateWrapper {
    Iterable& iterable;

public:
    EnumerateWrapper(Iterable& it) : iterable(it) {}

    struct iterator {
        std::size_t index;
        decltype(std::begin(iterable)) iter;

        bool operator!=(const iterator& other) const { return iter != other.iter; }
        void operator++() { ++index; ++iter; }
        auto operator*() const { return std::pair{index, *iter}; }
    };

    auto begin() { return iterator{0, std::begin(iterable)}; }
    auto end()   { return iterator{0, std::end(iterable)}; }
};

template <typename Iterable>
auto enumerate(Iterable& iterable) {
    return EnumerateWrapper(iterable);
}

bool is_burning_sequence(Graph const &g, int r, std::vector<int> &burning_sequence, Vertex &unburned);

template <class T>
void bucket_sort(std::vector<T> &v, int n, std::function<bool (T, T)> cmp, std::function<double (T)> val) {
    std::vector<std::vector<T>> B(n);

    for (int i = 0; i < n; i++) {
        B[i].reserve((double(n) / v.size()));
    }

    for (auto e : v) {
        B[floor(val(e) * n)].push_back(e);
    }

    unsigned count = 0;
    for (int i = 0; i < n; i++) {
        sort(B[i].begin(), B[i].end(), cmp);

        for (auto e : B[i]) {
            v[count] = e;
            count++;
        }
    }
}

}

#endif
