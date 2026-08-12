#ifndef CONTRACT_HPP
#define CONTRACT_HPP

#include <stdexcept>
#include <sstream>
#include <string>
#include <ranges>
#include <concepts>
#include <algorithm>
#include <vector>

// Enable/disable contracts via macro (define NO_CONTRACTS para desabilitar)
#if !defined(NO_CONTRACTS)
#define CONTRACT_CHECK(expr, kind) \
     do { \
         if (!(expr)) { \
             std::ostringstream oss; \
             oss << #kind << " failed: " << #expr << " at " << __FILE__ \
                 << ":" << __LINE__; \
             throw std::logic_error(oss.str()); \
         } \
     } while(0)
#else
#define CONTRACT_CHECK(expr, kind) (void)0
#endif

#define PRECONDITION(expr)  CONTRACT_CHECK(expr, "Precondition")
#define POSTCONDITION(expr) CONTRACT_CHECK(expr, "Postcondition")
#define INVARIANT(expr)     CONTRACT_CHECK(expr, "Invariant")


template <std::ranges::input_range R>
requires std::is_arithmetic_v<std::ranges::range_value_t<R>>
bool all_in_interval(const R& range,
                     const std::ranges::range_value_t<R>& a,
                     const std::ranges::range_value_t<R>& b)
{
    auto it = std::ranges::find_if(range, [&](const auto& x){
        return x < a || x > b;
    });
    return it == std::ranges::end(range);
}


template <std::ranges::input_range R, typename T>
bool all_elements_in(const R& collection, const std::vector<T>& allowed)
{
    for (const auto& x : collection) {
        if (std::find(allowed.begin(), allowed.end(), x) == allowed.end())
            return false;
    }
    return true;
}

#endif
