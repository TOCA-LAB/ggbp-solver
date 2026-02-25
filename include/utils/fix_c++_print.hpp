#pragma once

#include <iterator>
#include <map>
#include <ostream>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <vector>

// -----------------------------------------------------------------------------
//
//                            << for vectors 
//
// -----------------------------------------------------------------------------
template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
    os << "[";
    for (size_t i = 0; i < v.size(); ++i) {
        os << v[i];
        if (i + 1 < v.size()) os << ", ";
    }
    os << "]";
    return os;
}

// -----------------------------------------------------------------------------
//
//                   << for iterables that are not vectors
//
// -----------------------------------------------------------------------------

// Detect if T is iterable and NOT a string
template <typename T, typename = void>
struct is_iterable : std::false_type {};

template <typename T>
struct is_iterable<T, std::void_t<
    decltype(std::begin(std::declval<T>())),
    decltype(std::end(std::declval<T>()))
>> : std::true_type {};

template <>
struct is_iterable<std::string> : std::false_type {};     // avoid breaking operator<< for std::string

template <>
struct is_iterable<std::string_view> : std::false_type {};

template <typename T>
constexpr bool is_iterable_v = is_iterable<T>::value;

template <typename C>
std::enable_if_t<is_iterable_v<C>, std::ostream&>
operator<<(std::ostream& os, const C& container) {
    os << "{";
    bool first = true;
    for (const auto& elem : container) {
        if (!first) os << ", ";
        first = false;
        os << elem;
    }
    os << "}";
    return os;
}

// -----------------------------------------------------------------------------
//
//                            << for pairs
//
// -----------------------------------------------------------------------------
template <typename A, typename B>
std::ostream& operator<<(std::ostream& os, const std::pair<A, B>& p) {
    os << "(" << p.first << ", " << p.second << ")";
    return os;
}

// -----------------------------------------------------------------------------
//
//                            << for map
//
// -----------------------------------------------------------------------------
template <typename K, typename V, typename... Args>
std::ostream& operator<<(std::ostream& os,
                         const std::map<K, V, Args...>& m) {
    os << "{";
    bool first = true;
    for (const auto& [key, value] : m) {
        if (!first) os << ", ";
        first = false;
        os << key << ": " << value;
    }
    os << "}";
    return os;
}

// -----------------------------------------------------------------------------
//
//                       << for unordered_map
//
// -----------------------------------------------------------------------------
template <typename K, typename V, typename... Args>
std::ostream& operator<<(std::ostream& os,
                         const std::unordered_map<K, V, Args...>& m) {
    os << "{";
    bool first = true;
    for (const auto& [key, value] : m) {
        if (!first) os << ", ";
        first = false;
        os << key << ": " << value;
    }
    os << "}";
    return os;
}

