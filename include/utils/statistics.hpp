#ifndef STAT_HPP
#define STAT_HPP

#include <chrono>
#include <deque>
#include <iostream>
#include <string>

#ifndef HOTSTATS
#define HOTSTATS 0
#endif


namespace ggbp::stats {

#if HOTSTATS 

class ElapsedTimer {
public:
    using clock = std::chrono::steady_clock;

    ElapsedTimer(std::string label) noexcept
        : label_(std::move(label))  {}

    void start() noexcept {
        start_ = clock::now();
        running_ = true;
    }

    void stop() noexcept {
        if (running_) {
            auto end = clock::now();
            elapsed_ms_ +=
                std::chrono::duration_cast<std::chrono::milliseconds>(end - start_)
                    .count();
            running_ = false;
        }
    }

    void print(std::ostream& os = std::cerr) const {
        os << "[STATS]\t" << label_ << ": " << elapsed_ms_ << " ms\n";
    }

private:
    std::string label_;
    clock::time_point start_;
    long long elapsed_ms_ = 0;
    bool running_ = false;
};

template <typename T>
class Counter {
public:
    Counter(std::string label)
        : label_(std::move(label)) {}

    void add(T v) { value_ += v; }

    void print(std::ostream& os = std::cerr) const {
        os << "[STATS]\t" << label_ << ": " << value_ << "\n";
    }

private:
    std::string label_;
    T value_ = 0;
};


template <typename T>
class Value {
public:
    Value(std::string label)
        : label_(std::move(label)) {}

    void set(T v) { value_ = v; }

    void print(std::ostream& os = std::cerr) const {
        os << "[STATS]\t" << label_ << ": " << value_ << "\n";
    }

private:
    std::string label_;
    T value_ = T{};
};

/*
  MetricTrace<T>
  - armazena um histórico (lista) de valores do tipo T
  - max_entries == 0 -> sem limite
*/
template <typename T>
class Trace {
public:
    Trace(std::string label,
                std::size_t max_entries = 0)
        : label_(std::move(label)),
          max_entries_(max_entries) {}

    void push(const T& value) {
        if (max_entries_ > 0 && data_.size() == max_entries_) {
            data_.pop_front();
        }
        data_.push_back(value);
    }

    // C++17: rvalue overload
    void push(T&& value) {
        if (max_entries_ > 0 && data_.size() == max_entries_) {
            data_.pop_front();
        }
        data_.push_back(value);
    }

    void clear() { data_.clear(); }

    void print(std::ostream& os = std::cerr) const {
        os << "[STATS]\t" << label_ << ": [";
        bool first = true;
        for (const auto& v : data_) {
            if (!first) os << ", ";
            os << v;
            first = false;
        }
        os << "]\n";
    }

private:
    std::string label_;
    std::deque<T> data_;
    std::size_t max_entries_ = 0;
};


#else


class ElapsedTimer {
public:
    explicit ElapsedTimer(std::string = {}) noexcept {}
    inline void start() noexcept {}
    inline void stop() noexcept {}
    inline void print(std::ostream& = std::cerr) const noexcept {}
    inline const std::string& label() const noexcept { static std::string s; return s; }
};


template <typename T>
class Counter {
public:
    explicit Counter(std::string = {}) noexcept {}
    inline void add(T) noexcept {}
    inline void print(std::ostream& = std::cerr) const noexcept {}
};


template <typename T>
class Value {
public:
    explicit Value(std::string = {}) noexcept {}
    inline void set(T) noexcept {}
    inline void print(std::ostream& = std::cerr) const noexcept {}
};


template <typename T>
class Trace {
public:
    Trace(std::string = {}, std::size_t = 0, bool = false) noexcept {}
    inline void push(const T&) noexcept {}
    inline void push(T&&) noexcept {}
    inline void clear() noexcept {}
    inline void print(std::ostream& = std::cerr) const noexcept {}
};

#endif

} // namespace ggbp::stats



#endif
