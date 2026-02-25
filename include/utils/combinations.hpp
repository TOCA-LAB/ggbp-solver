// GPT - implementation

#ifndef COMBINATIONS_HPP
#define COMBINATIONS_HPP

#include <vector>
#include <cstddef>
#include <iterator>
#include <type_traits>


namespace ggbp::utils {

template<typename Container>
class Combinations {
public:
    using container_type = Container;
    using value_type = typename Container::value_type;
    using size_type = std::size_t;
    using subset_type = std::vector<value_type>;

    Combinations(const Container& c, size_type k) noexcept
        : data_(&c), n_(c.size()), k_(k) {}

    class iterator {
    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = subset_type;
        using difference_type = std::ptrdiff_t;
        using pointer = const subset_type*;
        using reference = const subset_type&;

        iterator() noexcept : data_(nullptr), n_(0), k_(0), end_(true) {}

        iterator(const Container* data, size_type n, size_type k, bool end_flag)
            : data_(data), n_(n), k_(k), end_(false)
        {
            if (k_ > n_) { end_ = true; return; }
            if (k_ == 0) {
                if (end_flag) end_ = true;
                else cur_.clear();
                return;
            }
            if (end_flag) {
                end_ = true;
                return;
            }
            cur_idx_.resize(k_);
            for (size_type i = 0; i < k_; ++i) cur_idx_[i] = i;
            build_cur_from_idx();
        }

        reference operator*() const noexcept { return cur_; }
        pointer operator->() const noexcept { return &cur_; }

        iterator& operator++() noexcept {
            advance();
            return *this;
        }

        iterator operator++(int) noexcept {
            iterator tmp = *this;
            advance();
            return tmp;
        }

        bool operator==(const iterator& other) const noexcept {
            if (end_ && other.end_) return true;
            if (end_ != other.end_) return false;
            // compare indices
            return cur_idx_ == other.cur_idx_;
        }

        bool operator!=(const iterator& other) const noexcept {
            return !(*this == other);
        }

    private:
        const Container* data_{nullptr};
        size_type n_{0};
        size_type k_{0};
        std::vector<size_type> cur_idx_;
        subset_type cur_;
        bool end_{true};

        void build_cur_from_idx() {
            cur_.clear();
            cur_.reserve(k_);
            for (size_type i = 0; i < k_; ++i) cur_.push_back((*data_)[cur_idx_[i]]);
        }

        void advance() noexcept {
            if (end_) return;
            if (k_ == 0) { end_ = true; return; }

            // find rightmost index that can be incremented
            int i = static_cast<int>(k_) - 1;
            for (; i >= 0; --i) {
                if (cur_idx_[static_cast<size_type>(i)] != static_cast<size_type>(i) + n_ - k_)
                    break;
            }
            if (i < 0) { end_ = true; return; }

            ++cur_idx_[static_cast<size_type>(i)];
            for (size_type j = static_cast<size_type>(i) + 1; j < k_; ++j)
                cur_idx_[j] = cur_idx_[j - 1] + 1;

            build_cur_from_idx();
        }
    };

    iterator begin() const noexcept { return iterator(data_, n_, k_, false); }
    iterator end()   const noexcept { return iterator(data_, n_, k_, true); }

private:
    const Container* data_;
    size_type n_;
    size_type k_;
};

} // namespace ggbp

#endif // GGBP_COMBINATIONS_VIEW_HPP

