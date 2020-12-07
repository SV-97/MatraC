#if !defined(RANGE_H)
#define RANGE_H

#include <iostream>
#include <cstdint> // int64_t

#include "util.h"

class Range
{
    /*
    Example:
        >>>auto r = Range(4, 15);
        >>>for (auto &&i : Range(4, 10))
        >>>{
        >>>    std::cout << i << " " << r[i] << std::endl;
        >>>}
        4 8
        5 9
        6 10
        7 11
        8 12
        9 13
    */
private:
    const int64_t START_;
    int64_t i_;
    const int64_t STOP_;
    const int64_t STEP_;

public:
    inline Range(int64_t stop) : Range(0, stop, 1) {}

    inline Range(int64_t start, int64_t stop, int64_t step = 1) : START_(start), i_(start), STOP_(stop), STEP_(step)
    {
        if (stop < start)
            PANIC("Invalid range: ", this->show());
    }

    inline String show() const { return string_format("[%d,%d) by %d", START_, STOP_, STEP_); }

    inline bool done() const
    {
        return !(START_ <= i_ && i_ < STOP_);
    }

    inline bool operator!=(const Range &rhs) const { return !this->done(); }
    inline int64_t &operator*() { return i_; }

    inline Range &operator++()
    {
        i_ += STEP_;
        return *this;
    }

    inline Range operator++(int i)
    {
        i_ += STEP_; // yes this isn't really a post increment
        return *this;
    }

    inline Range &operator--()
    {
        i_ -= STEP_;
        return *this;
    }

    inline Range operator--(int i)
    {
        i_ -= STEP_;
        return *this;
    }

    inline int64_t operator[](int64_t i) const
    {
        auto x = i_ + i * STEP_;
        if (x > STOP_)
            PANIC("Invalid index in Range: ", this->show());
        return x;
    }

    inline Range operator[](std::pair<int64_t, int64_t> idx) const
    {
        auto s = i_ + idx.first;
        auto e = i_ + idx.second;
        if (s > STOP_ || e > STOP_)
            PANIC("Invalid index in Range: ", this->show(), ' ', idx.first, ' ', idx.second);
        return Range(s, e, STEP_);
    }

    inline auto begin() const { return *this; }
    inline auto end() const { return *this; }
    inline auto rbegin() const { return *this; }
    inline auto rend() const { return *this; }
};

#endif // RANGE_H
