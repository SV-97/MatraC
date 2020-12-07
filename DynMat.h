#if !defined(DYN_MATRIX_H)
#define DYN_MATRIX_H

#include <memory>
#include <iostream>
#include <cstdarg> // va_start, va_end
#include <type_traits> // is_pointer
#include <vector>
#include <unordered_map>
#include <cassert> // assert
#include <cmath>

#include "util.h"
#include "Range.h"

#include "Matrix.h"

namespace internal
{
    /* Membuf has to be a template of kind * -> *, that if instantiated with T has a constructor of
    type (size_t, size_t) -> MemBuf<T>, as well as implementations of T operator[](size_t) const
        / T& operator[](size_t)
    */
    template <typename T, template <class> typename MemBuf>
    class AbstractDynMat
    {
    private:
        using Type = T;

        MemBuf<T> raw_; // raw data in row major order

    public:
        const size_t SIZE;
        const size_t ROWS_;
        const size_t COLS_;

        inline AbstractDynMat(size_t rows, size_t cols) : SIZE(rows * cols), ROWS_(rows), COLS_(cols), raw_(rows, cols) {}

        inline AbstractDynMat(size_t rows, size_t cols, T a11, ...) : AbstractDynMat(rows, cols)
        {
            raw_ = MemBuf<T>(rows, cols);
            va_list vl;
            va_start(vl, a11);
            raw_[0] = a11;
            for (auto &&i : Range(1, SIZE))
                raw_[i] = va_arg(vl, T);
            va_end(vl);
        }

        static AbstractDynMat identity(size_t rows, size_t cols)
        {
            auto m = AbstractDynMat(rows, cols);
            for (auto &&i : Range(cols))
            {
                m(i, i) = 1;
            }
            return m;
        }

        inline AbstractDynMat transpose()
        {
            auto m1 = *this;
            auto m2 = AbstractDynMat(COLS_, ROWS_);
            for (auto &&i : Range(ROWS_))
                for (auto &&j : Range(COLS_))
                    m2(j, i) = m1(i, j);
            return m2;
        };

        inline T &operator()(size_t i, size_t j)
        {
            if (i >= ROWS_)
                PANIC("Invalid Matrix index, tried to access row: ", i);
            if (j >= COLS_)
                PANIC("Invalid Matrix index, tried to access column: ", j);
            return raw_[j + i * COLS_];
        }

        template <template <class> typename MemBufOut = MemBuf, typename S = T>
        typename std::enable_if<!std::is_pointer<typename to_raw_pointer<S>::Raw>::value, AbstractDynMat<T *, MemBufOut>>::type
        slice(size_t start_row, size_t stop_row, size_t start_col, size_t stop_col, size_t step_row = 1, size_t step_col = 1)
        {
            auto m = AbstractDynMat<T *, MemBufOut>((stop_row - start_row + 1) / step_row, (stop_col - start_col + 1) / step_col);
            size_t i_m = 0;
            for (auto &&i : Range(start_row, stop_row + 1, step_row))
            {
                size_t j_m = 0;
                for (auto &&j : Range(start_col, stop_col + 1, step_col))
                {
                    m(i_m, j_m) = &this->operator()(i, j);
                    j_m += 1;
                }
                i_m += 1;
            }
            return m;
        }

        inline T operator()(size_t i, size_t j) const
        {
            if (i >= ROWS_)
                PANIC("Invalid Matrix index, tried to access row: ", i);
            if (j >= COLS_)
                PANIC("Invalid Matrix index, tried to access column: ", j);
            return raw_[j + i * COLS_];
        }

        inline T &operator[](size_t i)
        {
            if (i >= SIZE)
                PANIC("Invalid Matrix index, tried to access index: ", i);
            return raw_[i];
        }

        inline T operator[](size_t i) const
        {
            if (i >= SIZE)
                PANIC("Invalid Matrix index, tried to access index: ", i);
            return raw_[i];
        }

        // Matrix addition
        template <template <class> typename MemBufOther = MemBuf, template <class> typename MemBufOut = MemBuf, typename S = T>
        inline typename std::enable_if<!std::is_pointer<typename to_raw_pointer<S>::Raw>::value,
                                       AbstractDynMat<T, MemBufOut>>::type
        operator+(const AbstractDynMat<T, MemBufOther> &other) const
        {
            auto m3 = AbstractDynMat<T, MemBufOut>(ROWS_, COLS_);
            for (auto &&i : Range(ROWS_))
            {
                for (auto &&j : Range(COLS_))
                {
                    m3(i, j) = (*this)(i, j) + other(i, j);
                }
            }
            return m3;
        }


        // Matrix addition with static
        template <size_t ROWS, size_t COLS, template <class> typename MemBufOther = MemBuf, template <class> typename MemBufOut = MemBuf, typename S = T>
        inline typename std::enable_if<!std::is_pointer<typename to_raw_pointer<S>::Raw>::value,
                                       AbstractDynMat<T, MemBufOut>>::type
        operator+(const Mat<T, ROWS, COLS> &other) const
        {
            assert(ROWS == ROWS_ && COLS == COLS_);
            auto m3 = AbstractDynMat<T, MemBufOut>(ROWS_, COLS_);
            m3 += other;
            return m3;
        }

        // compound matrix addition with static
        template <size_t ROWS, size_t COLS, template <class> typename MemBufOther = MemBuf, typename S = T, typename Ptr = to_raw_pointer<S>>
        typename std::enable_if<std::is_pointer<typename Ptr::Raw>::value, AbstractDynMat<T, MemBuf> &>::type operator+=(const Mat<typename Ptr::Element, ROWS, COLS> &other) 
        {
            if (ROWS != ROWS_ || COLS != COLS_) {
                PANIC("Incompatible matrix dimensions: ", ROWS_,'x', COLS_, " + ", ROWS, 'x', COLS)
            }
            for (auto &&i : Range(SIZE))
            {
                *(this->raw_[i]) += other[i];
            }
            return *this;
        }

        // Matrix subtraction
        template <template <class> typename MemBufOther = MemBuf, template <class> typename MemBufOut = MemBuf, typename S = T>
        inline typename std::enable_if<!std::is_pointer<typename to_raw_pointer<S>::Raw>::value,
                                       AbstractDynMat<T, MemBufOut>>::type
        operator-(const AbstractDynMat<T, MemBufOther> &other) const
        {
            assert(ROWS_ == other.ROWS_ && COLS_ == other.COLS_);
            auto m3 = AbstractDynMat<T, MemBufOut>(ROWS_, COLS_);
            for (auto &&i : Range(ROWS_))
            {
                for (auto &&j : Range(COLS_))
                {
                    m3(i, j) = (*this)(i, j) - other(i, j);
                }
            }
            return m3;
        }

        template <typename U, template <class> typename MemBufOut = MemBuf>
        AbstractDynMat<U, MemBufOut> map(U f(T)) const
        {
            auto m = AbstractDynMat<U, MemBufOut>(ROWS_, COLS_);
            for (auto &&i : Range(ROWS_))
            {
                for (auto &&j : Range(COLS_))
                {
                    m(i, j) = f((*this)(i, j));
                }
            }
            return m;
        }

        template <template <class> typename MemBufOut, typename S = T, typename Ptr = to_raw_pointer<S>>
        inline typename std::enable_if<std::is_pointer<typename Ptr::Raw>::value,
                                       AbstractDynMat<typename Ptr::Element, MemBufOut>>::type
        to_owned()
        {
            return this->map<typename Ptr::Element>([](T x) { return *x; });
        }

        template <typename S = T>
        typename std::enable_if<std::is_pointer<typename to_raw_pointer<S>::Raw>::value, String>::type show()
        {
            auto s = String();
            for (auto &&i : Range(ROWS_))
            {
                for (auto &&j : Range(COLS_))
                    s.append(std::to_string(*((*this)(i, j))) + " ");
                s.append("\n");
            }
            s.pop_back();
            return s;
        };

        template <typename S = T>
        typename std::enable_if<!std::is_pointer<typename to_raw_pointer<S>::Raw>::value, String>::type show()
        {
            auto s = String();
            for (auto &&i : Range(ROWS_))
            {
                for (auto &&j : Range(COLS_))
                    s.append(std::to_string((*this)(i, j)) + " ");
                s.append("\n");
            }
            s.pop_back();
            return s;
        };

        template <template <class> typename MemBufOther = MemBuf, typename S = T, typename Ptr = to_raw_pointer<S>>
        typename std::enable_if<std::is_pointer<typename Ptr::Raw>::value, AbstractDynMat<T, MemBuf> &>::type
        operator=(const AbstractDynMat<typename Ptr::Element, MemBufOther> other)
        {
            for (auto &&i : Range(SIZE))
            {
                *(this->raw_[i]) = other[i];
            }
            return *this;
        }

        template <typename S = T, typename Ptr = to_raw_pointer<S>>
        typename std::enable_if<std::is_pointer<typename Ptr::Raw>::value, AbstractDynMat<T, MemBuf> &>::type
        operator=(const typename Ptr::Element scalar)
        {
            for (auto &&i : Range(SIZE))
            {
                *(this->raw_[i]) = scalar;
            }
            return *this;
        }

        inline auto begin() const { return raw_.begin(); }
        inline auto end() const { return raw_.end(); }
    };

    // Scalar multiplication
    template <typename T, template <class> typename MemBuf, template <class> typename MemBufOut = MemBuf>
    inline typename std::enable_if<!std::is_pointer<typename to_raw_pointer<T>::Raw>::value,
                                   AbstractDynMat<T, MemBufOut>>::type
    operator*(const T factor, const AbstractDynMat<T, MemBuf> &mat)
    {
        auto m3 = AbstractDynMat<T, MemBufOut>(mat.ROWS_, mat.COLS_);
        for (auto &&i : Range(mat.ROWS_))
        {
            for (auto &&j : Range(mat.COLS_))
            {
                m3(i, j) = factor * mat(i, j);
            }
        }
        return m3;
    }

    // Matrix multiplication
    template <typename T, template <class> typename MemBuf, template <class> typename MemBufOther = MemBuf, template <class> typename MemBufOut = MemBuf>
    inline typename std::enable_if<!std::is_pointer<typename to_raw_pointer<T>::Raw>::value,
                                   AbstractDynMat<T, MemBufOut>>::type
    operator*(const AbstractDynMat<T, MemBuf> &self, const AbstractDynMat<T, MemBufOther> &other)
    {
        assert(self.COLS_ == other.ROWS_);
        auto m3 = AbstractDynMat<T, MemBufOut>(self.ROWS_, other.COLS_);
        for (auto &&i : Range(self.ROWS_))
        {
            auto j = 0;
            for (auto &&j : Range(other.COLS_))
            {
                auto sum = 0.0;
                for (auto &&k : Range(self.COLS_))
                {
                    sum += self(static_cast<size_t>(i), static_cast<size_t>(k)) *
                           other(static_cast<size_t>(k), static_cast<size_t>(j));
                }
                m3(i, j) = sum;
            }
        }
        return m3;
    }

    // Dot product
    template <typename T, template <class> typename MemBuf, template <class> typename MemBufOther = MemBuf>
    inline typename std::enable_if<!std::is_pointer<typename to_raw_pointer<T>::Raw>::value,
                                   T>::type
    operator*(const AbstractDynMat<T, MemBuf> &self, const AbstractDynMat<T, MemBufOther> &other)
    {
        assert(self.COLS_ == other.ROWS_ && self.ROWS_ == 1 && other.COLS_ == 1);
        auto scal_prod = 0.0;
        for (auto &&i : Range(self.COLS_))
        {
            scal_prod += self[static_cast<size_t>(i)] * other[static_cast<size_t>(i)];
        }
        return scal_prod;
    }
} // namespace internal

template <typename T>
class SparseBuffer
{
private:
    std::unordered_map<size_t, T> raw_;
    std::vector<size_t> potentially_zero_;
    size_t cnt_;
    size_t max_size_;
    size_t treshold_;

public:
    inline SparseBuffer(size_t rows, size_t cols) : raw_(), potentially_zero_(), cnt_(), max_size_(rows * cols), treshold_(ceil(0.05 * max_size_)) {}
    inline T operator[](size_t i) const
    {
        try
        {
            return raw_.at(i);
        }
        catch (const std::exception &e)
        {
            return T();
        }
    }
    inline T &operator[](size_t i)
    {
        cnt_++;
        if (cnt_ > treshold_)
        {
            auto zero = T();
            for (auto &&x : potentially_zero_)
            {
                if (raw_[x] == zero)
                {
                    raw_.erase(x);
                }
            }
            potentially_zero_.clear();
            cnt_ = 0;
        }
        potentially_zero_.push_back(i);
        return raw_[i];
    }

    inline auto begin() const { return raw_.begin(); }
    inline auto end() const { return raw_.end(); }
};

template <typename T>
class DynBuffer
{
private:
    std::vector<T> raw_;

public:
    inline DynBuffer(size_t rows, size_t cols) : raw_(rows * cols) {}
    inline T operator[](size_t i) const
    {
        return raw_[i];
    }
    inline T &operator[](size_t i)
    {
        return raw_[i];
    }
};

template <typename T>
using DynMat = internal::AbstractDynMat<T, DynBuffer>;

template <typename T>
using SparseMat = internal::AbstractDynMat<T, SparseBuffer>;

template <typename T>
typename std::enable_if<!std::is_pointer<typename to_raw_pointer<T>::Raw>::value, String>::type show_sparse(const SparseMat<T>& m)
{
    auto s = String();
    for (auto &&x : m)
    {
        // we know that the index x.first = j + i * COLS_ with row index i and column index j < COLS_
        auto j = x.first % m.COLS_;
        auto i = (x.first - j) / m.COLS_;
        assert(j + i*m.COLS_ == x.first);
        s.append(/*std::to_string(x.first) + */'(' + std::to_string(i) + ',' + std::to_string(j) + ')' + std::to_string(x.second) + " ");
    }
    return s;
};

#endif // DYN_MATRIX_H
