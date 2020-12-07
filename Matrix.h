#if !defined(MATRIX_H)
#define MATRIX_H

#include <memory>
#include <iostream>
#include <cstdarg> //va_start, va_end
#include <type_traits> // is_pointer
#include <cmath>

#include "util.h"
#include "Range.h"

template <typename T, size_t _ROWS, size_t _COLS>
class Mat
{
private:
    using Type = T;
    T raw_[_ROWS * _COLS]; // raw data in row major order

public:
    constexpr static const size_t SIZE = _ROWS * _COLS;
    static const size_t ROWS = _ROWS;
    static const size_t COLS = _COLS;

    inline Mat() : raw_() {}

    inline Mat(T a11, ...) : raw_()
    {
        va_list vl;
        va_start(vl, a11);
        raw_[0] = a11;
        for (auto &&i : Range(1, SIZE))
            raw_[i] = va_arg(vl, T);
        va_end(vl);
    }

    inline T * as_raw_mut() {
        return raw_;
    }

    inline const T * as_raw() const {
        return raw_;
    }

    static Mat identity()
    {
        auto m = Mat();
        for (auto &&i : Range(COLS))
        {
            m(i, i) = 1;
        }
        return m;
    }

    inline Mat transpose()
    {
        auto m1 = *this;
        auto m2 = Mat();
        for (auto &&i : Range(ROWS))
            for (auto &&j : Range(COLS))
                m2(j, i) = m1(i, j);
        return m2;
    };

    inline T &operator()(size_t i, size_t j)
    {
        if (i >= ROWS)
            PANIC("Invalid Matrix index, tried to access row: ", i);
        if (j >= COLS)
            PANIC("Invalid Matrix index, tried to access column: ", j);
        return raw_[j + i * COLS];
    }

    template <size_t START_ROW, size_t STOP_ROW, size_t START_COL, size_t STOP_COL, size_t STEP_ROW = 1, size_t STEP_COL = 1, typename S = T>
    typename std::enable_if<!std::is_pointer<typename to_raw_pointer<S>::Raw>::value,
                            Mat<T *, (STOP_ROW - START_ROW + 1) / STEP_ROW, (STOP_COL - START_COL + 1) / STEP_COL>>::type
    slice()
    {
        auto m = Mat<T *, (STOP_ROW - START_ROW + 1) / STEP_ROW, (STOP_COL - START_COL + 1) / STEP_COL>();
        size_t i_m = 0;
        for (auto &&i : Range(START_ROW, STOP_ROW + 1, STEP_ROW))
        {
            size_t j_m = 0;
            for (auto &&j : Range(START_COL, STOP_COL + 1, STEP_COL))
            {
                m(i_m, j_m) = &this->operator()(i, j);
                j_m += 1;
            }
            i_m += 1;
        }
        return m;
    }

    inline T
    operator()(size_t i, size_t j) const
    {
        if (i >= ROWS)
            PANIC("Invalid Matrix index, tried to access row: ", i);
        if (j >= COLS)
            PANIC("Invalid Matrix index, tried to access column: ", j);
        return raw_[j + i * COLS];
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

    /// Indices of elements in column `column`
    inline static Range column_indices(size_t column)
    {
        if (column >= COLS)
            PANIC("Invalid Matrix index, tried to access column: ", column);
        return Range(column, SIZE, COLS);
    }

    /// Indices of elements in row `row`
    inline static Range row_indices(size_t row)
    {
        if (row >= ROWS)
            PANIC("Invalid Matrix index, tried to access row: ", row);
        auto start = row * COLS;
        return Range(start, start + COLS);
    }

    inline Mat<T, ROWS, 1> column(size_t column)
    {
        auto m = Mat<T, ROWS, 1>();
        auto i = 0;
        for (auto &&idx : this->column_indices(column))
        {
            m[i] = (*this)[idx];
            i++;
        }
        return m;
    }

    inline Mat<T, 1, COLS> row(size_t row)
    {
        auto m = Mat<T, 1, COLS>();
        auto i = 0;
        for (auto &&idx : this->row_indices(row))
        {
            m[i] = (*this)[idx];
            i++;
        }
        return m;
    }

    // Matrix addition
    template <typename S = T>
    inline typename std::enable_if<!std::is_pointer<typename to_raw_pointer<S>::Raw>::value,
                                   Mat<T, ROWS, COLS>>::type
    operator+(const Mat<T, ROWS, COLS> &other) const
    {
        auto m3 = Mat<T, ROWS, COLS>();
        for (auto &&i : Range(ROWS))
        {
            for (auto &&j : Range(COLS))
            {
                m3(i, j) = (*this)(i, j) + other(i, j);
            }
        }
        return m3;
    }

    // Matrix subtraction
    template <typename S = T>
    inline typename std::enable_if<!std::is_pointer<typename to_raw_pointer<S>::Raw>::value,
                                   Mat<T, ROWS, COLS>>::type
    operator-(const Mat<T, ROWS, COLS> &other) const
    {
        auto m3 = Mat<T, ROWS, COLS>();
        for (auto &&i : Range(ROWS))
        {
            for (auto &&j : Range(COLS))
            {
                m3(i, j) = (*this)(i, j) - other(i, j);
            }
        }
        return m3;
    }

    template <typename U>
    Mat<U, ROWS, COLS> map(U f(T)) const
    {
        auto m = Mat<U, ROWS, COLS>();
        for (auto &&i : Range(ROWS))
        {
            for (auto &&j : Range(COLS))
            {
                m(i, j) = f((*this)(i, j));
            }
        }
        return m;
    }

    template <typename S = T, typename Ptr = to_raw_pointer<S>>
    inline typename std::enable_if<std::is_pointer<typename Ptr::Raw>::value,
                                   Mat<typename Ptr::Element, ROWS, COLS>>::type
    to_owned() const
    {
        return this->map<typename Ptr::Element>([](T x) { return *x; });
    }

    template <typename S = T>
    typename std::enable_if<std::is_pointer<typename to_raw_pointer<S>::Raw>::value, String>::type show() const
    {
        auto s = String();
        for (auto &&i : Range(ROWS))
        {
            for (auto &&j : Range(COLS))
                s.append(std::to_string(*((*this)(i, j))) + " ");
            s.append("\n");
        }
        s.pop_back();
        return s;
    };

    template <typename S = T>
    typename std::enable_if<!std::is_pointer<typename to_raw_pointer<S>::Raw>::value, String>::type show() const
    {
        auto s = String();
        for (auto &&i : Range(ROWS))
        {
            for (auto &&j : Range(COLS))
                s.append(std::to_string((*this)(i, j)) + " ");
            s.append("\n");
        }
        s.pop_back();
        return s;
    };

    // compound matrix multiplication with static matrix no pointer
    template <typename S = T>
    typename std::enable_if<!std::is_pointer<typename to_raw_pointer<S>::Raw>::value, Mat<T, ROWS, COLS> &>::type
    operator*=(const Mat<T, ROWS, COLS> &other) 
    {
        *this = (*this) * other;
        return *this;
    }

    // compound matrix division with scalar
    template <typename S = T>
    typename std::enable_if<!std::is_pointer<typename to_raw_pointer<S>::Raw>::value, Mat<T, ROWS, COLS> &>::type
    operator/=(const T &other) 
    {
        for (auto &&i : Range(SIZE))
        {
            this->raw_[i] /= other;
        }
        return *this;
    }

    template <typename S = T, typename Ptr = to_raw_pointer<S>>
    typename std::enable_if<std::is_pointer<typename Ptr::Raw>::value, Mat<T, ROWS, COLS> &>::type
    operator=(const Mat<typename Ptr::Element, ROWS, COLS> other)
    {
        for (auto &&i : Range(SIZE))
        {
            *(this->raw_[i]) = other[i];
        }
        return *this;
    }

    template <typename S = T, typename Ptr = to_raw_pointer<S>>
    typename std::enable_if<std::is_pointer<typename Ptr::Raw>::value, Mat<T, ROWS, COLS> &>::type
    operator=(const typename Ptr::Element scalar)
    {
        for (auto &&i : Range(SIZE))
        {
            *(this->raw_[i]) = scalar;
        }
        return *this;
    }

    // Frobenius norm
    template <typename S = T>
    inline typename std::enable_if<!std::is_pointer<typename to_raw_pointer<S>::Raw>::value, T>::type
    frobenius_norm() const
    {
        auto sum = 0;
        for (auto &&i : Range(ROWS))
        {
            for (auto &&j : Range(COLS))
            {
                sum += (*this)(i, j);
            }
        }
        return sqrt(sum);
    }    

    inline auto begin() const { return raw_; }
    inline auto end() const { return raw_; }
};

// Scalar multiplication
template <typename T, size_t ROWS, size_t COLS>
inline typename std::enable_if<!std::is_pointer<typename to_raw_pointer<T>::Raw>::value,
                               Mat<T, ROWS, COLS>>::type
operator*(const T factor, const Mat<T, ROWS, COLS> &mat)
{
    return mat.map([&factor](T x) { return factor * x; });
    /*
    auto m3 = Mat<T, ROWS, COLS>();
    for (auto &&i : Range(ROWS))
    {
        for (auto &&j : Range(COLS))
        {
            m3(i, j) = factor * mat(i, j);
        }
    }
    return m3;
    */
}

// Scalar division
template <typename T, size_t ROWS, size_t COLS>
inline typename std::enable_if<!std::is_pointer<typename to_raw_pointer<T>::Raw>::value,
                               Mat<T, ROWS, COLS>>::type
operator/(const Mat<T, ROWS, COLS> &mat, const T divisor)
{
    auto m3 = Mat<T, ROWS, COLS>();
    for (auto &&i : Range(ROWS))
    {
        for (auto &&j : Range(COLS))
        {
            m3(i, j) = mat(i, j) / divisor;
        }
    }
    return m3;
}

// Matrix multiplication
template <typename T, size_t ROWS, size_t COLS, size_t COLS2>
inline typename std::enable_if<!std::is_pointer<typename to_raw_pointer<T>::Raw>::value,
                               Mat<T, ROWS, COLS2>>::type
operator*(const Mat<T, ROWS, COLS> &self, const Mat<T, COLS, COLS2> &other)
{
    auto m3 = Mat<T, ROWS, COLS2>();
    for (auto &&i : Range(ROWS))
    {
        auto j = 0;
        for (auto &&j : Range(COLS2))
        {
            auto sum = 0.0;
            for (auto &&k : Range(COLS))
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
template <typename T, size_t VEC_LENGTH>
inline typename std::enable_if<!std::is_pointer<typename to_raw_pointer<T>::Raw>::value,
                               T>::type
operator*(const Mat<T, 1, VEC_LENGTH> &self, const Mat<T, VEC_LENGTH, 1> &other)
{
    auto scal_prod = 0.0;
    for (auto &&i : Range(VEC_LENGTH))
    {
        scal_prod += self[static_cast<size_t>(i)] * other[static_cast<size_t>(i)];
    }
    return scal_prod;
}

using Mat4x4 = Mat<double, 4, 4>;

#endif // MATRIX_H
