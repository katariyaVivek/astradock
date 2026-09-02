#pragma once

#include "math/vector3.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace astradock::math {

// Fixed-size dense matrix with compile-time row/column dimensions.
//
// M13 introduces the first covariance-bearing subsystem in AstraDock. A Kalman
// filter needs matrix shapes beyond 3x3 (6x6 translational covariances now,
// larger error-state covariances later) but does not justify a generic linear
// algebra dependency or dynamic allocation. This class provides the minimum
// dense arithmetic required by the estimation module:
//
//   addition / subtraction / scalar multiplication,
//   compatible-dimension multiplication,
//   transpose,
//   symmetrization P <- (P + P^T)/2,
//   diagonal extraction,
//   Cholesky factorization and SPD linear solves.
//
// Storage is row-major std::array<double, R*C>: no heap allocation, identical
// operation ordering between runs, therefore fully deterministic.
//
// Units and frames are intentionally not encoded here. Covariance blocks mix
// units (m^2, m^2/s^2, m^2/s); the containing estimation types document the
// state ordering that gives every row/column its meaning.
template <std::size_t R, std::size_t C>
class Matrix {
public:
    // Constructs an all-zero matrix.
    constexpr Matrix() noexcept : data_{} {}

    // Constructs from row-major values, e.g. Matrix<2, 2>({1, 0, 0, 1}).
    constexpr explicit Matrix(const std::array<double, R * C>& values) noexcept
        : data_(values) {}

    [[nodiscard]] static constexpr Matrix zero() noexcept {
        return Matrix{};
    }

    [[nodiscard]] static constexpr Matrix identity()
        requires(R == C)
    {
        Matrix result = zero();
        for (std::size_t i = 0; i < R; ++i) {
            result.data_[i * C + i] = 1.0;
        }
        return result;
    }

    [[nodiscard]] static constexpr std::size_t rows() noexcept { return R; }
    [[nodiscard]] static constexpr std::size_t cols() noexcept { return C; }

    // Element access with bounds checking, mirroring Matrix3 conventions.
    [[nodiscard]] constexpr double operator()(std::size_t row, std::size_t col) const {
        if (row >= R || col >= C) {
            throw std::out_of_range("Matrix index out of range");
        }
        return data_[row * C + col];
    }

    [[nodiscard]] constexpr double& operator()(std::size_t row, std::size_t col) {
        if (row >= R || col >= C) {
            throw std::out_of_range("Matrix index out of range");
        }
        return data_[row * C + col];
    }

    [[nodiscard]] constexpr Matrix operator+(const Matrix& other) const noexcept {
        Matrix result;
        for (std::size_t i = 0; i < R * C; ++i) {
            result.data_[i] = data_[i] + other.data_[i];
        }
        return result;
    }

    [[nodiscard]] constexpr Matrix operator-(const Matrix& other) const noexcept {
        Matrix result;
        for (std::size_t i = 0; i < R * C; ++i) {
            result.data_[i] = data_[i] - other.data_[i];
        }
        return result;
    }

    [[nodiscard]] constexpr Matrix operator*(double scalar) const noexcept {
        Matrix result;
        for (std::size_t i = 0; i < R * C; ++i) {
            result.data_[i] = data_[i] * scalar;
        }
        return result;
    }

    [[nodiscard]] constexpr Matrix<C, R> transpose() const noexcept {
        std::array<double, C * R> transposed{};
        for (std::size_t r = 0; r < R; ++r) {
            for (std::size_t c = 0; c < C; ++c) {
                transposed[c * R + r] = data_[r * C + c];
            }
        }
        return Matrix<C, R>(transposed);
    }

    // Compatible-dimension matrix product: this (R x C) times other (C x K).
    template <std::size_t K>
    [[nodiscard]] constexpr Matrix<R, K> operator*(const Matrix<C, K>& other) const noexcept {
        Matrix<R, K> result;
        for (std::size_t r = 0; r < R; ++r) {
            for (std::size_t k = 0; k < K; ++k) {
                double sum = 0.0;
                for (std::size_t c = 0; c < C; ++c) {
                    sum += data_[r * C + c] * other.data_[c * K + k];
                }
                result.data_[r * K + k] = sum;
            }
        }
        return result;
    }

    // Symmetrized copy: (M + M^T)/2. Only meaningful for square matrices.
    // Used after covariance operations where floating-point roundoff can
    // introduce tiny asymmetric components; see Lesson 013 for the numerical
    // rationale and why arbitrary entry clamping is forbidden instead.
    [[nodiscard]] constexpr Matrix symmetrized() const noexcept
        requires(R == C)
    {
        Matrix result;
        for (std::size_t r = 0; r < R; ++r) {
            for (std::size_t c = r; c < R; ++c) {
                const double average = 0.5 * (data_[r * C + c] + data_[c * C + r]);
                result.data_[r * C + c] = average;
                result.data_[c * C + r] = average;
            }
        }
        return result;
    }

    // Maximum absolute difference between M and M^T (symmetry defect measure).
    [[nodiscard]] double symmetry_error() const noexcept
        requires(R == C)
    {
        double worst = 0.0;
        for (std::size_t r = 0; r < R; ++r) {
            for (std::size_t c = r + 1; c < R; ++c) {
                const double diff = data_[r * C + c] - data_[c * C + r];
                const double magnitude = diff < 0.0 ? -diff : diff;
                if (magnitude > worst) {
                    worst = magnitude;
                }
            }
        }
        return worst;
    }

    // Copies the main diagonal into a column vector. Requires a square matrix.
    [[nodiscard]] constexpr Matrix<R, 1> diagonal() const noexcept
        requires(R == C)
    {
        Matrix<R, 1> result;
        for (std::size_t i = 0; i < R; ++i) {
            result.data_[i] = data_[i * C + i];
        }
        return result;
    }

    [[nodiscard]] bool all_finite() const noexcept {
        for (const double value : data_) {
            if (!std::isfinite(value)) {
                return false;
            }
        }
        return true;
    }

    // Raw row-major storage access (read-only), for tests and interop.
    [[nodiscard]] constexpr const std::array<double, R * C>& data() const noexcept {
        return data_;
    }

template <std::size_t R2, std::size_t C2>
friend class Matrix;

private:
    std::array<double, R * C> data_;
};

template <std::size_t R>
using ColVector = Matrix<R, 1>;

template <std::size_t R, std::size_t C>
[[nodiscard]] constexpr Matrix<R, C> operator*(
    double scalar,
    const Matrix<R, C>& matrix) noexcept {
    return matrix * scalar;
}

// Elementwise approximate equality for scalars is provided by
// math/vector3.hpp (approximately_equal); reuse it here.

// Maximum absolute elementwise difference between two same-shape matrices.
template <std::size_t R, std::size_t C>
[[nodiscard]] double max_abs_difference(
    const Matrix<R, C>& lhs,
    const Matrix<R, C>& rhs) noexcept {
    double worst = 0.0;
    for (std::size_t r = 0; r < R; ++r) {
        for (std::size_t c = 0; c < C; ++c) {
            const double diff = lhs(r, c) - rhs(r, c);
            const double magnitude = diff < 0.0 ? -diff : diff;
            if (magnitude > worst) {
                worst = magnitude;
            }
        }
    }
    return worst;
}

// Extracts one column of a matrix as a column vector.
template <std::size_t R, std::size_t C>
[[nodiscard]] constexpr Matrix<R, 1> extract_column(
    const Matrix<R, C>& matrix,
    std::size_t col) {
    if (col >= C) {
        throw std::out_of_range("Matrix column index out of range");
    }
    Matrix<R, 1> result;
    for (std::size_t r = 0; r < R; ++r) {
        result(r, 0) = matrix(r, col);
    }
    return result;
}

// Cholesky factorization A = L * L^T of a symmetric positive definite matrix.
// Returns the lower-triangular factor L. Throws std::domain_error when the
// input contains non-finite values and std::runtime_error when a non-positive
// pivot reveals that A is not positive definite.
//
// This is the numerically stable workhorse behind every covariance solve in
// the M13 filter: innovation-covariance solves S x = b and PSD diagnostics.
// A failed factorization is reported loudly rather than silently repaired.
template <std::size_t N>
[[nodiscard]] Matrix<N, N> cholesky_lower(const Matrix<N, N>& a) {
    if (!a.all_finite()) {
        throw std::domain_error("Cholesky input must contain only finite values");
    }

    Matrix<N, N> lower;
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            double sum = a(i, j);
            for (std::size_t k = 0; k < j; ++k) {
                sum -= lower(i, k) * lower(j, k);
            }
            if (i == j) {
                if (!(sum > 0.0)) {
                    throw std::runtime_error(
                        "Cholesky factorization requires a positive definite matrix");
                }
                lower(i, j) = std::sqrt(sum);
            } else {
                lower(i, j) = sum / lower(j, j);
            }
        }
    }
    return lower;
}

// Solves the SPD system A x = b via Cholesky factorization without forming an
// explicit inverse. Forward substitution with L followed by back substitution
// with L^T costs O(N^2) per right-hand side and is more stable than computing
// A^{-1} explicitly.
template <std::size_t N>
[[nodiscard]] Matrix<N, 1> solve_spd(const Matrix<N, N>& a, const Matrix<N, 1>& b) {
    if (!b.all_finite()) {
        throw std::domain_error("SPD solve right-hand side must contain only finite values");
    }

    const Matrix<N, N> lower = cholesky_lower(a);

    // Forward substitution: L y = b.
    Matrix<N, 1> y;
    for (std::size_t i = 0; i < N; ++i) {
        double sum = b(i, 0);
        for (std::size_t k = 0; k < i; ++k) {
            sum -= lower(i, k) * y(k, 0);
        }
        y(i, 0) = sum / lower(i, i);
    }

    // Back substitution: L^T x = y.
    Matrix<N, 1> x;
    for (std::size_t ii = 0; ii < N; ++ii) {
        const std::size_t i = N - 1 - ii;
        double sum = y(i, 0);
        for (std::size_t k = i + 1; k < N; ++k) {
            sum -= lower(k, i) * x(k, 0);
        }
        x(i, 0) = sum / lower(i, i);
    }
    return x;
}

// Solves A X = B for several stacked right-hand sides (columns of B).
// Equivalent to applying solve_spd column-by-column.
template <std::size_t N, std::size_t M>
[[nodiscard]] Matrix<N, M> solve_spd_multi(
    const Matrix<N, N>& a,
    const Matrix<N, M>& b) {
    Matrix<N, M> result;
    for (std::size_t col = 0; col < M; ++col) {
        const Matrix<N, 1> rhs_column = extract_column(b, col);
        const Matrix<N, 1> solution_column = solve_spd(a, rhs_column);
        for (std::size_t row = 0; row < N; ++row) {
            result(row, col) = solution_column(row, 0);
        }
    }
    return result;
}

// Converts a Vector3 into a 3x1 column vector (frame semantics unchanged).
[[nodiscard]] constexpr Matrix<3, 1> to_column(const Vector3& vector) noexcept {
    return Matrix<3, 1>(std::array<double, 3>{vector.x(), vector.y(), vector.z()});
}

// Converts a 3x1 column vector back into a Vector3 (frame semantics unchanged).
[[nodiscard]] constexpr Vector3 to_vector3(const Matrix<3, 1>& column) noexcept {
    return Vector3{column(0, 0), column(1, 0), column(2, 0)};
}

}  // namespace astradock::math
