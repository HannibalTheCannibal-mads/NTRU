#include "ntru/convpoly.hpp"
#include <stdexcept>
#include <vector>
#include <iostream>

// Constructor
ConvPoly::ConvPoly(std::vector<int> coef_, int N_) {
    if (N_ == -1) {
        N = static_cast<int>(coef_.size());
        coef = coef_;
    } else {
        N = N_;
        coef = coef_;
        coef.resize(N, 0);
    }
}

// String representation
std::ostream& operator<<(std::ostream& os, const ConvPoly& p) {
    os << "ConvPoly[";
    for (size_t i = 0; i < p.coef.size(); ++i) {
        os << p.coef[i];
        if (i != p.coef.size() - 1) os << ", ";
    }
    os << "]";
    return os;
}

// Addition
ConvPoly ConvPoly::operator+(const ConvPoly& other) const {
    if (N == other.N) {
        std::vector<int> result(N);
        for (int i = 0; i < N; ++i)
            result[i] = coef[i] + other.coef[i];
        return ConvPoly(result);
    }
    throw std::invalid_argument("ConvPoly: N mismatch in addition");
}

ConvPoly ConvPoly::operator+(int other) const {
    std::vector<int> result = coef;
    result[0] += other;
    return ConvPoly(result, N);
}

ConvPoly operator+(int lhs, const ConvPoly& rhs) {
    return rhs + lhs;
}

// Equality
bool ConvPoly::operator==(const ConvPoly& other) const {
    return coef == other.coef;
}

bool ConvPoly::operator!=(const ConvPoly& other) const {
    return !(*this == other);
}

// Negation
ConvPoly ConvPoly::operator-() const {
    std::vector<int> result(N);
    for (int i = 0; i < N; ++i)
        result[i] = -coef[i];
    return ConvPoly(result, N);
}

// Subtraction
ConvPoly ConvPoly::operator-(const ConvPoly& other) const {
    return *this + (-other);
}

ConvPoly ConvPoly::operator-(int other) const {
    return *this + (-other);
}

ConvPoly operator-(int lhs, const ConvPoly& rhs) {
    return (-rhs) + lhs;
}

// Multiplication (convolution)
ConvPoly ConvPoly::operator*(const ConvPoly& other) const {
    if (N == other.N) {
        std::vector<int> result(N, 0);
        for (int k = 0; k < N; ++k) {
            int s = 0;
            for (int i = 0; i < N; ++i)
                s += coef[i] * other.coef[(k - i + N) % N];
            result[k] = s;
        }
        return ConvPoly(result, N);
    }
    throw std::invalid_argument("ConvPoly: N mismatch in multiplication");
}

ConvPoly ConvPoly::operator*(int other) const {
    std::vector<int> result(N);
    for (int i = 0; i < N; ++i)
        result[i] = coef[i] * other;
    return ConvPoly(result, N);
}

ConvPoly operator*(int lhs, const ConvPoly& rhs) {
    return rhs * lhs;
}
// This code defines the ConvPoly class, which represents a polynomial with coefficients in a finite field.
// It includes methods for polynomial arithmetic (addition, subtraction, multiplication),
// equality checks, and a method for negation. The class also provides a string representation for easy debugging.
// The ConvPoly class is used in the context of NTRU encryption, where polynomials are manipulated
// in a finite field defined by the parameters of the NTRU algorithm.