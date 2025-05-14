#ifndef NTRU_POLYMODQ_HPP
#define NTRU_POLYMODQ_HPP

#include <vector>
#include <iostream>
#include <stdexcept>
#include "convpoly.hpp" // For ConvPoly return type

class PolyModQ {
public:
    std::vector<int> coef;
    int degree;
    int q;
    void trim();
    bool is_zero() const;
    std::pair<PolyModQ, PolyModQ> divmod(const PolyModQ& divisor) const;
    PolyModQ(const std::vector<int>& coef_ = {0}, int q_ = 3);

    bool operator==(const PolyModQ& other) const;
    bool operator!=(const PolyModQ& other) const;
    PolyModQ inverse(int Nval = 0, bool debug = false) const;
    PolyModQ operator-() const;
    PolyModQ operator+(const PolyModQ& other) const;
    PolyModQ operator+(int other) const;
    friend PolyModQ operator+(int lhs, const PolyModQ& rhs);

    PolyModQ operator-(const PolyModQ& other) const;
    PolyModQ operator-(int other) const;
    friend PolyModQ operator-(int lhs, const PolyModQ& rhs);

    PolyModQ operator*(const PolyModQ& other) const;
    PolyModQ operator*(int other) const;
    friend PolyModQ operator*(int lhs, const PolyModQ& rhs);

    friend std::ostream& operator<<(std::ostream& os, const PolyModQ& p);

    ConvPoly centerLift() const;
};

#endif // NTRU_POLYMODQ_HPP
// This code defines the PolyModQ class, which represents a polynomial modulo q.
// It includes methods for polynomial arithmetic (addition, subtraction, multiplication),
// equality checks, and a method for center lifting the polynomial.