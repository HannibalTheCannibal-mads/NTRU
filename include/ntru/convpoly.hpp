#ifndef NTRU_CONVPOLY_HPP
#define NTRU_CONVPOLY_HPP

#include <vector>
#include <iostream>

class ConvPoly {
public:
    std::vector<int> coef;
    int N;

    // Constructors
    ConvPoly(std::vector<int> coef_ = {0}, int N_ = -1);

    // String representation
    friend std::ostream& operator<<(std::ostream& os, const ConvPoly& p);

    // Addition
    ConvPoly operator+(const ConvPoly& other) const;
    ConvPoly operator+(int other) const;
    friend ConvPoly operator+(int lhs, const ConvPoly& rhs);

    // Equality
    bool operator==(const ConvPoly& other) const;
    bool operator!=(const ConvPoly& other) const;

    // Negation
    ConvPoly operator-() const;

    // Subtraction
    ConvPoly operator-(const ConvPoly& other) const;
    ConvPoly operator-(int other) const;
    friend ConvPoly operator-(int lhs, const ConvPoly& rhs);

    // Multiplication (convolution)
    ConvPoly operator*(const ConvPoly& other) const;
    ConvPoly operator*(int other) const;
    friend ConvPoly operator*(int lhs, const ConvPoly& rhs);
};

#endif // NTRU_CONVPOLY_HPP
