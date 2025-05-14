#ifndef NTRU_CONVMODQ_HPP
#define NTRU_CONVMODQ_HPP

#include "polymodq.hpp"
#include "convpoly.hpp"
#include <vector>
#include <iostream>

// Forward declaration of invModQ (if it's in utils)
int invModQ(int d, int q);

class ConvModQ : public PolyModQ {
public:
    int N;

    ConvModQ();
    ConvModQ(const std::vector<int>& coef, int q = 3, int N_ = -1);
    friend std::ostream& operator<<(std::ostream& os, const ConvModQ& p);
    ConvPoly centerLift() const;
    ConvModQ operator+(const ConvModQ& other) const;
    ConvModQ operator-(const ConvModQ& other) const;
    ConvModQ operator*(const ConvModQ& other) const;
    ConvModQ operator*(const ConvPoly& other) const; // If you implement this
    ConvModQ operator*(int other) const;
    friend ConvModQ operator*(int lhs, const ConvModQ& rhs);

    ConvModQ operator/(const ConvModQ& other) const;
    ConvModQ operator/(int other) const;

    ConvModQ modQ(int q_) const;
    ConvModQ inverse(int Nval = -1, bool debug = false) const;
};

#endif // NTRU_CONVMODQ_HPP