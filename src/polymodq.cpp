#include "ntru/polymodq.hpp"
#include "ntru/convpoly.hpp"
#include <vector>
#include <iostream>
#include <algorithm>
#include <stdexcept>

// Constructor
PolyModQ::PolyModQ(const std::vector<int>& coef_, int q_) : q(q_) {
    std::vector<int> temp(coef_.size());
    std::transform(coef_.begin(), coef_.end(), temp.begin(),
                   [this](int x) { return ((x % q + q) % q); });

    int index = 0;
    for (int i = static_cast<int>(temp.size()) - 1; i >= 0; --i) {
        if (temp[i] != 0) {
            index = static_cast<int>(temp.size()) - i - 1;
            break;
        }
    }
    coef.assign(temp.begin(), temp.end() - index);
    degree = static_cast<int>(coef.size()) - 1;
}

// Equality
bool PolyModQ::operator==(const PolyModQ& other) const {
    return degree == other.degree && coef == other.coef;
}

bool PolyModQ::operator!=(const PolyModQ& other) const {
    return !(*this == other);
}

// ostream operator<<
std::ostream& operator<<(std::ostream& os, const PolyModQ& p) {
    os << "PolyModQ([";
    for (size_t i = 0; i < p.coef.size(); ++i) {
        os << p.coef[i];
        if (i != p.coef.size() - 1) os << ", ";
    }
    os << "], " << p.q << ")";
    return os;
}

// Negation
PolyModQ PolyModQ::operator-() const {
    std::vector<int> neg(coef.size());
    for (size_t i = 0; i < coef.size(); ++i)
        neg[i] = ((-coef[i] % q) + q) % q;
    return PolyModQ(neg, q);
}

// Addition
PolyModQ PolyModQ::operator+(const PolyModQ& other) const {
    if (q != other.q)
        throw std::invalid_argument("PolyModQ: modulus q mismatch in addition");

    if (degree > other.degree) {
        std::vector<int> padded = other.coef;
        padded.resize(coef.size(), 0);
        std::vector<int> result(coef.size());
        for (size_t i = 0; i < coef.size(); ++i)
            result[i] = (coef[i] + padded[i]) % q;
        return PolyModQ(result, q);
    } else if (degree < other.degree) {
        return other + *this;
    } else {
        std::vector<int> result(coef.size());
        for (size_t i = 0; i < coef.size(); ++i)
            result[i] = (coef[i] + other.coef[i]) % q;
        return PolyModQ(result, q);
    }
}

PolyModQ PolyModQ::operator+(int other) const {
    std::vector<int> result = coef;
    if (!result.empty())
        result[0] = (result[0] + other) % q;
    return PolyModQ(result, q);
}

PolyModQ operator+(int lhs, const PolyModQ& rhs) {
    return rhs + lhs;
}

// Subtraction
PolyModQ PolyModQ::operator-(const PolyModQ& other) const {
    return *this + (-other);
}

PolyModQ PolyModQ::operator-(int other) const {
    return *this + (-other);
}

PolyModQ operator-(int lhs, const PolyModQ& rhs) {
    return (-rhs) + lhs;
}
// Cyclic convolution for NTRU: multiplies two polynomials mod x^N-1 and mod q
static std::vector<int> cyclic_convolution(const std::vector<int>& a, const std::vector<int>& b, int q) {
    int N = a.size();
    std::vector<int> result(N, 0);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            int idx = (i + j) % N;
            result[idx] = (result[idx] + a[i] * b[j]) % q;
            if (result[idx] < 0) result[idx] += q;
        }
    }
    return result;
}
// Multiplication: cyclic convolution for NTRU
PolyModQ PolyModQ::operator*(const PolyModQ& other) const {
    if (q != other.q)
        throw std::invalid_argument("PolyModQ: modulus q mismatch in multiplication");
    if (coef.size() != other.coef.size())
        throw std::invalid_argument("PolyModQ: size mismatch in cyclic multiplication");

    std::vector<int> result = cyclic_convolution(coef, other.coef, q);
    return PolyModQ(result, q);
}

PolyModQ PolyModQ::operator*(int other) const {
    std::vector<int> result(coef.size());
    for (size_t i = 0; i < coef.size(); ++i)
        result[i] = (coef[i] * other) % q;
    return PolyModQ(result, q);
}

PolyModQ operator*(int lhs, const PolyModQ& rhs) {
    return rhs * lhs;
}

ConvPoly PolyModQ::centerLift() const {
    std::vector<int> lifted(coef.size());
    for (size_t i = 0; i < coef.size(); ++i) {
        int x = coef[i];
        // Center  lift: map x in [0, q-1] to [-q/2, q/2)
        if (x > q / 2) x -= q;
        lifted[i] = x;
    }
    return ConvPoly(lifted, static_cast<int>(lifted.size()));
}
/*
int main() {
    // Example polynomials: 2x^2 + 3x + 4 and x^2 + 1, modulus 5
    std::vector<int> coef1 = {4, 3, 2}; // 2x^2 + 3x + 4
    std::vector<int> coef2 = {1, 0, 1}; // x^2 + 1
    int q = 5;

    PolyModQ p1(coef1, q);
    PolyModQ p2(coef2, q);

    std::cout << "p1: " << p1 << std::endl;
    std::cout << "p2: " << p2 << std::endl;

    // Addition
    PolyModQ sum = p1 + p2;
    std::cout << "p1 + p2: " << sum << std::endl;

    // Subtraction
    PolyModQ diff = p1 - p2;
    std::cout << "p1 - p2: " << diff << std::endl;

    // Multiplication
    PolyModQ prod = p1 * p2;
    std::cout << "p1 * p2: " << prod << std::endl;

    // Scalar multiplication
    PolyModQ scalar = p1 * 3;
    std::cout << "p1 * 3: " << scalar << std::endl;

    // Center lift
    ConvPoly lifted = p1.centerLift();
    std::cout << "centerLift(p1): " << lifted << std::endl;

    return 0;
}*/
/*Sample output:
kevin_abraham@kevinabraham:~/NTRU$ g++ -Iinclude src/polymodq.cpp src/convpoly.cpp -o bin/polymodq_test
kevin_abraham@kevinabraham:~/NTRU$ ./bin/polymodq_test
p1: PolyModQ([4, 3, 2], 5)
p2: PolyModQ([1, 0, 1], 5)
p1 + p2: PolyModQ([0, 3, 3], 5)
p1 - p2: PolyModQ([3, 3, 1], 5)
p1 * p2: PolyModQ([4, 3, 1, 3, 2], 5)
p1 * 3: PolyModQ([2, 4, 1], 5)
centerLift(p1): ConvPoly[-1, -2, 2]
kevin_abraham@kevinabraham:~/NTRU$ 
*/