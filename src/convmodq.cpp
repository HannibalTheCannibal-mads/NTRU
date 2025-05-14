#include "ntru/convmodq.hpp"
#include "ntru/polymodq.hpp"
#include "ntru/convpoly.hpp"
#include "ntru/utils.hpp"      // for invModQ
#include <algorithm>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <string>
#include <bitset>
#include <numeric>
#include <cassert>
#include <iterator>
#include <sstream>
#include <iomanip>
#include <openssl/buffer.h>
#include <openssl/evp.h>
#include <openssl/err.h>
#include <openssl/bio.h>
// Default constructor
ConvModQ::ConvModQ() : PolyModQ(), N(0) {}
// At the top of your file, after includes:
std::vector<int> cyclic_convolution(const std::vector<int>& a, const std::vector<int>& b, int q) {
    int N = a.size();
    std::vector<int> result(N, 0);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            int idx = (i + j) % N;
            result[idx] = (result[idx] + a[i] * b[j]) % q;
            if (result[idx] < 0) result[idx] += q; // handle negatives
        }
    }
    return result;
}

// Main constructor
ConvModQ::ConvModQ(const std::vector<int>& coef, int q, int N_)
    : PolyModQ({}, q), N(N_ == -1 ? static_cast<int>(coef.size()) : N_) {
    std::vector<int> temp = coef;
    for (auto& x : temp) x = ((x % q) + q) % q;
    temp.resize(N, 0); // Always pad to N
    this->coef = temp;
    int deg = N-1;
    while (deg > 0 && this->coef[deg] == 0) --deg;
    this->degree = deg;
    this->q = q;
}
// Equality operator

// ostream operator<<
std::ostream& operator<<(std::ostream& os, const ConvModQ& p) {
    os << "ConvModQ([";
    for (size_t i = 0; i < p.coef.size(); ++i) {
        os << p.coef[i];
        if (i != p.coef.size() - 1) os << ", ";
    }
    os << "], " << p.q << ")";
    return os;
}

// Convolution (cyclic) multiplication
ConvModQ ConvModQ::operator*(const ConvModQ& other) const {
    if (N != other.N) throw std::invalid_argument("ConvModQ: N mismatch in multiplication");
    std::vector<int> coefs(N, 0);
    for (int k = 0; k < N; ++k) {
        int s = 0;
        for (int i = 0; i < N; ++i)
            s += this->coef[i] * other.coef[(k - i + N) % N];
        coefs[k] = ((s % q) + q) % q; // Always non-negative
    }
    return ConvModQ(coefs, q, N);
}

ConvModQ ConvModQ::operator*(int other) const {
    std::vector<int> result(this->coef.size());
    for (size_t i = 0; i < this->coef.size(); ++i)
        result[i] = ((this->coef[i] * other) % q + q) % q;
    return ConvModQ(result, q, N);
}


ConvModQ operator*(int lhs, const ConvModQ& rhs) {
    return rhs * lhs;
}

ConvModQ ConvModQ::operator/(const ConvModQ& other) const {
    return *this * other.inverse();
}

ConvModQ ConvModQ::operator/(int other) const {
    int otherinv = invModQ(other, q);
    if (otherinv == -1)
        throw std::runtime_error("Not invertible mod q");
    return *this * otherinv;
}

ConvModQ ConvModQ::modQ(int q_) const {
    return ConvModQ(this->coef, q_, N);
}
ConvModQ ConvModQ::operator*(const ConvPoly& other) const {
    // Convert ConvPoly to ConvModQ with the same modulus and degree
    ConvModQ other_modq(other.coef, this->q, this->N);
    return (*this) * other_modq;
}
ConvPoly ConvModQ::centerLift() const {
    std::vector<int> lifted(coef.size());
    for (size_t i = 0; i < coef.size(); ++i) {
        int x = coef[i];
        // Center lift: map x in [0, q-1] to [-q/2, q/2)
        if (x > q / 2) x -= q;
        lifted[i] = x;
    }
    return ConvPoly(lifted, static_cast<int>(lifted.size()));
}
// Remove leading zeros and update degree
void PolyModQ::trim() {
    while (!coef.empty() && coef.back() == 0) {
        coef.pop_back();
    }
    degree = coef.empty() ? -1 : static_cast<int>(coef.size()) - 1;
}

// Return true if the polynomial is zero
bool PolyModQ::is_zero() const {
    return coef.empty() || std::all_of(coef.begin(), coef.end(), [](int x){return x==0;});
}
// Extended Euclidean Algorithm for division
ConvModQ ConvModQ::inverse(int Nval,bool debug) const {
    int orig_q = this->q;
    int N = Nval;
    if (N == -1) N = this->N; // Use object's N if not provided

    int q = this->q;
    if (q == 2048) q = 2;
    const int FAIL = 100000;
    int i = 0;
    std::vector<ConvModQ> quotients;

    // q(x) = x^N - 1
    std::vector<int> q_poly(N + 1, 0);
    q_poly[0] = -1;
    q_poly[N] = 1;
    ConvModQ qx(q_poly, q, N);

    ConvModQ k(std::vector<int>{0}, q, N);
    ConvModQ b(this->coef, q, N);
    ConvModQ r = qx;
    // If the polynomial is the identity, its inverse is itself
    if (coef.size() == N && coef[0] == 1 && std::all_of(coef.begin() + 1, coef.end(), [](int x){ return x == 0; })) {
        return *this;
    }
        // If the polynomial is a constant (all other coefficients are zero)
    if (coef.size() == N && std::all_of(coef.begin() + 1, coef.end(), [](int x){ return x == 0; })) {
        int c = coef[0];
        int c_inv = invModQ(c, q); // Implement this for modular inverse
        if (c_inv == -1) {
            // Not invertible
            return ConvModQ(std::vector<int>(N, 0), q, N);
        }
        std::vector<int> inv_coef(N, 0);
        inv_coef[0] = c_inv;
        return ConvModQ(inv_coef, q, N);
    }

    int bdinv = invModQ(b.coef[b.degree], q);
    if (bdinv == -1) {
        return ConvModQ(std::vector<int>(N, 0), q, N);; // No inverse
    }

    while (r.degree >= b.degree && i < FAIL) {
        int factor = (r.coef[r.degree] * bdinv) % q;
        if (factor < 0) factor += q;
        std::vector<int> kp_coef(r.degree - b.degree + 1, 0);
        kp_coef.back() = factor;
        ConvModQ kp(kp_coef, q, N);
        k = k + kp;
        r = r - kp * b;
        i++;
    }
    quotients.push_back(k);

    while (!r.is_zero() && i < FAIL) {
        qx = b;
        b = r;
        k = ConvModQ(std::vector<int>(N + 1, 0), q, N);
        r = qx;
        bdinv = invModQ(b.coef[b.degree], q);
        if (bdinv == -1) {
            return ConvModQ(std::vector<int>(N, 0), q, N);;
        }
        while (r.degree >= b.degree && !r.is_zero() && i < FAIL) {
            int factor = (r.coef[r.degree] * bdinv) % q;
            if (factor < 0) factor += q;
            std::vector<int> kp_coef(r.degree - b.degree + 1, 0);
            kp_coef.back() = factor;
            ConvModQ kp(kp_coef, q, N);
            k = k + kp;
            r = r - kp * b;
            i++;
        }
        quotients.push_back(k);
        i++;
    }
    if (i >= FAIL) {
        // Optionally, keep this line for user feedback:
        // std::cout << "Failed to generate inverse in " << FAIL << " steps, stopping." << std::endl;
        return ConvModQ(std::vector<int>(N, 0), q, N);;
    }

    // x and y sequences
    std::vector<ConvModQ> x, y;
    x.push_back(ConvModQ(std::vector<int>{0}, q, N));
    x.push_back(ConvModQ(std::vector<int>{1}, q, N));
    y.push_back(ConvModQ(std::vector<int>{1}, q, N));
    y.push_back(ConvModQ(std::vector<int>{0}, q, N));
    for (size_t index = 0; index < quotients.size(); ++index) {
        x.push_back(quotients[index] * x[index + 1] + x[index]);
        y.push_back(quotients[index] * y[index + 1] + y[index]);
    }

    ConvModQ tinv;
    if (q == 2) {
        int n = 2;
        q = 2048;
        tinv = ConvModQ(x[x.size() - 2].coef, q, N);
        while (n <= 2048) {
            tinv = 2 * tinv - (*this) * tinv * tinv;
            n *= 2;
        }
        return tinv;
    }
    tinv = ConvModQ(x[x.size() - 2].coef, q, N);
    tinv = q * tinv - (*this) * tinv * tinv;
    return 2 * tinv;
}

// Function to generate a random trinary polynomial
// with specified number of 1s and -1s
ConvModQ ConvModQ::operator+(const ConvModQ& other) const {
    std::vector<int> result_coef(N, 0);
    for (int i = 0; i < N; ++i)
        result_coef[i] = (this->coef[i] + other.coef[i]) % q;
    return ConvModQ(result_coef, q, N);
}

ConvModQ ConvModQ::operator-(const ConvModQ& other) const {
    std::vector<int> result_coef(N, 0);
    for (int i = 0; i < N; ++i)
        result_coef[i] = ((this->coef[i] - other.coef[i]) % q + q) % q;
    return ConvModQ(result_coef, q, N);
}
    
// Helper to print a polynomial
void print_poly(const std::vector<int>& coef, const std::string& label) {
    std::cout << label << ": ";
    for (int c : coef) std::cout << c << " ";
    std::cout << std::endl;
}


// Generates a random trinary polynomial with `ones` +1s, `neg_ones` -1s, rest 0
ConvModQ random_trinary_convmodq(int N, int q, int ones, int neg_ones) {
    std::vector<int> arr;
    arr.insert(arr.end(), ones, 1);
    arr.insert(arr.end(), neg_ones, -1);
    arr.insert(arr.end(), N - ones - neg_ones, 0);
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(arr.begin(), arr.end(), g);
    return ConvModQ(arr, q, N);
}

int main() {
    const int N = 11;
    const int q = 41;
    const int p = 3;
    const int df = 2; // 2 +1s, 2 -1s, rest 0
    const int d_g = N / 3; // 3
    const int d_r = df;
    // 2 x +1, 2 x -1, rest 0
    std::vector<int> f_coef = {1, 1, -1, 0, 0, 0, -1, 0, 0, 0, 0}; // length N=11
    ConvModQ f(f_coef, q, N);
    ConvModQ finv = f.inverse(-1); // Should work for a correct inverse implementation
        if (finv.coef.empty()) {
            std::cerr << "Inverse not found!" << std::endl;
            return 1;
        }
    // --- Print f and its inverse ---

    print_poly(f.coef, "Private key f");
    print_poly(finv.coef, "Inverse f^-1 mod q");

    // --- Generate random g ---
    ConvModQ g = random_trinary_convmodq(N, q, d_g, d_g);
    print_poly(g.coef, "Random g");

    // --- Compute public key h = p * f^-1 * g mod q ---
    // (Assume operator* and scalar multiplication are defined)
    // h = (p * finv) * g mod q
    std::vector<int> pfinv_coef(N);
    for (int i = 0; i < N; ++i)
        pfinv_coef[i] = (p * finv.coef[i]) % q;
    ConvModQ pfinv(pfinv_coef, q, N);

    ConvModQ h = pfinv * g;
    print_poly(h.coef, "Public key h");

    // --- Verify inverse ---
    ConvModQ prod = f * finv;
    print_poly(prod.coef, "f * f^-1");

    bool is_identity = (prod.coef.size() > 0 && (prod.coef[0] % q + q) % q == 1);
    for (size_t i = 1; i < prod.coef.size(); ++i)
        if (((prod.coef[i] % q) + q) % q != 0) is_identity = false;
    if (is_identity)
        std::cout << "Inverse verified: f * f^-1 == 1 mod (x^N-1, q)" << std::endl;
    else
        std::cout << "Inverse check failed!" << std::endl;

    // --- Generate ephemeral r for encryption (optional) ---
    ConvModQ r = random_trinary_convmodq(N, q, d_r, d_r);
    print_poly(r.coef, "Ephemeral r");

    return 0;
}