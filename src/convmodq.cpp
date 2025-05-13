#include "ntru/convmodq.hpp"
#include "ntru/polymodq.hpp"
#include "ntru/convpoly.hpp"
#include "ntru/utils.hpp"      // for invModQ
#include <algorithm>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>
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
    std::vector<int> temp(coef.size());
    std::transform(coef.begin(), coef.end(), temp.begin(),
                   [q](int x) { return ((x % q + q) % q); });
    if (N_ == -1) {
        N = static_cast<int>(temp.size());
        this->coef = temp;
    } else {
        N = N_;
        this->coef = temp;
        this->coef.resize(N, 0);
    }
    this->degree = static_cast<int>(this->coef.size()) - 1;
    this->q = q;
}

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

// Convolution multiplication
ConvModQ ConvModQ::operator*(const ConvModQ& other) const {
    if (N != other.N) throw std::invalid_argument("ConvModQ: N mismatch in multiplication");
    std::vector<int> coefs(N, 0);
    for (int k = 0; k < N; ++k) {
        int s = 0;
        for (int i = 0; i < N; ++i)
            s += this->coef[i] * other.coef[(k - i + N) % N];
        coefs[k] = s % q;
    }
    return ConvModQ(coefs, q, N);
}

// TODO: Implement if needed
// ConvModQ ConvModQ::operator*(const ConvPoly& other) const { ... }

ConvModQ ConvModQ::operator*(int other) const {
    std::vector<int> result(this->coef.size());
    for (size_t i = 0; i < this->coef.size(); ++i)
        result[i] = (this->coef[i] * other) % q;
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

ConvModQ ConvModQ::inverse(int N_, bool debug) const {
    int FAIL = 100000;
    int i = 0;
    int Nval = (N_ == -1) ? N : N_;
    int work_q = q;
    int orig_q = q;

    if (debug) std::cout << "[DEBUG] Starting inverse computation\n";

    // If q == 2048, temporarily set q = 2 for the EEA part (as in your Python)
    if (work_q == 2048) {
        if (debug) std::cout << "[DEBUG] work_q changed from 2048 to 2 for EEA\n";
        work_q = 2;
    }

    // Extended Euclidean Algorithm setup
    std::vector<int> qpoly_vec(Nval, 0);
    qpoly_vec[0] = -1;
    qpoly_vec[Nval - 1] = 1;
    PolyModQ qpoly(qpoly_vec, work_q);
    std::vector<int> bvec = this->coef;
    while (bvec.size() < Nval) bvec.push_back(0);
    PolyModQ bpoly(bvec, work_q);

    PolyModQ kpoly(std::vector<int>{0}, work_q);
    while (kpoly.coef.size() < Nval) kpoly.coef.push_back(0);
    PolyModQ rpoly = qpoly;
    std::vector<PolyModQ> quotients;

    if (debug) {
        std::cout << "[DEBUG] Initial qpoly: " << qpoly << "\n";
        std::cout << "[DEBUG] Initial bpoly: " << bpoly << "\n";
    }

    int bdinv = invModQ(bpoly.coef.back(), work_q);
    if (debug) std::cout << "[DEBUG] Initial bdinv: " << bdinv << "\n";
    if (bdinv == -1) {
        if (debug) std::cout << "[DEBUG] bdinv == -1, returning empty ConvModQ\n";
        return ConvModQ({}, orig_q, Nval);
    }

    while (rpoly.degree >= bpoly.degree && i < FAIL) {
        if (debug) {
            std::cout << "[DEBUG] Top of outer loop, i = " << i << "\n";
            std::cout << "[DEBUG] rpoly: " << rpoly << "\n";
            std::cout << "[DEBUG] bpoly: " << bpoly << "\n";
        }
        std::vector<int>& rcoef = rpoly.coef;
        std::vector<int> kpvec(rpoly.degree - bpoly.degree, 0);
        kpvec.push_back(rcoef[rpoly.degree] * bdinv % work_q);
        PolyModQ kp(kpvec, work_q);
        while (kp.coef.size() < Nval) kp.coef.push_back(0);
        if (debug) std::cout << "[DEBUG] kp: " << kp << "\n";
        kpoly = kpoly + kp;
        while (kpoly.coef.size() < Nval) kpoly.coef.push_back(0);
        rpoly = rpoly - kp * bpoly;
        while (rpoly.coef.size() < Nval) rpoly.coef.push_back(0);
        if (debug) std::cout << "[DEBUG] Updated rpoly: " << rpoly << "\n";
        i++;
    }
    quotients.push_back(kpoly);

    if (debug) {
        std::cout << "[DEBUG] After first quotient: " << kpoly << "\n";
        std::cout << "[DEBUG] qpoly = " << qpoly << " = " << bpoly << "*" << kpoly << " + " << rpoly << "\n";
    }

    while (!(rpoly == PolyModQ()) && i < FAIL) {
        if (debug) std::cout << "[DEBUG] Entering main EEA loop, i = " << i << "\n";
        qpoly = bpoly;
        while (qpoly.coef.size() < Nval) qpoly.coef.push_back(0);
        bpoly = rpoly;
        while (bpoly.coef.size() < Nval) bpoly.coef.push_back(0);
        kpoly = PolyModQ(std::vector<int>(Nval , 0), work_q);
        while (kpoly.coef.size() < Nval) kpoly.coef.push_back(0);
        rpoly = qpoly;
        while (rpoly.coef.size() < Nval) rpoly.coef.push_back(0);

        bdinv = invModQ(bpoly.coef.back(), work_q);
        if (debug) std::cout << "[DEBUG] bdinv in loop: " << bdinv << "\n";
        if (bdinv == -1) {
            if (debug) std::cout << "[DEBUG] bdinv == -1 in loop, returning empty ConvModQ\n";
            return ConvModQ({}, orig_q, Nval);
        }
        std::cout << "bpoly.size=" << bpoly.coef.size() << ", kp.size=" << kpoly.coef.size() << std::endl;
        while (rpoly.degree >= bpoly.degree && !(rpoly == PolyModQ()) && i < FAIL) {
            if (debug) {
                std::cout << "[DEBUG] Inner loop, i = " << i << "\n";
                std::cout << "[DEBUG] rpoly: " << rpoly << "\n";
                std::cout << "[DEBUG] bpoly: " << bpoly << "\n";
            }
            std::vector<int>& rcoef = rpoly.coef;
            std::vector<int> kpvec(rpoly.degree - bpoly.degree, 0);
            kpvec.push_back(rcoef[rpoly.degree] * bdinv % work_q);
            PolyModQ kp(kpvec, work_q);
            // --- Pad kp to Nval ---
            while (kp.coef.size() < Nval) kp.coef.push_back(0);

            if (debug) std::cout << "[DEBUG] kp: " << kp << "\n";
            kpoly = kpoly + kp;
            while (kpoly.coef.size() < Nval) kpoly.coef.push_back(0);

            rpoly = rpoly - kp * bpoly;
            while (rpoly.coef.size() < Nval) rpoly.coef.push_back(0);

            if (debug) std::cout << "[DEBUG] Updated rpoly: " << rpoly << "\n";
            i++;

        }
        std::cout << "bpoly.size=" << bpoly.coef.size() << ", kp.size=" << kpoly.coef.size() << std::endl;
        quotients.push_back(kpoly);
        i++;
        if (debug) {
            std::cout << "[DEBUG] After quotient: " << kpoly << "\n";
            std::cout << "[DEBUG] qpoly = " << qpoly << " = " << bpoly << "*" << kpoly << " + " << rpoly << "\n";
        }
    }
    if (i >= FAIL) {
        std::cerr << "[DEBUG] Failed to generate inverse in " << FAIL << " steps, stopping.\n";
        return ConvModQ({}, orig_q, Nval);
    }
    std::cout << "bpoly.size=" << bpoly.coef.size() << ", kp.size=" << kpoly.coef.size() << std::endl;
    // Back substitution as in Python
    if (debug) std::cout << "[DEBUG] Starting back substitution\n";
    std::vector<PolyModQ> x{PolyModQ({0}, work_q), PolyModQ({1}, work_q)};
    std::vector<PolyModQ> y{PolyModQ({1}, work_q), PolyModQ({0}, work_q)};
    for (size_t index = 0; index < quotients.size(); ++index) {
        x.push_back(quotients[index] * x[index + 1] + x[index]);
        y.push_back(quotients[index] * y[index + 1] + y[index]);
        while (x.back().coef.size() < Nval) x.back().coef.push_back(0);
        while (y.back().coef.size() < Nval) y.back().coef.push_back(0);
        if (debug) {
            std::cout << "[DEBUG] Backsub index " << index << ": x = " << x.back() << ", y = " << y.back() << "\n";
        }
    }

    // ... (your EEA and backsubstitution code above)
    std::cout << "bpoly.size=" << bpoly.coef.size() << ", kp.size=" << kpoly.coef.size() << std::endl;
    // --- Normalization step ---
    // The GCD is in bpoly after the last EEA loop (since rpoly == 0)
    bpoly.trim();
    if (bpoly.degree != 0 || bpoly.coef.empty() || bpoly.coef[0] == 0) {
        if (debug) std::cout << "[DEBUG] GCD is not a unit, no inverse exists.\n";
        return ConvModQ({}, orig_q, Nval);
    }
    int gcd_const = bpoly.coef[0];
    int gcd_inv = invModQ(gcd_const, work_q);
    if (gcd_inv == -1) {
        if (debug) std::cout << "[DEBUG] GCD constant not invertible, no inverse exists.\n";
        return ConvModQ({}, orig_q, Nval);
    }

    // Normalize the inverse from EEA
    PolyModQ inverse_poly = x[x.size() - 2];
    for (auto& c : inverse_poly.coef) {
        c = (c * gcd_inv) % work_q;
        if (c < 0) c += work_q;
    }
    inverse_poly.trim();
    while (inverse_poly.coef.size() < Nval) inverse_poly.coef.push_back(0);
    // --- Newton iteration for q == 2048 (special NTRU case) ---
    if (work_q == 2) {
        int n = 2;
        ConvModQ tinv(inverse_poly.coef, orig_q, Nval);
        if (debug) std::cout << "[DEBUG] Newton iteration for q=2048\n";
        while (n <= 2048) {
            tinv = ConvModQ((2 * tinv - (*this) * tinv * tinv).coef, orig_q, Nval);
            n *= 2;
            if (debug) std::cout << "[DEBUG] Newton step n=" << n << ", tinv=" << tinv << "\n";
        }
        // Final normalization for Newton result
        ConvModQ product = (*this) * tinv;
        int c = product.coef[0];
        int c_inv = invModQ(c, orig_q);
        ConvModQ tinv_true = tinv * c_inv;
        tinv_true.trim();
        if (debug) std::cout << "[DEBUG] Returning Newton-normalized inverse: " << tinv_true << "\n";
        return tinv_true;
    }
    std::cout << "bpoly.size=" << bpoly.coef.size() << ", kp.size=" << kpoly.coef.size() << std::endl;
    // --- Final normalization for general case ---
    ConvModQ f_inv(inverse_poly.coef, orig_q, Nval);
    // Now ensure f * f_inv = [1, 0, ...] (not just [c, 0, ...])
    ConvModQ product = (*this) * f_inv;
    int c = product.coef[0];
    int c_inv = invModQ(c, orig_q);
    ConvModQ f_inv_true = f_inv * c_inv;
    f_inv_true.trim();
    while (f_inv_true.coef.size() < Nval) f_inv_true.coef.push_back(0);
    if (debug) std::cout << "[DEBUG] Returning fully normalized inverse: " << f_inv_true << "\n";
    return f_inv_true;


}


// Helper to print a polynomial
void print_poly(const std::vector<int>& coef, const std::string& label) {
    std::cout << label << ": ";
    for (int c : coef) std::cout << c << " ";
    std::cout << std::endl;
}

// Generate a random trinary polynomial of degree N, mod q
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
    int N = 11;      // Small N for demonstration; use your real N for real keys
    int q = 41;      // Small q for demonstration; use your real q for real keys
    ConvModQ f, finv;
    do {
        f = random_trinary_convmodq(N, q, 4, 3);
        finv = f.inverse(-1, true);
    } while ( finv.coef.empty() || std::all_of(finv.coef.begin(), finv.coef.end(), [](int x){ return x == 0; }));

    print_poly(finv.coef, "f^-1");

    // Multiply and check if the result is the identity
    ConvModQ prod = f * finv;
    print_poly(prod.coef, "f * f^-1");

    // Check if prod is [1, 0, 0, ..., 0] mod q
    bool is_identity = (prod.coef.size() > 0 && (prod.coef[0] % q + q) % q == 1);
    for (size_t i = 1; i < prod.coef.size(); ++i)
        if (((prod.coef[i] % q) + q) % q != 0) is_identity = false;

    if (is_identity)
        std::cout << "Inverse verified: f * f^-1 == 1 mod (x^N-1, q)" << std::endl;
    else
        std::cout << "Inverse check failed!" << std::endl;

    return 0;
}
