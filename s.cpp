#include <iostream>        // cout, endl
#include <vector>          // vector
#include <algorithm>       // shuffle, transform, max, min
#include <initializer_list>// initializer_list
#include <string>          // string, ostringstream
#include <bitset>          // bitset
#include <iterator>        // back_inserter, etc.
#include <cmath>           // ceil, abs
#include <random>          // random_device, mt19937
#include <numeric>         // accumulate
#include <sstream>         // ostringstream
#include <iomanip>         // setw, setprecision
#include <openssl/evp.h>   // For base64 encoding
#include <chrono>          // timing
#include <stdexcept>       // exceptions
#include <cassert>
#include <openssl/buffer.h>



int invModQ(int d, int q) {
    for (int i = 1; i < q; ++i) {
        if ((d * i) % q == 1) {
            return i;
        }
    }
    return -1; // Use -1 to indicate "not found"
}

class ConvPoly {
public:
    std::vector<int> coef;
    int N;

    // Constructors
    ConvPoly(std::vector<int> coef_ = {0}, int N_ = -1) {
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
    friend std::ostream& operator<<(std::ostream& os, const ConvPoly& p) {
        os << "ConvPoly[";
        for (size_t i = 0; i < p.coef.size(); ++i) {
            os << p.coef[i];
            if (i != p.coef.size() - 1) os << ", ";
        }
        os << "]";
        return os;
    }

    // Addition
    ConvPoly operator+(const ConvPoly& other) const {
        if (N == other.N) {
            std::vector<int> result(N);
            for (int i = 0; i < N; ++i)
                result[i] = coef[i] + other.coef[i];
            return ConvPoly(result);
        }
        throw std::invalid_argument("ConvPoly: N mismatch in addition");
    }
    ConvPoly operator+(int other) const {
        std::vector<int> result = coef;
        result[0] += other;
        return ConvPoly(result, N);
    }
    friend ConvPoly operator+(int lhs, const ConvPoly& rhs) {
        return rhs + lhs;
    }

    // Equality
    bool operator==(const ConvPoly& other) const {
        return coef == other.coef;
    }
    bool operator!=(const ConvPoly& other) const {
        return !(*this == other);
    }

    // Negation
    ConvPoly operator-() const {
        std::vector<int> result(N);
        for (int i = 0; i < N; ++i)
            result[i] = -coef[i];
        return ConvPoly(result, N);
    }

    // Subtraction
    ConvPoly operator-(const ConvPoly& other) const {
        return *this + (-other);
    }
    ConvPoly operator-(int other) const {
        return *this + (-other);
    }
    friend ConvPoly operator-(int lhs, const ConvPoly& rhs) {
        return (-rhs) + lhs;
    }

    // Multiplication (convolution)
    ConvPoly operator*(const ConvPoly& other) const {
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
    ConvPoly operator*(int other) const {
        std::vector<int> result(N);
        for (int i = 0; i < N; ++i)
            result[i] = coef[i] * other;
        return ConvPoly(result, N);
    }
    friend ConvPoly operator*(int lhs, const ConvPoly& rhs) {
        return rhs * lhs;
    }
};

class PolyModQ {
public:
    std::vector<int> coef;
    int degree;
    int q;

    PolyModQ(const std::vector<int>& coef_ = {0}, int q_ = 3) : q(q_) {
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

    bool operator==(const PolyModQ& other) const {
        return degree == other.degree && coef == other.coef;
    }
    friend std::ostream& operator<<(std::ostream& os, const PolyModQ& p) {
    os << "PolyModQ([";
    for (size_t i = 0; i < p.coef.size(); ++i) {
        os << p.coef[i];
        if (i != p.coef.size() - 1) os << ", ";
    }
    os << "], " << p.q << ")";
    return os;
    }


    bool operator!=(const PolyModQ& other) const {
        return !(*this == other);
    }

    PolyModQ operator-() const {
        std::vector<int> neg(coef.size());
        for (size_t i = 0; i < coef.size(); ++i)
            neg[i] = ((-coef[i] % q) + q) % q;
        return PolyModQ(neg, q);
    }

    PolyModQ operator+(const PolyModQ& other) const {
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

    PolyModQ operator+(int other) const {
        std::vector<int> result = coef;
        if (!result.empty())
            result[0] = (result[0] + other) % q;
        return PolyModQ(result, q);
    }

    friend PolyModQ operator+(int lhs, const PolyModQ& rhs) {
        return rhs + lhs;
    }

    PolyModQ operator-(const PolyModQ& other) const {
        return *this + (-other);
    }

    PolyModQ operator-(int other) const {
        return *this + (-other);
    }

    friend PolyModQ operator-(int lhs, const PolyModQ& rhs) {
        return (-rhs) + lhs;
    }

    PolyModQ operator*(const PolyModQ& other) const {
        if (q != other.q)
            throw std::invalid_argument("PolyModQ: modulus q mismatch in multiplication");

        std::vector<int> result(degree + other.degree + 1, 0);
        for (size_t i = 0; i < coef.size(); ++i) {
            for (size_t j = 0; j < other.coef.size(); ++j) {
                result[i + j] = (result[i + j] + coef[i] * other.coef[j]) % q;
            }
        }
        return PolyModQ(result, q);
    }

    PolyModQ operator*(int other) const {
        std::vector<int> result(coef.size());
        for (size_t i = 0; i < coef.size(); ++i)
            result[i] = (coef[i] * other) % q;
        return PolyModQ(result, q);
    }

    friend PolyModQ operator*(int lhs, const PolyModQ& rhs) {
        return rhs * lhs;
    }

    ConvPoly centerLift() const;
};


class ConvModQ : public PolyModQ {
public:
    int N;
    ConvModQ() : PolyModQ(), N(0) {}
    ConvModQ(const std::vector<int>& coef, int q = 3, int N_ = -1)
        : PolyModQ({}, q), N(N_ == -1 ? static_cast<int>(coef.size()) : N_) {
        std::vector<int> temp(coef.size());
        std::transform(coef.begin(), coef.end(), temp.begin(),
                       [this](int x) { return ((x % q + q) % q); });
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

    friend std::ostream& operator<<(std::ostream& os, const ConvModQ& p) {
        os << "ConvModQ([";
        for (size_t i = 0; i < p.coef.size(); ++i) {
            os << p.coef[i];
            if (i != p.coef.size() - 1) os << ", ";
        }
        os << "], " << p.q << ")";
        return os;
    }

    // Convolution multiplication
    ConvModQ operator*(const ConvModQ& other) const {
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

    ConvModQ operator*(const ConvPoly& other) const; // You can implement this if needed

    ConvModQ operator*(int other) const {
        std::vector<int> result(this->coef.size());
        for (size_t i = 0; i < this->coef.size(); ++i)
            result[i] = (this->coef[i] * other) % q;
        return ConvModQ(result, q, N);
    }

    friend ConvModQ operator*(int lhs, const ConvModQ& rhs) {
        return rhs * lhs;
    }

    // Division by polynomial or integer
    ConvModQ operator/(const ConvModQ& other) const {
        return *this * other.inverse();
    }

    ConvModQ operator/(int other) const {
        int otherinv = invModQ(other, q);
        if (otherinv == -1)
            throw std::runtime_error("Not invertible mod q");
        return *this * otherinv;
    }

    // modQ (change modulus)
    ConvModQ modQ(int q_) const {
        return ConvModQ(this->coef, q_, N);
    }

    // Inverse using extended Euclidean algorithm
    ConvModQ inverse(int N_ = -1, bool debug = false) const {
        int FAIL = 100000;
        int i = 0;
        int Nval = (N_ == -1) ? N : N_;
        std::vector<PolyModQ> quotients;

        PolyModQ qpoly(std::vector<int>(Nval + 1, 0), q);
        qpoly.coef[0] = -1;
        qpoly.coef[Nval] = 1;
        qpoly.degree = Nval;

        PolyModQ kpoly(std::vector<int>{0}, q);
        PolyModQ bpoly(this->coef, q);
        PolyModQ rpoly = qpoly;

        int bdinv = invModQ(bpoly.coef.back(), q);
        if (bdinv == -1) return ConvModQ({}, q, Nval);

        while (rpoly.degree >= bpoly.degree && i < FAIL) {
            std::vector<int>& rcoef = rpoly.coef;
            std::vector<int> kpvec(rpoly.degree - bpoly.degree + 1, 0);
            kpvec.back() = rcoef[rpoly.degree] * bdinv % q;
            PolyModQ kp(kpvec, q);
            kpoly = kpoly + kp;
            rpoly = rpoly - kp * bpoly;
            i++;
        }
        quotients.push_back(kpoly);

        if (debug) {
            std::cout << qpoly << " = " << bpoly << "*" << kpoly << " + " << rpoly << std::endl;
        }

        while (!(rpoly == PolyModQ()) && i < FAIL) {
            qpoly = bpoly;
            bpoly = rpoly;
            kpoly = PolyModQ(std::vector<int>(Nval + 1, 0), q);
            rpoly = qpoly;

            bdinv = invModQ(bpoly.coef.back(), q);
            if (bdinv == -1) return ConvModQ({}, q, Nval);

            while (rpoly.degree >= bpoly.degree && !(rpoly == PolyModQ()) && i < FAIL) {
                std::vector<int>& rcoef = rpoly.coef;
                std::vector<int> kpvec(rpoly.degree - bpoly.degree + 1, 0);
                kpvec.back() = rcoef[rpoly.degree] * bdinv % q;
                PolyModQ kp(kpvec, q);
                kpoly = kpoly + kp;
                rpoly = rpoly - kp * bpoly;
                i++;
            }
            quotients.push_back(kpoly);
            i++;
            if (debug) {
                std::cout << qpoly << " = " << bpoly << "*" << kpoly << " + " << rpoly << std::endl;
            }
        }
        if (i >= FAIL) {
            std::cerr << "Failed to generate inverse in " << FAIL << " steps, stopping." << std::endl;
            return ConvModQ({}, q, Nval);
        }
        std::vector<PolyModQ> x{PolyModQ({0}, q), PolyModQ({1}, q)};
        std::vector<PolyModQ> y{PolyModQ({1}, q), PolyModQ({0}, q)};
        for (size_t index = 0; index < quotients.size(); ++index) {
            x.push_back(quotients[index] * x[index + 1] + x[index]);
            y.push_back(quotients[index] * y[index + 1] + y[index]);
        }
        if (q == 2) {
            int n = 2;
            int oldq = q;
            const_cast<ConvModQ*>(this)->q = 2048;
            ConvModQ tinv(x[x.size() - 2].coef, this->q, Nval);
            while (n <= 2048) {
                tinv = ConvModQ((2 * tinv - (*this) * tinv * tinv).coef, this->q, this->N);
                n *= 2;
            }
            const_cast<ConvModQ*>(this)->q = oldq;
            return tinv;
        }
        ConvModQ tinv(x[x.size() - 2].coef, q, Nval);
        tinv = q * tinv - (*this) * tinv * tinv;
        return 2 * tinv;
    }
};



class NTRUParams {
public:
    int N, d_f, d_g, d_r, q, p;

    NTRUParams(int k, const std::string& choice = "speed") {
        if (k == 112) {
            if (choice == "space") {
                N = 401; d_f = 113;
            } else if (choice == "hybrid") {
                N = 541; d_f = 49;
            } else if (choice == "speed") {
                N = 659; d_f = 38;
            }
        } else if (k == 128) {
            if (choice == "space") {
                N = 449; d_f = 134;
            } else if (choice == "hybrid") {
                N = 613; d_f = 55;
            } else if (choice == "speed") {
                N = 761; d_f = 42;
            }
        } else if (k == 192) {
            if (choice == "space") {
                N = 677; d_f = 157;
            } else if (choice == "hybrid") {
                N = 887; d_f = 81;
            } else if (choice == "speed") {
                N = 1087; d_f = 63;
            }
        } else if (k == 256) {
            if (choice == "space") {
                N = 1087; d_f = 120;
            } else if (choice == "hybrid") {
                N = 1171; d_f = 106;
            } else if (choice == "speed") {
                N = 1499; d_f = 79;
            }
        } else {
            throw std::runtime_error("Not implemented. :(");
        }
        q = 2048;
        p = 3;
        d_g = N / 3;
        d_r = d_f;
    }
};



class NTRUKey {
public:
    NTRUParams ring;
    ConvModQ f, g, finvq, finvp, h;

    NTRUKey(const NTRUParams& ring_param = NTRUParams(256, "speed"),
            const ConvModQ* f_ = nullptr,
            const ConvModQ* g_ = nullptr)
        : ring(ring_param)
    {
        // Generate f if not provided
        if (f_ == nullptr) {
            f = randomTrinary(ring.d_f + 1, ring.d_f);
        } else {
            f = *f_;
        }

        // Generate g if not provided
        if (g_ == nullptr) {
            g = randomTrinary(ring.d_g, ring.d_g);
        } else {
            g = *g_;
        }

        // Compute inverses
        finvq = f.inverse();
        while (finvq.coef.empty()) {
            std::cout << "finv was None. Retrying." << std::endl;
            f = randomTrinary(ring.d_f + 1, ring.d_f);
            finvq = f.inverse();
        }

        finvp = ConvModQ(f.centerLift().coef, ring.p, ring.N).inverse();
        while (finvq.coef.empty() || finvp.coef.empty()) {
            std::cout << "finv was None. Retrying." << std::endl;
            f = randomTrinary(ring.d_f + 1, ring.d_f);
            finvq = f.inverse();
            if (finvq.coef.empty()) continue;
            finvp = ConvModQ(f.coef, ring.p, ring.N).inverse();
        }

        h = finvq * g;
    }

    ConvModQ randomTrinary(int d1, int d2) {
        std::vector<int> arr;
        arr.insert(arr.end(), d1, 1);
        arr.insert(arr.end(), d2, -1);
        arr.insert(arr.end(), ring.N - d1 - d2, 0);

        // Shuffle using random device
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(arr.begin(), arr.end(), g);

        return ConvModQ(arr, ring.q, ring.N);
    }

    std::pair<NTRUParams, ConvModQ> publicKey() const {
        return std::make_pair(ring, h);
    }
};




std::string strToBin(const std::string& m) {
    std::string result;
    for (char c : m) {
        std::bitset<8> bits(static_cast<unsigned char>(c));
        result += bits.to_string();
    }
    return result;
}




std::string binToStr(const std::vector<int>& b) {
    std::string result;
    for (size_t i = 0; i + 8 <= b.size(); i += 8) {
        std::string byteStr;
        for (size_t j = 0; j < 8; ++j)
            byteStr += (b[i + j] ? '1' : '0');
        char c = static_cast<char>(std::stoi(byteStr, nullptr, 2));
        result += c;
    }
    return result;
}


template <typename T>
std::vector<std::vector<T>> chunk(size_t N, const std::vector<T>& iter) {
    std::vector<std::vector<T>> chunks;
    size_t total = iter.size();
    size_t numChunks = static_cast<size_t>(std::ceil(static_cast<double>(total) / N));
    for (size_t i = 0; i < numChunks; ++i) {
        size_t start = i * N;
        size_t end = std::min(start + N, total);
        chunks.emplace_back(iter.begin() + start, iter.begin() + end);
    }
    return chunks;
}






std::pair<std::vector<ConvModQ>, int> NTRUEncrypt(const NTRUParams& ring, const ConvModQ& pub, const std::string& m) {
    std::string bin = strToBin(m);
    std::vector<int> mvec;
    for (char c : bin) mvec.push_back(c - '0');
    std::vector<ConvPoly> m_chunks;
    int n;
    if (mvec.size() > static_cast<size_t>(ring.N)) {
        auto msplit = chunk(ring.N, mvec);
        n = static_cast<int>(msplit.back().size());
        for (const auto& chunk : msplit)
            m_chunks.emplace_back(chunk, ring.N);
    } else {
        n = static_cast<int>(mvec.size());
        m_chunks.emplace_back(mvec, ring.N);
    }
    std::vector<ConvModQ> menc;
    for (const auto& mc : m_chunks)
        menc.push_back(NTRUBlockEncrypt(ring, pub, mc));
    return std::make_pair(menc, n);
}


ConvModQ NTRUBlockEncrypt(const NTRUParams& ring, const ConvModQ& h, const ConvPoly& m) {
    std::vector<int> rvars;
    rvars.insert(rvars.end(), ring.d_r, 1);
    rvars.insert(rvars.end(), ring.d_r, -1);
    rvars.insert(rvars.end(), ring.N - 2 * ring.d_r, 0);

    // Shuffle
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(rvars.begin(), rvars.end(), g);

    ConvModQ r(rvars, ring.q, ring.N);
    ConvModQ c = ring.p * r * h + m;
    return c;
}

ConvPoly NTRUBlockDecrypt(const NTRUKey& key, const ConvModQ& c) {
    ConvModQ a = key.f * c;
    ConvPoly aprime = a.centerLift();
    ConvModQ m = key.finvp * aprime;
    return m.centerLift();
}


std::string NTRUDecrypt(const NTRUKey& key, const std::vector<ConvModQ>& cl, int n) {
    if (cl.size() == 1) {
        ConvPoly m = NTRUBlockDecrypt(key, cl[0]);
        std::vector<int> mcoef(m.coef.begin(), m.coef.begin() + n);
        return binToStr(mcoef);
    }
    std::vector<std::vector<int>> mlist;
    for (size_t i = 0; i + 1 < cl.size(); ++i)
        mlist.push_back(NTRUBlockDecrypt(key, cl[i]).coef);
    std::vector<int> last = NTRUBlockDecrypt(key, cl.back()).coef;
    last.resize(n);
    mlist.push_back(last);

    // Flatten mlist
    std::vector<int> m;
    for (const auto& v : mlist)
        m.insert(m.end(), v.begin(), v.end());

    return binToStr(m);
}


// Helper: Base64 encoding using OpenSSL
std::string base64_encode(const std::vector<unsigned char>& data) {
    BIO *bio, *b64;
    BUF_MEM *bufferPtr;
    b64 = BIO_new(BIO_f_base64());
    bio = BIO_new(BIO_s_mem());
    bio = BIO_push(b64, bio);

    // No newlines - PEM requires them, but we'll split manually
    BIO_set_flags(bio, BIO_FLAGS_BASE64_NO_NL);
    BIO_write(bio, data.data(), data.size());
    BIO_flush(bio);
    BIO_get_mem_ptr(bio, &bufferPtr);

    std::string result(bufferPtr->data, bufferPtr->length);
    BIO_free_all(bio);
    return result;
}

// Helper: Split string into lines of given width
std::string wrap_lines(const std::string& input, size_t width = 64) {
    std::ostringstream oss;
    for (size_t i = 0; i < input.size(); i += width) {
        oss << input.substr(i, width) << '\n';
    }
    return oss.str();
}

// Main function
std::string to_pem(const std::vector<int>& coeffs, const std::string& header = "NTRU PUBLIC KEY", int bytes_per_coeff = 2) {
    // Serialize coefficients
    std::vector<unsigned char> binary_blob;
    for (int c : coeffs) {
        for (int i = bytes_per_coeff - 1; i >= 0; --i) {
            binary_blob.push_back(static_cast<unsigned char>((c >> (8 * i)) & 0xFF));
        }
    }

    // Base64 encode
    std::string b64_encoded = base64_encode(binary_blob);

    // Split into lines of 64 characters (PEM convention)
    std::string b64_lines = wrap_lines(b64_encoded, 64);

    // Assemble PEM block
    std::ostringstream pem;
    pem << "-----BEGIN " << header << "-----\n"
        << b64_lines
        << "-----END " << header << "-----";
    return pem.str();
}




// Helper: serialize coefficients to binary blob (signed 2 bytes per coeff)
std::vector<unsigned char> coeffs_to_blob(const std::vector<int>& coeffs) {
    std::vector<unsigned char> blob;
    for (int c : coeffs) {
        // Signed 2-byte big-endian
        blob.push_back((c >> 8) & 0xFF);
        blob.push_back(c & 0xFF);
    }
    return blob;
}

int main() {
    std::cout << "Generating key" << std::endl;
    NTRUKey key(NTRUParams(256, "speed"));

    // Print PEM-encoded public and private key components
    std::cout << "PEM-encoded public key:\n";
    std::cout << to_pem(key.h.coef, "NTRU PUBLIC KEY") << std::endl;

    std::cout << "\nPEM-encoded private key f:\n";
    std::cout << to_pem(key.f.coef, "NTRU PRIVATE KEY COMPONENT_F") << std::endl;

    std::cout << "\nPEM-encoded private key g:\n";
    std::cout << to_pem(key.g.coef, "NTRU PRIVATE KEY COMPONENT_G") << std::endl;

    std::string morig = "2QhSx0PlIU1qNGbpz+76J7rEtVRbb8CahNeymYrbQns=";
    std::cout << "\nEncrypting message" << std::endl;

    // Encrypt once to show details
    auto enc = NTRUEncrypt(key.ring, key.h, morig);
    std::cout << "\nOriginal message: " << morig << std::endl;

    // Serialize ciphertext coefficients into a binary blob (first ciphertext block)
    std::vector<unsigned char> ciphertext_blob = coeffs_to_blob(enc.first[0].coef);

    // Base64 encode the binary blob
    std::string ciphertext_base64 = base64_encode(ciphertext_blob);

    std::cout << "Ciphertext (Base64 encoded): " << ciphertext_base64 << std::endl;

    // Decrypt to verify
    std::string m = NTRUDecrypt(key, enc.first, enc.second);
    std::cout << "Decrypted message: " << m << std::endl;

    // Timing as before
    
    int N = 100;
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i) {
        auto enc = NTRUEncrypt(key.ring, key.h, morig);
        std::string m = NTRUDecrypt(key, enc.first, enc.second);
    }
    auto end = std::chrono::high_resolution_clock::now();
    double total = std::chrono::duration<double>(end - start).count();
    double timeper = total / N;
    std::cout << "\nTotal time for " << N << " cycles: " << std::fixed << std::setprecision(6) << total << " seconds" << std::endl;
    std::cout << "Time per Encryption/Decryption cycle: " << timeper << " seconds" << std::endl;
     assert(m == morig);
    
    return 0;
}

