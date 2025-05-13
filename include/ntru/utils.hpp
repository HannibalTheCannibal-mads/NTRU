#ifndef NTRU_UTILS_HPP
#define NTRU_UTILS_HPP

#include <vector>
#include <string>
#include <sstream>
#include <cstddef>
#include <cmath>

// Modular inverse
int invModQ(int d, int q);

// Serialize coefficients to binary blob
std::vector<unsigned char> coeffs_to_blob(const std::vector<int>& coeffs);

// PEM encoding helpers
std::string to_pem(const std::vector<int>& coeffs, const std::string& header = "NTRU PUBLIC KEY", int bytes_per_coeff = 2);
std::string wrap_lines(const std::string& input, size_t width = 64);

// Base64 encoding
std::string base64_encode(const std::vector<unsigned char>& data);

// Chunk a vector into pieces of size N
template <typename T>
std::vector<std::vector<T>> chunk(size_t N, const std::vector<T>& iter);

// String/bit conversion helpers
std::string strToBin(const std::string& m);
std::string binToStr(const std::vector<int>& b);

// Template implementations must be in the header:
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

#endif // NTRU_UTILS_HPP
// This code provides utility functions for modular arithmetic, serialization, PEM encoding,
// base64 encoding, and chunking vectors. It also includes string conversion functions for binary