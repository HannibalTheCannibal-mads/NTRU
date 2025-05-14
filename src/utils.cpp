#include "ntru/utils.hpp"
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <bitset>
#include <algorithm>
#include <openssl/bio.h>
#include <openssl/evp.h>
#include <openssl/buffer.h>

// Modular inverse
int invModQ(int a, int q) {
    int t = 0, newt = 1;
    int r = q, newr = a;
    while (newr != 0) {
        int quotient = r / newr;
        int temp = t;
        t = newt;
        newt = temp - quotient * newt;
        temp = r;
        r = newr;
        newr = temp - quotient * newr;
    }
    if (r > 1) return -1; // No inverse
    if (t < 0) t += q;
    return t;
}
// Modular inverse for negative numbers
int invModQ_neg(int a, int q) {
    int t = 0, newt = 1;
    int r = q, newr = a;
    while (newr != 0) {
        int quotient = r / newr;
        int temp = t;
        t = newt;
        newt = temp - quotient * newt;
        temp = r;
        r = newr;
        newr = temp - quotient * newr;
    }
    if (r > 1) return -1; // No inverse
    if (t < 0) t += q;
    return t % q;
}
// Modular inverse for positive numbers
int invModQ_pos(int a, int q) {
    int t = 0, newt = 1;
    int r = q, newr = a;
    while (newr != 0) {
        int quotient = r / newr;
        int temp = t;
        t = newt;
        newt = temp - quotient * newt;
        temp = r;
        r = newr;
        newr = temp - quotient * newr;
    }
    if (r > 1) return -1; // No inverse
    if (t < 0) t += q;
    return t % q;
}
// Serialize coefficients to binary blob (signed 2 bytes per coeff)
std::vector<unsigned char> coeffs_to_blob(const std::vector<int>& coeffs) {
    std::vector<unsigned char> blob;
    for (int c : coeffs) {
        // Signed 2-byte big-endian
        blob.push_back((c >> 8) & 0xFF);
        blob.push_back(c & 0xFF);
    }
    return blob;
}

// PEM encoding helpers
std::string to_pem(const std::vector<int>& coeffs, const std::string& header, int bytes_per_coeff) {
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

// Split string into lines of given width
std::string wrap_lines(const std::string& input, size_t width) {
    std::ostringstream oss;
    for (size_t i = 0; i < input.size(); i += width) {
        oss << input.substr(i, width) << '\n';
    }
    return oss.str();
}

// Base64 encoding using OpenSSL
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

// Template function chunk<T> is implemented in the header (utils.hpp)

// String to binary string
std::string strToBin(const std::string& m) {
    std::string result;
    for (char c : m) {
        std::bitset<8> bits(static_cast<unsigned char>(c));
        result += bits.to_string();
    }
    return result;
}

// Binary vector to string
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
/*
int main() {
    // Test invModQ
    int d = 7, q = 26;
    int inv = invModQ(d, q);
    std::cout << "invModQ(" << d << ", " << q << ") = " << inv << std::endl;

    // Test coeffs_to_blob
    std::vector<int> coeffs = {300, -21, 42, 1023};
    std::vector<unsigned char> blob = coeffs_to_blob(coeffs);
    std::cout << "coeffs_to_blob: ";
    for (unsigned char c : blob)
        std::cout << std::hex << (int)c << " ";
    std::cout << std::dec << std::endl;

    // Test to_pem
    std::string pem = to_pem(coeffs, "NTRU PUBLIC KEY", 2);
    std::cout << "PEM encoding:\n" << pem << std::endl;

    // Test strToBin and binToStr
    std::string message = "Hi";
    std::string bin = strToBin(message);
    std::cout << "strToBin(\"" << message << "\") = " << bin << std::endl;

    // Convert binary string to vector<int> for binToStr
    std::vector<int> binvec;
    for (char ch : bin) binvec.push_back(ch - '0');
    std::string recovered = binToStr(binvec);
    std::cout << "binToStr(...) = " << recovered << std::endl;

    return 0;
}*/
/* Sample output:
kevin_abraham@kevinabraham:~/NTRU$ g++ -Iinclude src/utils.cpp  -lssl -lcrypto -o bin/utils_test
kevin_abraham@kevinabraham:~/NTRU$ ./bin/utils_test
invModQ(7, 26) = 15
coeffs_to_blob: 1 2c ff eb 0 2a 3 ff 
PEM encoding:
-----BEGIN NTRU PUBLIC KEY-----
ASz/6wAqA/8=
-----END NTRU PUBLIC KEY-----
strToBin("Hi") = 0100100001101001
binToStr(...) = Hi
kevin_abraham@kevinabraham:~/NTRU$ 
*/