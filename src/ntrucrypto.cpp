#include <random>
#include <algorithm>
#include <utility>
#include <vector>
#include <string>
#include "ntru/ntrucrypto.hpp"
#include "ntru/ntrukey.hpp"
#include "ntru/ntruparams.hpp"
#include "ntru/convmodq.hpp"
#include "ntru/convpoly.hpp"
#include "ntru/utils.hpp"

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
    ConvModQ c((ring.p * r * h + ConvModQ(m.coef, ring.q, ring.N)).coef, ring.q, ring.N);
    return c;
}

ConvPoly NTRUBlockDecrypt(const NTRUKey& key, const ConvModQ& c) {
    ConvModQ a = key.f * c;
    ConvPoly aprime = a.centerLift();
    ConvModQ m = key.finvp * aprime;
    return m.centerLift();
}
// This function is a placeholder for the actual NTRU decryption process.
// It takes a key and a ciphertext, and returns the decrypted message.
// The actual decryption logic would depend on the NTRU algorithm specifics.
// The function NTRUBlockEncrypt is responsible for encrypting a block of data using the NTRU algorithm.
// It takes the NTRU parameters, a public key, and a message, and returns the ciphertext.
#include "ntru/utils.hpp"

int main() {
    // 1. Set up NTRU parameters and generate a keypair
    NTRUParams params(112, "speed");
    NTRUKey key(params);

    // 2. Get the public key
    auto pub = key.publicKey();

    // --- Print PRIVATE KEY in PEM format ---
    std::cout << "\n--- NTRU PRIVATE KEY (f) ---\n";
    std::cout << to_pem(key.f.coef, "NTRU PRIVATE KEY", 2) << "\n";

    std::cout << "\n--- NTRU PRIVATE KEY (g) ---\n";
    std::cout << to_pem(key.g.coef, "NTRU PRIVATE KEY G", 2) << "\n";

    // You may also want to print finvq and finvp if desired
    std::cout << to_pem(key.finvq.coef, "NTRU PRIVATE KEY FINVQ", 2) << "\n";
    std::cout << to_pem(key.finvp.coef, "NTRU PRIVATE KEY FINVP", 2) << "\n";

    // --- Print PUBLIC KEY in PEM format ---
    std::cout << "\n--- NTRU PUBLIC KEY (h) ---\n";
    std::cout << to_pem(pub.second.coef, "NTRU PUBLIC KEY", 2) << "\n";

    // 3. Message to encrypt
    std::string message = "Hallo Hannibal";

    // 4. Encrypt the message
    auto enc_result = NTRUEncrypt(params, pub.second, message);
    std::vector<ConvModQ> ciphertext = enc_result.first;
    int n = enc_result.second;

    std::cout << "Original message: " << message << std::endl;
    std::cout << "Ciphertext blocks: " << ciphertext.size() << std::endl;

    // 5. Decrypt the message
    std::string decrypted = NTRUDecrypt(key, ciphertext, n);

    std::cout << "Decrypted message: " << decrypted << std::endl;

    // 6. Check if successful
    if (message == decrypted) {
        std::cout << "Encryption and decryption successful!" << std::endl;
    } else {
        std::cout << "Decryption failed." << std::endl;
    }
    std::cout << "h coefficients: ";
    for (int c : pub.second.coef) std::cout << c << " ";
    std::cout << std::endl;
    std::cout << "f: ";
    for (int c : key.f.coef) std::cout << c << " ";
    std::cout << "\ng: ";
    for (int c : key.g.coef) std::cout << c << " ";
    std::cout << "\nfinvq: ";
    for (int c : key.finvq.coef) std::cout << c << " ";
    std::cout << "\n";

    return 0;
}
