#ifndef NTRU_NTRUKEY_HPP
#define NTRU_NTRUKEY_HPP

#include "ntruparams.hpp"
#include "convmodq.hpp"
#include <utility>

class NTRUKey {
public:
    NTRUParams ring;
    ConvModQ f, g, finvq, finvp, h;

    NTRUKey(const NTRUParams& ring_param = NTRUParams(256, "speed"),
            const ConvModQ* f_ = nullptr,
            const ConvModQ* g_ = nullptr);

    ConvModQ randomTrinary(int d1, int d2);

    std::pair<NTRUParams, ConvModQ> publicKey() const;
};

#endif // NTRU_NTRUKEY_HPP
