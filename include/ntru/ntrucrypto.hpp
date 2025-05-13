#ifndef NTRU_NTRUCRYPTO_HPP
#define NTRU_NTRUCRYPTO_HPP

#include "ntruparams.hpp"
#include "ntrukey.hpp"
#include "convmodq.hpp"
#include "convpoly.hpp"
#include <vector>
#include <string>
#include <utility>

std::string NTRUDecrypt(const NTRUKey& key, const std::vector<ConvModQ>& cl, int n);

std::pair<std::vector<ConvModQ>, int> NTRUEncrypt(
    const NTRUParams& ring,
    const ConvModQ& pub,
    const std::string& m);

ConvModQ NTRUBlockEncrypt(
    const NTRUParams& ring,
    const ConvModQ& h,
    const ConvPoly& m);

ConvPoly NTRUBlockDecrypt(
    const NTRUKey& key,
    const ConvModQ& c);

#endif // NTRU_NTRUCRYPTO_HPP
