#ifndef NTRU_NTRUPARAMS_HPP
#define NTRU_NTRUPARAMS_HPP

#include <string>
#include <stdexcept>

class NTRUParams {
public:
    int N, d_f, d_g, d_r, q, p;

    NTRUParams(int k, const std::string& choice = "speed");
};

#endif // NTRU_NTRUPARAMS_HPP
