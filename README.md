# NTRU Polynomial Inversion in C++

## Overview

This project implements polynomial arithmetic and modular inversion in the ring \(\mathbb{Z}_q[x]/(x^N - 1)\), which is fundamental to the NTRU public-key cryptosystem. It features efficient cyclic polynomial multiplication, random trinary polynomial generation, and the Extended Euclidean Algorithm for computing polynomial inverses modulo \(q\) and \(x^N - 1\). The implementation includes detailed debug output to trace the inversion process step-by-step, making it suitable for educational purposes, research, and cryptographic experimentation.

## Features

- Polynomial convolution modulo \(x^N - 1\) and modulo \(q\)
- Extended Euclidean Algorithm for polynomial inversion
- Random trinary polynomial generation with coefficients in \(\{-1,0,1\}\)
- Debugging output for detailed step tracing
- Modular, clean C++ codebase

## Project Structure
```bash
├── bin
│   ├── ntrucrypto_test
│   ├── ntrukey_test
│   ├── ntruparams_test
│   ├── polymodq_test
│   └── utils_test
├── CMakeLists.txt
├── data
├── include
│   └── ntru
│       ├── convmodq.hpp
│       ├── convpoly.hpp
│       ├── ntrucrypto.hpp
│       ├── ntrukey.hpp
│       ├── ntruparams.hpp
│       ├── polymodq.hpp
│       └── utils.hpp
├── LICENSE
├── obj
├── README.md
├── s.cpp
├── src
│   ├── convmodq.cpp
│   ├── convpoly.cpp
│   ├── ntrucrypto.cpp
│   ├── ntrukey.cpp
│   ├── ntruparams.cpp
│   ├── polymodq.cpp
│   └── utils.cpp
├── test
│   └── test_ntru.cpp
└── test_convmodq_inverse 
```


## References

- [NTRUEncrypt: The NTRU Public Key Cryptosystem](https://www.ntru.org/f/hps98.pdf)
- [Wikipedia: NTRUEncrypt](https://en.wikipedia.org/wiki/NTRUEncrypt)
- [IETF NTRU Key Encapsulation draft](https://www.ietf.org/archive/id/draft-fluhrer-cfrg-ntru-02.html)

---

*For questions, contributions, or feedback, feel free to open an issue or submit a pull request.*
