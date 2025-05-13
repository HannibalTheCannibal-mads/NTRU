#include <string>
#include <iostream>
#include <stdexcept>
#include "ntru/ntruparams.hpp"

NTRUParams::NTRUParams(int k, const std::string& choice) {
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
/*
int main() {
    try {
        // Test all valid parameter sets
        int levels[] = {112, 128, 192, 256};
        std::string choices[] = {"space", "hybrid", "speed"};

        for (int k : levels) {
            for (const std::string& choice : choices) {
                NTRUParams params(k, choice);
                std::cout << "NTRUParams(" << k << ", \"" << choice << "\"):\n";
                std::cout << "  N = " << params.N << "\n";
                std::cout << "  d_f = " << params.d_f << "\n";
                std::cout << "  d_g = " << params.d_g << "\n";
                std::cout << "  d_r = " << params.d_r << "\n";
                std::cout << "  q = " << params.q << "\n";
                std::cout << "  p = " << params.p << "\n";
                std::cout << std::endl;
            }
        }

        // Test invalid parameter
        try {
            NTRUParams bad(42, "speed");
        } catch (const std::runtime_error& e) {
            std::cout << "Expected error: " << e.what() << std::endl;
        }
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
        return 1;
    }
    return 0;
}*/
// This code defines the NTRUParams class, which is used to store parameters for the NTRU encryption algorithm.
// It includes a constructor that initializes the parameters based on the specified security level (k)
// and choice (space, hybrid, or speed). The parameters include N (the polynomial degree), d_f (the degree of f),
// d_g (the degree of g), d_r (the degree of r), q (the modulus), and p (the small modulus).
// The constructor throws an exception if the specified parameters are not implemented.

/*Sample output:
kevin_abraham@kevinabraham:~/NTRU$ g++ -Iinclude src/ntruparams.cpp -o bin/ntruparams_test
kevin_abraham@kevinabraham:~/NTRU$ ./bin/ntruparams_test
NTRUParams(112, "space"):
  N = 401
  d_f = 113
  d_g = 133
  d_r = 113
  q = 2048
  p = 3

NTRUParams(112, "hybrid"):
  N = 541
  d_f = 49
  d_g = 180
  d_r = 49
  q = 2048
  p = 3

NTRUParams(112, "speed"):
  N = 659
  d_f = 38
  d_g = 219
  d_r = 38
  q = 2048
  p = 3

NTRUParams(128, "space"):
  N = 449
  d_f = 134
  d_g = 149
  d_r = 134
  q = 2048
  p = 3

NTRUParams(128, "hybrid"):
  N = 613
  d_f = 55
  d_g = 204
  d_r = 55
  q = 2048
  p = 3

NTRUParams(128, "speed"):
  N = 761
  d_f = 42
  d_g = 253
  d_r = 42
  q = 2048
  p = 3

NTRUParams(192, "space"):
  N = 677
  d_f = 157
  d_g = 225
  d_r = 157
  q = 2048
  p = 3

NTRUParams(192, "hybrid"):
  N = 887
  d_f = 81
  d_g = 295
  d_r = 81
  q = 2048
  p = 3

NTRUParams(192, "speed"):
  N = 1087
  d_f = 63
  d_g = 362
  d_r = 63
  q = 2048
  p = 3

NTRUParams(256, "space"):
  N = 1087
  d_f = 120
  d_g = 362
  d_r = 120
  q = 2048
  p = 3

NTRUParams(256, "hybrid"):
  N = 1171
  d_f = 106
  d_g = 390
  d_r = 106
  q = 2048
  p = 3

NTRUParams(256, "speed"):
  N = 1499
  d_f = 79
  d_g = 499
  d_r = 79
  q = 2048
  p = 3

Expected error: Not implemented. :(
kevin_abraham@kevinabraham:~/NTRU$ 
*/