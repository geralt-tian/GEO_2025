// Faithful division protocol implementation based on Algorithm 3
#include "geometric_perspective_protocols.h"

void GeometricPerspectiveProtocols::division_faithful(int32_t dim, uint64_t *input, uint64_t *output,
                                            uint64_t divisor, uint32_t bw) {
    // Implementation of Algorithm 3: Faithful division protocol
    // Input: P0 and P1 hold [x]^l and a public d
    // Output: P0 and P1 output [[int(x)/d]]^l

    uint64_t mask_bw = (bw == 64 ? -1 : ((1ULL << bw) - 1));
    uint64_t L = 1ULL << bw;  // L = 2^bw

    // Step 1: Let l_d = floor(log d)
    int l_d = 0;
    uint64_t temp_d = divisor;
    while (temp_d > 1) {
        l_d++;
        temp_d >>= 1;
    }

    // Mask for l_d+1 bits
    uint64_t mask_ld = ((l_d + 1) == 64 ? -1 : ((1ULL << (l_d + 1)) - 1));

    // Check if divisor is a power of 2 - can use optimized truncation
    bool is_power_of_2 = (divisor != 0) && ((divisor & (divisor - 1)) == 0);

    if (is_power_of_2) {
        // For powers of 2, use truncation protocol
        new_truncate(dim, input, output, l_d, bw);
        return;
    }

    // For general divisors, use full protocol
    uint64_t *mw_output = new uint64_t[dim];
    uint64_t *X1 = new uint64_t[dim];  // Will hold quotient part from lookup
    uint64_t *l_epsilon = new uint64_t[dim];  // Will hold remainder info

    // Step 6: Invoke MW protocol to learn [MW]^2
    // MW classifies input into ranges for table lookup
    mw(dim, input, mw_output, bw, 2);

    if (party == sci::ALICE) {  // P0 in the algorithm
        // Steps 2-5: Build lookup tables T^d and T^epsilon for P0
        uint64_t **T_d = new uint64_t*[dim];  // Quotient table
        uint64_t **T_epsilon = new uint64_t*[dim];  // Remainder table

        for (int i = 0; i < dim; i++) {
            T_d[i] = new uint64_t[4];  // 4 possible MW values (0,1,2,3)
            T_epsilon[i] = new uint64_t[4];

            // P0 sets T^d[j] = 0 and T^epsilon[j] = 0 for all j
            for (int j = 0; j < 4; j++) {
                T_d[i][j] = 0;
                T_epsilon[i][j] = 0;
            }
        }

        // Step 7: P0 invokes LUT(T^d, [MW]^2) to learn [X1]^l
        aux->lookup_table<uint64_t>(T_d, nullptr, nullptr, dim, 2, bw);

        // Step 8: P0 invokes LUT(T^epsilon, [MW]^2) to learn [l_epsilon]^(l_d+1)
        aux->lookup_table<uint64_t>(T_epsilon, nullptr, nullptr, dim, 2, l_d + 1);

        // Clean up tables
        for (int i = 0; i < dim; i++) {
            delete[] T_d[i];
            delete[] T_epsilon[i];
        }
        delete[] T_d;
        delete[] T_epsilon;

        // P0's share of X1 and l_epsilon is 0 (from lookup tables)
        for (int i = 0; i < dim; i++) {
            X1[i] = 0;
            l_epsilon[i] = 0;
        }

        // Step 9: P0 computes temp = (x0 mod d) + l_epsilon - d
        uint64_t *temp_val = new uint64_t[dim];
        for (int i = 0; i < dim; i++) {
            uint64_t x0_mod_d = input[i] % divisor;
            // Compute in larger bitwidth to handle overflow properly
            int64_t temp_signed = (int64_t)x0_mod_d + (int64_t)l_epsilon[i] - (int64_t)divisor;
            temp_val[i] = temp_signed & mask_ld;
        }

        // Step 10: DReLU to detect carry (if temp < 0, e = 0, else e = 1)
        uint8_t *e_bool = new uint8_t[dim];
        // Using MSB to check sign - if MSB is 1, temp is negative
        aux->MSB(temp_val, e_bool, dim, l_d + 1);

        // Invert for DReLU (DReLU returns 0 if negative, 1 if positive)
        for (int i = 0; i < dim; i++) {
            e_bool[i] = 1 - e_bool[i];
        }

        // Step 11: B2A conversion
        uint64_t *e_arith = new uint64_t[dim];
        aux->B2A(e_bool, e_arith, dim, bw);

        // Step 12: P0 outputs floor(x0/d) + [X1] + [e]
        for (int i = 0; i < dim; i++) {
            uint64_t x0_div_d = input[i] / divisor;  // Local division of P0's share
            output[i] = (x0_div_d + X1[i] + e_arith[i]) & mask_bw;
        }

        delete[] temp_val;
        delete[] e_bool;
        delete[] e_arith;

    } else { // party == BOB (P1 in the algorithm)
        // Steps 2-5: Build lookup tables T^d and T^epsilon for P1
        uint64_t **T_d = new uint64_t*[dim];
        uint64_t **T_epsilon = new uint64_t*[dim];

        for (int i = 0; i < dim; i++) {
            T_d[i] = new uint64_t[4];
            T_epsilon[i] = new uint64_t[4];

            // P1 sets T^d[j] = floor((x1 - j*L)/d) and
            // T^epsilon[j] = (x1 - j*L) mod d
            for (int j = 0; j < 3; j++) {  // j ∈ {0,1,2}
                // Compute (x1 - j*L) with proper wrapping
                int64_t val_signed = (int64_t)input[i] - (int64_t)(j * L);
                uint64_t val = val_signed & mask_bw;

                T_d[i][j] = (val / divisor) & mask_bw;
                T_epsilon[i][j] = (val % divisor) & mask_ld;
            }
            T_d[i][3] = 0;  // j=3 case (overflow)
            T_epsilon[i][3] = 0;
        }

        // Step 7: P1 invokes LUT to learn [X1]^l
        aux->lookup_table<uint64_t>(nullptr, mw_output, X1, dim, 2, bw);

        // Step 8: P1 invokes LUT to learn [l_epsilon]^(l_d+1)
        aux->lookup_table<uint64_t>(nullptr, mw_output, l_epsilon, dim, 2, l_d + 1);

        // Clean up tables (they were already used by aux)
        for (int i = 0; i < dim; i++) {
            delete[] T_d[i];
            delete[] T_epsilon[i];
        }
        delete[] T_d;
        delete[] T_epsilon;

        // Step 9: P1 computes temp = l_epsilon - d
        uint64_t *temp_val = new uint64_t[dim];
        for (int i = 0; i < dim; i++) {
            int64_t temp_signed = (int64_t)l_epsilon[i] - (int64_t)divisor;
            temp_val[i] = temp_signed & mask_ld;
        }

        // Step 10: DReLU to detect carry
        uint8_t *e_bool = new uint8_t[dim];
        aux->MSB(temp_val, e_bool, dim, l_d + 1);

        // Invert for DReLU
        for (int i = 0; i < dim; i++) {
            e_bool[i] = 1 - e_bool[i];
        }

        // Step 11: B2A conversion
        uint64_t *e_arith = new uint64_t[dim];
        aux->B2A(e_bool, e_arith, dim, bw);

        // Step 12: P1 outputs [X1] + [e]
        for (int i = 0; i < dim; i++) {
            output[i] = (X1[i] + e_arith[i]) & mask_bw;
        }

        delete[] temp_val;
        delete[] e_bool;
        delete[] e_arith;
    }

    delete[] mw_output;
    delete[] X1;
    delete[] l_epsilon;
}