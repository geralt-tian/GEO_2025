// Test file for matrix multiplication function
// This file tests the matrix_multiplication function which computes matrix multiplication
// of two secret-shared matrices using secure multiparty computation  
// Including accuracy verification, communication measurement, and timing

#include <cstdint>
#include <cstdio>
#include <iostream>
#include <chrono>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <limits>
#include <random>
#include <cstring>
#include <iomanip>

#include "LinearOT/linear-ot.h"
#include "utils/emp-tool.h"
#include "FloatingPoint/floating-point.h"
#include "FloatingPoint/fp-math.h"
#include "BuildingBlocks/aux-protocols.h"
#include "BuildingBlocks/geometric_perspective_protocols.h"
#include "BuildingBlocks/truncation.h"
#include "Math/math-functions.h"
#include "Millionaire/equality.h"
#include "Millionaire/millionaire.h"
#include "Millionaire/millionaire_with_equality.h"

using namespace sci;
using namespace std;

int party, port = 32000;
string address = "127.0.0.1";
IOPack *iopack;
OTPack *otpack;
LinearOT *mult;
GeometricPerspectiveProtocols *gp;
AuxProtocols *aux;

// Test parameters for matrix multiplication
int32_t rows_A = 1;      // Matrix A rows (100)
int32_t common_dim = 2;  // Common dimension (200) 
int32_t cols_B = 2;      // Matrix B columns (300)
int32_t bwA = 16;          // Input A bit width
int32_t bwB = 32;          // Input B bit width  
int32_t bwC = bwA + bwB;   // Output bit width
bool signed_arithmetic = true;
bool signed_B = true;
MultMode mode = MultMode::None;

uint64_t mask_A = (bwA == 64 ? -1 : ((1ULL << bwA) - 1));
uint64_t mask_B = (bwB == 64 ? -1 : ((1ULL << bwB) - 1));
uint64_t mask_C = (bwC == 64 ? -1 : ((1ULL << bwC) - 1));

// Function to convert to signed representation
int64_t to_signed(uint64_t val, int bw) {
    uint64_t sign_mask = 1ULL << (bw - 1);
    if (val & sign_mask) {
        return (int64_t)(val | (~((1ULL << bw) - 1)));
    }
    return (int64_t)val;
}

// Function to compute plaintext matrix multiplication for verification
void compute_matrix_mult_plaintext(uint64_t *matA, uint64_t *matB, uint64_t *result, 
                                   int32_t rows_A, int32_t common_dim, int32_t cols_B, 
                                   int32_t bwA, int32_t bwB, int32_t bwC,
                                   bool signed_arithmetic, bool signed_B) {
    // matA: rows_A x common_dim
    // matB: common_dim x cols_B  
    // result: rows_A x cols_B
    for (int i = 0; i < rows_A; i++) {
        for (int j = 0; j < cols_B; j++) {
            int64_t sum = 0;
            for (int k = 0; k < common_dim; k++) {
                if (signed_arithmetic) {
                    int64_t a = to_signed(matA[i * common_dim + k], bwA);
                    int64_t b = signed_B ? to_signed(matB[k * cols_B + j], bwB) : (int64_t)matB[k * cols_B + j];
                    sum += a * b;
                } else {
                    uint64_t a = matA[i * common_dim + k] & mask_A;
                    uint64_t b = matB[k * cols_B + j] & mask_B;
                    sum += a * b;
                }
            }
            result[i * cols_B + j] = (uint64_t)sum & mask_C;
        }
    }
}

// Function to reconstruct secret shares
void reconstruct_shares(uint64_t *share_alice, uint64_t *share_bob, uint64_t *result, 
                        int32_t dim, int32_t bw) {
    uint64_t mask = (bw == 64 ? -1 : ((1ULL << bw) - 1));
    for (int i = 0; i < dim; i++) {
        result[i] = (share_alice[i] + share_bob[i]) & mask;
    }
}

// Function to generate random shares for secret sharing
void generate_shares(uint64_t *original, uint64_t *share_alice, uint64_t *share_bob,
                     int32_t dim, int32_t bw) {
    uint64_t mask = (bw == 64 ? -1 : ((1ULL << bw) - 1));
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint64_t> dis(0, mask);
    
    for (int i = 0; i < dim; i++) {
        share_alice[i] = dis(gen) & mask;
        share_bob[i] = (original[i] - share_alice[i]) & mask;
    }
}

int main(int argc, char **argv) {
    ArgMapping amap;
    amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
    amap.arg("p", port, "Port Number");
    amap.arg("ip", address, "IP Address of server (ALICE)");
    amap.arg("rows_A", rows_A, "Matrix A rows");
    amap.arg("common_dim", common_dim, "Common dimension");
    amap.arg("cols_B", cols_B, "Matrix B columns");
    amap.arg("bwA", bwA, "Input A bit width");
    amap.arg("bwB", bwB, "Input B bit width");
    amap.arg("signed", signed_arithmetic, "Use signed arithmetic (0/1)");
    amap.arg("signedB", signed_B, "Treat B as signed (0/1)");
    amap.parse(argc, argv);

    // Update dependent parameters
    bwC = bwA + bwB;  // Standard output bitwidth for multiplication
    mask_A = (bwA == 64 ? -1 : ((1ULL << bwA) - 1));
    mask_B = (bwB == 64 ? -1 : ((1ULL << bwB) - 1));
    mask_C = (bwC == 64 ? -1 : ((1ULL << bwC) - 1));

    cout << "=== Matrix Multiplication Test Suite ===" << endl;
    cout << "Testing matrix multiplication of secret-shared matrices" << endl;
    cout << "Party: " << (party == sci::ALICE ? "ALICE" : "BOB") << endl;
    cout << "Matrix A dimensions: " << rows_A << " x " << common_dim << endl;
    cout << "Matrix B dimensions: " << common_dim << " x " << cols_B << endl;
    cout << "Result matrix dimensions: " << rows_A << " x " << cols_B << endl;
    cout << "Input A bit width: " << bwA << endl;
    cout << "Input B bit width: " << bwB << endl;
    cout << "Output bit width: " << bwC << endl;
    cout << "Signed arithmetic: " << (signed_arithmetic ? "Yes" : "No") << endl;
    cout << "Signed B: " << (signed_B ? "Yes" : "No") << endl;

    // Initialize communication
    iopack = new IOPack(party, port, address);
    otpack = new OTPack(iopack, party);
    mult = new LinearOT(party, iopack, otpack);
    gp = new GeometricPerspectiveProtocols(party, iopack, otpack);
    aux = new AuxProtocols(party, iopack, otpack);

    // Calculate matrix sizes
    int32_t size_A = rows_A * common_dim;     // Matrix A size
    int32_t size_B = common_dim * cols_B;     // Matrix B size  
    int32_t size_C = common_dim * cols_B;         // Result matrix size

    // Allocate arrays
    uint64_t *inA = new uint64_t[size_A];      // Matrix A shares
    uint64_t *inB = new uint64_t[size_B];      // Matrix B shares
    uint64_t *outC = new uint64_t[size_C];     // Output shares

    // Arrays for verification (only used by ALICE for coordination)
    uint64_t *test_matA = new uint64_t[size_A];
    uint64_t *test_matB = new uint64_t[size_B];
    uint64_t *expected_result = new uint64_t[size_C];
    uint64_t *actual_result = new uint64_t[size_C];
    uint64_t *shareA_A = new uint64_t[size_A];
    uint64_t *shareB_A = new uint64_t[size_B];
    uint64_t *shareA_B = new uint64_t[size_A];
    uint64_t *shareB_B = new uint64_t[size_B];
    uint64_t *shareC_A = new uint64_t[size_C];
    uint64_t *shareC_B = new uint64_t[size_C];

    // Generate test data with fixed seed for reproducibility
    std::mt19937 gen(42);
    std::uniform_int_distribution<uint64_t> matA_dis(1, mask_A >> 1);  // Smaller range for signed
    std::uniform_int_distribution<uint64_t> matB_dis(1, mask_B >> 1);

    cout << "Generating test data for " << rows_A << "x" << common_dim << " and " 
         << common_dim << "x" << cols_B << " matrices..." << endl;

    if (party == sci::ALICE) {
        // ALICE generates the original matrices and creates shares
        for (int i = 0; i < size_A; i++) {
            test_matA[i] = matA_dis(gen);
        }
        for (int i = 0; i < size_B; i++) {
            test_matB[i] = matB_dis(gen);
        }

        // Generate shares
        generate_shares(test_matA, shareA_A, shareA_B, size_A, bwA);
        generate_shares(test_matB, shareB_A, shareB_B, size_B, bwB);

        // Send BOB's shares
        iopack->io->send_data(shareA_B, size_A * sizeof(uint64_t));
        iopack->io->send_data(shareB_B, size_B * sizeof(uint64_t));

        // ALICE uses her shares
        memcpy(inA, shareA_A, size_A * sizeof(uint64_t));
        memcpy(inB, shareB_A, size_B * sizeof(uint64_t));

        // Compute expected result for verification
        compute_matrix_mult_plaintext(test_matA, test_matB, expected_result, 
                                      rows_A, common_dim, cols_B, bwA, bwB, bwC, 
                                      signed_arithmetic, signed_B);

        cout << "Test matrices generated and shares distributed." << endl;
        cout << "Sample original values (first row of A, first column of B):" << endl;
        cout << "MatA[0]: ";
        for (int i = 0; i < std::min(5, common_dim); i++) {
            if (signed_arithmetic) {
                cout << to_signed(test_matA[i], bwA) << " ";
            } else {
                cout << test_matA[i] << " ";
            }
        }
        cout << endl;
        cout << "MatB[:][0]: ";
        for (int i = 0; i < std::min(5, common_dim); i++) {
            if (signed_B) {
                cout << to_signed(test_matB[i * cols_B], bwB) << " ";
            } else {
                cout << test_matB[i * cols_B] << " ";
            }
        }
        cout << endl;
    } else {
        // BOB receives his shares
        iopack->io->recv_data(inA, size_A * sizeof(uint64_t));
        iopack->io->recv_data(inB, size_B * sizeof(uint64_t));

        cout << "Received shares from ALICE." << endl;
    }

    // Synchronize before testing
    iopack->io->flush();

    cout << "Starting matrix multiplication computation..." << endl;
    auto start_time = chrono::high_resolution_clock::now();

    // Execute matrix_multiplication function
    mult->matrix_multiplication(rows_A, common_dim, cols_B, inA, inB, outC, bwA, bwB, bwC, 
                                signed_arithmetic, signed_B, false, mode, nullptr, nullptr);
    
    
    for (int i = 0; i < 2; i++) {
        printf("inA[%d] = %llu\n", i, inA[i]);
    }

    for (int i = 0; i < 4; i++) {
        printf("inB[%d] = %llu\n", i, inB[i]);
        printf("outC[%d] = %llu\n", i, outC[i]);
    }
    
    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end_time - start_time);

    cout << "Matrix multiplication computation completed." << endl;
    cout << "Computation time: " << duration.count() << " microseconds" << endl;

    // Verification phase
    if (party == sci::ALICE) {
        // Send output shares to BOB for reconstruction
        iopack->io->send_data(outC, size_C * sizeof(uint64_t));
        memcpy(shareC_A, outC, size_C * sizeof(uint64_t));

        // Receive BOB's output shares
        iopack->io->recv_data(shareC_B, size_C * sizeof(uint64_t));

        // Reconstruct the result
        reconstruct_shares(shareC_A, shareC_B, actual_result, size_C, bwC);

        // Verify correctness
        int errors = 0;
        int max_errors_to_show = 10;
        cout << "\n=== Verification Results ===" << endl;
        cout << "Comparing MPC result with expected plaintext result..." << endl;

        for (int i = 0; i < size_C; i++) {
            if (actual_result[i] != expected_result[i]) {
                errors++;
                if (errors <= max_errors_to_show) {
                    cout << "ERROR at index " << i << ": ";
                    cout << "Expected " << expected_result[i] 
                         << ", Got " << actual_result[i] << endl;
                }
            }
        }

        if (errors == 0) {
            cout << "✓ SUCCESS: All " << size_C << " elements match expected results!" << endl;
        } else {
            cout << "✗ FAILURE: " << errors << " out of " << size_C 
                 << " elements are incorrect (" << (100.0 * errors / size_C) << "%)" << endl;
        }

        // Show some sample results
        cout << "\nSample results (first 5 elements of result matrix C):" << endl;
        cout << "Position | Expected | Actual | Status" << endl;
        cout << "---------|----------|--------|--------" << endl;
        for (int i = 0; i < std::min(5, size_C); i++) {
            int row = i / cols_B;
            int col = i % cols_B;
            cout << "C[" << std::setw(2) << row << "," << std::setw(2) << col << "] | ";
            cout << std::setw(8) << expected_result[i] << " | ";
            cout << std::setw(6) << actual_result[i] << " | ";
            cout << (actual_result[i] == expected_result[i] ? "OK" : "ERR") << endl;
        }

    } else {
        // BOB receives ALICE's output shares and sends his own
        iopack->io->recv_data(shareC_A, size_C * sizeof(uint64_t));
        iopack->io->send_data(outC, size_C * sizeof(uint64_t));
        
        cout << "Output shares exchanged for verification." << endl;
    }

    // Communication statistics
    cout << "\n=== Communication Statistics ===" << endl;
    cout << "Total bytes sent: " << iopack->io->counter << endl;
    cout << "Bytes per element: " << (double)iopack->io->counter / size_C << endl;

    // Cleanup
    if (inA) delete[] inA;
    if (inB) delete[] inB;
    if (outC) delete[] outC;
    if (test_matA) delete[] test_matA;
    if (test_matB) delete[] test_matB;
    if (expected_result) delete[] expected_result;
    if (actual_result) delete[] actual_result;
    if (shareA_A) delete[] shareA_A;
    if (shareB_A) delete[] shareB_A;
    if (shareA_B) delete[] shareA_B;
    if (shareB_B) delete[] shareB_B;
    if (shareC_A) delete[] shareC_A;
    if (shareC_B) delete[] shareC_B;

    if (aux) delete aux;
    if (gp) delete gp;
    if (mult) delete mult;
    if (otpack) delete otpack;
    if (iopack) delete iopack;

    cout << "\n=== Test Completed ===" << endl;
    return 0;
}
