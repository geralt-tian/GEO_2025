/*
 * Test for matrix_vector_unsigned_mul
 * Note: Unlike crossterm, inA and inB are SHARED values, outC is also shared
 * All true values are positive numbers
 */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <random>
#include <chrono>
#include <iomanip>
#include <cstring>

#include "LinearOT/linear-ot.h"
#include "utils/emp-tool.h"
#include "FloatingPoint/floating-point.h"
#include "FloatingPoint/fp-math.h"
#include "BuildingBlocks/aux-protocols.h"
#include "BuildingBlocks/geometric_perspective_protocols.h"
#include "NonLinear/relu-ring.h"
#include "NonLinear/relu-field.h"
#include "Math/math-functions.h"
#include "Millionaire/millionaire.h"
#include "Millionaire/millionaire_with_equality.h"

using namespace sci;
using namespace std;

int party, port = 32000;
string address = "127.0.0.1";
IOPack *iopack;
OTPack *otpack;
LinearOT *prod;
GeometricPerspectiveProtocols *gp;
AuxProtocols *aux;

// Test parameters for matrix_vector_unsigned_mul
int32_t m = 1;          // Matrix rows (and vector length)
int32_t n = 100;         // Matrix columns
int32_t bwA = 16;        // Matrix bit width
int32_t bwB = 32;        // Vector bit width  
int32_t bwC = bwA + bwB;       // Output bit width

uint64_t mask_A = (1ULL << bwA) - 1;
uint64_t mask_B = (1ULL << bwB) - 1;
uint64_t mask_C = (1ULL << bwC) - 1;

int main(int argc, char **argv) {
    ArgMapping amap;
    amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
    amap.arg("p", port, "Port Number");
    amap.arg("ip", address, "IP Address of server (ALICE)");
    amap.parse(argc, argv);

    cout << "======================================" << endl;
    cout << "=== matrix_vector_unsigned_mul Test Suite ===" << endl;
    cout << "Testing matrix-vector unsigned multiplication" << endl;
    cout << "Party: " << (party == sci::ALICE ? "ALICE" : "BOB") << endl;
    cout << "Matrix dimensions: " << m << " x " << n << endl;
    cout << "Vector length: " << m << endl;
    cout << "Matrix bit width: " << bwA << endl;
    cout << "Vector bit width: " << bwB << endl;
    cout << "Output bit width: " << bwC << endl;
    cout << "======================================" << endl;

    // Initialize communication
    iopack = new IOPack(party, port, address);
    otpack = new OTPack(iopack, party);
    prod = new LinearOT(party, iopack, otpack);
    gp = new GeometricPerspectiveProtocols(party, iopack, otpack);
    aux = new AuxProtocols(party, iopack, otpack);

    // Allocate arrays for shares
    uint64_t *inA = new uint64_t[m * n];    // Matrix input shares (m x n)
    uint64_t *inB = new uint64_t[m];        // Vector input shares (m)
    uint64_t *outC = new uint64_t[m * n];   // Output shares (m x n)

    // Generate test data with fixed seed for reproducibility
    std::mt19937 gen(1234);
    
    // Generate positive cleartext values (ensure all are positive)
    uint64_t *cleartext_matrix = new uint64_t[m * n];
    uint64_t *cleartext_vector = new uint64_t[m];
    
    // 在模环上，前50%的值被解释为正数：[0, 2^{n-1} - 1]
    // 我们使用 [1, 2^{n-1} - 1] 确保都是正数且不为0
    uint64_t max_positive_A = (1ULL << (bwA - 1)) - 1;  // 2^{bwA-1} - 1
    uint64_t max_positive_B = (1ULL << (bwB - 1)) - 1;  // 2^{bwB-1} - 1
    
    std::uniform_int_distribution<uint64_t> mat_dis(1, max_positive_A);    // Matrix: 1 to 2^{bwA-1} - 1
    std::uniform_int_distribution<uint64_t> vec_dis(1, max_positive_B);    // Vector: 1 to 2^{bwB-1} - 1

    cout << "Generating positive cleartext data..." << endl;
    cout << "Matrix range: [1, " << max_positive_A << "] (positive half of " << bwA << "-bit ring)" << endl;
    cout << "Vector range: [1, " << max_positive_B << "] (positive half of " << bwB << "-bit ring)" << endl;

    for (int i = 0; i < m * n; i++) {
        cleartext_matrix[i] = mat_dis(gen);
        // 确保在环的正数部分：[1, 2^{bwA-1} - 1]
        assert(cleartext_matrix[i] > 0 && cleartext_matrix[i] <= max_positive_A);
    }
    
    for (int i = 0; i < m; i++) {
        cleartext_vector[i] = vec_dis(gen);
        // 确保在环的正数部分：[1, 2^{bwB-1} - 1]
        assert(cleartext_vector[i] > 0 && cleartext_vector[i] <= max_positive_B);
    }

    // Create shares: cleartext = share_alice + share_bob (mod 2^bw)
    uint64_t *matrix_share_alice = new uint64_t[m * n];
    uint64_t *matrix_share_bob = new uint64_t[m * n];
    uint64_t *vector_share_alice = new uint64_t[m];
    uint64_t *vector_share_bob = new uint64_t[m];

    cout << "Creating secret shares..." << endl;
    
    // 标准的additive secret sharing：alice的share完全随机，bob的share保证重构正确性
    for (int i = 0; i < m * n; i++) {
        matrix_share_alice[i] = gen() & mask_A;  // Alice的share完全随机
        matrix_share_bob[i] = (cleartext_matrix[i] - matrix_share_alice[i]) & mask_A;  // Bob的share确保重构正确
        
        // 验证重构正确性
        uint64_t reconstructed = (matrix_share_alice[i] + matrix_share_bob[i]) & mask_A;
        assert(reconstructed == cleartext_matrix[i]);
    }
    
    for (int i = 0; i < m; i++) {
        vector_share_alice[i] = gen() & mask_B;  // Alice的share完全随机
        vector_share_bob[i] = (cleartext_vector[i] - vector_share_alice[i]) & mask_B;  // Bob的share确保重构正确
        
        // 验证重构正确性
        uint64_t reconstructed = (vector_share_alice[i] + vector_share_bob[i]) & mask_B;
        assert(reconstructed == cleartext_vector[i]);
    }

    // Assign shares to each party
    if (party == sci::ALICE) {
        memcpy(inA, matrix_share_alice, m * n * sizeof(uint64_t));
        memcpy(inB, vector_share_alice, m * sizeof(uint64_t));
    } else {
        memcpy(inA, matrix_share_bob, m * n * sizeof(uint64_t));
        memcpy(inB, vector_share_bob, m * sizeof(uint64_t));
    }

    // Print sample cleartext data for verification
    if (party == sci::ALICE) {
        cout << "Cleartext Matrix A (first 3x3 elements):" << endl;
        for (int i = 0; i < min(3, m); i++) {
            for (int j = 0; j < min(3, n); j++) {
                printf("%4lu ", cleartext_matrix[i * n + j]);
            }
            printf("\n");
        }
        
        cout << "Cleartext Vector B (first 3 elements): ";
        for (int i = 0; i < min(3, m); i++) {
            printf("%lu ", cleartext_vector[i]);
        }
        printf("\n");
    }

    cout << "Starting matrix_vector_unsigned_mul computation..." << endl;
    
    // Record communication start
    size_t comm_start = iopack->io->counter;
    auto start_time = chrono::high_resolution_clock::now();

    // Call the matrix_vector_unsigned_mul function
    gp->matrix_vector_unsigned_mul(m, n, inA, inB, outC, bwA, bwB, bwC);

    auto end_time = chrono::high_resolution_clock::now();
    size_t comm_end = iopack->io->counter;
    
    auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
 
    cout << "Computation completed in " << duration.count() << " ms" << endl;

    // Share results for accuracy testing
    uint64_t *result_alice = new uint64_t[m * n];
    uint64_t *result_bob = new uint64_t[m * n];
    
    if (party == sci::ALICE) {
        iopack->io->send_data(outC, m * n * sizeof(uint64_t));
        iopack->io->recv_data(result_bob, m * n * sizeof(uint64_t));
        memcpy(result_alice, outC, m * n * sizeof(uint64_t));
    } else {
        iopack->io->recv_data(result_alice, m * n * sizeof(uint64_t));
        iopack->io->send_data(outC, m * n * sizeof(uint64_t));
        memcpy(result_bob, outC, m * n * sizeof(uint64_t));
    }

    // Accuracy analysis (only on one party to avoid duplicate output)
    if (party == sci::ALICE) {
        cout << "\n=== Accuracy Analysis ===" << endl;
        
        double total_error = 0.0;
        double max_error = 0.0;
        int correct_count = 0;
        int total_tests = m * n;

        cout << "Reconstructing and verifying results..." << endl;
        
        // Calculate expected results and compare
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                int idx = i * n + j;
                
                // Reconstruct the result
                uint64_t reconstructed = (result_alice[idx] + result_bob[idx]) & mask_C;
                
                                 // Calculate expected result: matrix[i][j] * vector[i] (真正的矩阵-向量乘法)
                uint64_t expected = (cleartext_matrix[idx] * cleartext_vector[i]) & mask_C;
                
                // Calculate error
                uint64_t error = (reconstructed > expected) ? 
                    (reconstructed - expected) : (expected - reconstructed);
                
                total_error += error;
                if (error > max_error) max_error = error;
                if (error == 0) correct_count++;
                
                // Print first few results for verification
                if (i < 3 && j < 3) {
                    printf("Element [%d][%d]: matrix=%lu * vector=%lu = expected=%lu, got=%lu, error=%lu\n",
                           i, j, cleartext_matrix[idx], cleartext_vector[i], 
                           expected, reconstructed, error);
                }
            }
        }
        
        double avg_error = total_error / total_tests;
        double accuracy = (double)correct_count / total_tests * 100.0;
        
        cout << fixed << setprecision(2);
        cout << "Total tests: " << total_tests << endl;
        cout << "Correct results: " << correct_count << " (" << accuracy << "%)" << endl;
        cout << "Average error: " << avg_error << endl;
        cout << "Maximum error: " << max_error << endl;
        
        if (accuracy >= 99.0) {
            cout << "✅ Test PASSED - High accuracy achieved!" << endl;
        } else if (accuracy >= 95.0) {
            cout << "⚠️  Test MARGINAL - Acceptable accuracy" << endl;
        } else {
            cout << "❌ Test FAILED - Low accuracy" << endl;
        }
    }

    // Performance and communication analysis
    size_t communication_bytes = comm_end - comm_start;
    
    cout << "\n=== Performance Analysis ===" << endl;
    cout << "Execution time: " << duration.count() << " ms" << endl;
    cout << "Communication: " << communication_bytes << " bytes" << endl;
    cout << "Communication per element: " << (double)communication_bytes / (m * n) << " bytes" << endl;
    
    // Theoretical analysis
    cout << "\n=== Protocol Analysis ===" << endl;
    cout << "Expected operations:" << endl;
    cout << "- Matrix elements: " << m * n << endl;
    cout << "- Vector elements: " << m << endl;
    cout << "- Output elements: " << m * n << endl;
    cout << "- Bit width combinations: bwA=" << bwA << ", bwB=" << bwB << ", bwC=" << bwC << endl;
    
    // Estimate expected communication for unsigned multiplication
    // This depends on the specific protocol implementation
    int expected_rounds = bwA + bwB;  // Typical for multiplication protocols
    size_t expected_comm = m * n * expected_rounds * 16;  // Rough estimate
    cout << "Estimated communication: ~" << expected_comm << " bytes" << endl;
    
    double comm_efficiency = (double)expected_comm / communication_bytes;
    cout << "Communication efficiency: " << fixed << setprecision(2) << comm_efficiency << "x" << endl;

    // Cleanup
    delete[] inA;
    delete[] inB;
    delete[] outC;
    delete[] cleartext_matrix;
    delete[] cleartext_vector;
    delete[] matrix_share_alice;
    delete[] matrix_share_bob;
    delete[] vector_share_alice;
    delete[] vector_share_bob;
    delete[] result_alice;
    delete[] result_bob;

    delete iopack;
    delete otpack;
    delete prod;
    delete gp;
    delete aux;

    cout << "\n=== Test completed successfully ===" << endl;
    return 0;
}
