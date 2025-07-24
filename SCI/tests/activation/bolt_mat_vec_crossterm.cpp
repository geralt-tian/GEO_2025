// Test file for matrix_vector_crossterm function
// This file tests the matrix_vector_crossterm function which computes matrix-vector crossterm
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
LinearOT *prod;
GeometricPerspectiveProtocols *gp;
AuxProtocols *aux;

// Test parameters for matrix_vector_crossterm
int32_t m = 50;         // Matrix rows (and vector length)
int32_t n = 20;        // Matrix columns
int32_t bwA = 20;       // Matrix element bit width
int32_t bwB = 20;       // Vector element bit width  
int32_t bwC = bwA + bwB;  // Output bit width

uint64_t mask_A = (bwA == 64 ? -1 : ((1ULL << bwA) - 1));
uint64_t mask_B = (bwB == 64 ? -1 : ((1ULL << bwB) - 1));
uint64_t mask_C = (bwC == 64 ? -1 : ((1ULL << bwC) - 1));

int main(int argc, char **argv) {
    ArgMapping amap;
    amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
    amap.arg("p", port, "Port Number");
    amap.arg("ip", address, "IP Address of server (ALICE)");
    amap.arg("m", m, "Number of matrix rows (and vector length)");
    amap.arg("n", n, "Number of matrix columns");
    amap.arg("bwA", bwA, "Matrix element bit width");
    amap.arg("bwB", bwB, "Vector element bit width");
    amap.parse(argc, argv);

    cout << "=== matrix_vector_crossterm Test Suite ===" << endl;
    cout << "Testing matrix-vector crossterm computation" << endl;
    cout << "Party: " << (party == sci::ALICE ? "ALICE" : "BOB") << endl;
    cout << "Matrix dimensions: " << m << " x " << n << endl;
    cout << "Vector length: " << m << endl;
    cout << "Matrix bit width: " << bwA << endl;
    cout << "Vector bit width: " << bwB << endl;
    cout << "Output bit width: " << bwC << endl;

    // Initialize communication
    iopack = new IOPack(party, port, address);
    otpack = new OTPack(iopack, party);
    prod = new LinearOT(party, iopack, otpack);
    gp = new GeometricPerspectiveProtocols(party, iopack, otpack);
    aux = new AuxProtocols(party, iopack, otpack);

    // Allocate arrays
    uint64_t *inA = new uint64_t[m * n];    // Matrix input (m x n)
    uint64_t *inB = new uint64_t[m];        // Vector input (m)
    uint64_t *outC = new uint64_t[m * n];   // Output matrix (m x n)

    // Generate test data with fixed seed for reproducibility
    std::mt19937 gen(4);
    std::uniform_int_distribution<uint64_t> mat_dis(1, mask_A);    // Matrix elements
    std::uniform_int_distribution<uint64_t> vec_dis(1, mask_B);    // Vector elements

    cout << "Generating test data for " << m << "x" << n << " matrix and " << m << "-element vector..." << endl;

    // Store original test values for later verification
    uint64_t *test_matrix = new uint64_t[m * n];
    uint64_t *test_vector = new uint64_t[m];
    
    // Generate clear text inputs for crossterm computation
    // Note: In crossterm, inA and inB are NOT shared values, outC is shared
    for (int i = 0; i < m * n; i++) {
        test_matrix[i] = mat_dis(gen);
        inA[i] = test_matrix[i];  // Both parties have the same clear matrix
    }
    
    for (int i = 0; i < m; i++) {
        test_vector[i] = vec_dis(gen);
        inB[i] = test_vector[i];  // Both parties have the same clear vector
    }

    cout << "Matrix A (first 5x5 elements):" << endl;
    for (int i = 0; i < min(5, m); i++) {
        for (int j = 0; j < min(5, n); j++) {
            printf("%4lu ", inA[i * n + j]);
        }
        printf("\n");
    }
    
    cout << "Vector B (first 5 elements): ";
    for (int i = 0; i < min(5, m); i++) {
        printf("%lu ", inB[i]);
    }
    printf("\n");

    cout << "Starting matrix_vector_crossterm computation..." << endl;
    
    // Record communication start
    size_t comm_start = iopack->io->counter;
    auto start_time = chrono::high_resolution_clock::now();

    // Call the matrix_vector_crossterm function
    if (party == sci::ALICE) {
        gp->matrix_vector_crossterm(m, n, inA, nullptr, outC, bwA, bwB);
    } else {
        gp->matrix_vector_crossterm(m, n, nullptr, inB, outC, bwA, bwB);
    }

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
        int correct_results = 0;
        int total_elements = m * n;
        
        cout << "Verifying crossterm computation: matrix[i][j] * vector[i] for all i,j" << endl;
        
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                int idx = i * n + j;
                
                // Reconstruct the actual result from shares
                uint64_t actual_result = (result_alice[idx] + result_bob[idx]) & mask_C;
                
                // Expected result: matrix[i][j] * vector[i] (crossterm)
                uint64_t expected_result = (test_matrix[idx] * test_vector[i]) & mask_C;
                
                // Calculate errors
                double absolute_error = fabs((double)actual_result - (double)expected_result);
                total_error += absolute_error;
                
                if (absolute_error > max_error) max_error = absolute_error;
                
                // Consider result correct if error is 0
                if (actual_result == expected_result) correct_results++;
                
                // Print details for first few results
                if (i < 3 && j < 5) {
                    printf("Element [%d][%d]: mat=%lu, vec=%lu, actual=%lu, expected=%lu, error=%.0f\n", 
                           i, j, test_matrix[idx], test_vector[i], actual_result, expected_result, absolute_error);
                }
            }
        }
        
        cout << "\n=== Summary Statistics ===" << endl;
        printf("Total elements: %d\n", total_elements);
        printf("Correct results: %d (%.2f%%)\n", 
               correct_results, 100.0 * correct_results / total_elements);
        printf("Average absolute error: %.4f\n", total_error / total_elements);
        printf("Maximum absolute error: %.4f\n", max_error);
        
        cout << "\n=== Performance Statistics ===" << endl;
        printf("Communication: %zu bytes\n", comm_end - comm_start);
        printf("Communication per element: %.2f bytes\n", 
               (double)(comm_end - comm_start) / total_elements);
        printf("Computation time: %ld ms\n", duration.count());
        printf("Throughput: %.2f matrix_elements/sec\n", 
               1000.0 * total_elements / duration.count());
        
        cout << "\n=== matrix_vector_crossterm Specific Analysis ===" << endl;
        printf("Matrix dimensions: %d x %d\n", m, n);
        printf("Vector length: %d\n", m);
        printf("Matrix bit width: %d bits\n", bwA);
        printf("Vector bit width: %d bits\n", bwB);
        printf("Output bit width: %d bits\n", bwC);
        printf("Function tested: matrix[i][j] * vector[i] for all i,j (crossterm)\n");
        printf("Total OT rounds: %d (one for each bit of vector elements)\n", bwB);
        printf("Total AES operations: %d matrix rows * %d OT rounds\n", m, bwB);
        
        // Expected communication analysis
        int values_per_block = 128 / bwA;
        int blocks_per_row = (n + values_per_block - 1) / values_per_block;
        int total_blocks = m * blocks_per_row;
        printf("Blocks per matrix row: %d\n", blocks_per_row);
        printf("Total blocks transmitted: %d per round\n", total_blocks);
        printf("Expected comm per round: %d bytes (2 encrypted matrices + %d OT)\n", 
               total_blocks * 16 * 2, m * 16);
    }

    // Cleanup
    delete[] inA;
    delete[] inB;
    delete[] outC;
    delete[] result_alice;
    delete[] result_bob;
    delete[] test_matrix;
    delete[] test_vector;
    
    delete gp;
    delete aux;
    delete prod;
    delete otpack;
    delete iopack;

    cout << "\nmatrix_vector_crossterm test completed successfully!" << endl;
    
    return 0;
}
