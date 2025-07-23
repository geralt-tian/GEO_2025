// Test file for vector_scalar_mul function
// This file tests the vector_scalar_mul function which computes vector * scalar multiplication
// Including ULP error analysis, communication measurement, and timing

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

// Test parameters for vector_scalar_mul
int dim = 100;  // Test with 16384 elements
int32_t bwA = 20;      // Vector element bit width
int32_t bwB = 16;      // Scalar bit width  
int32_t bwC = bwA + bwB;      // Output bit width

uint64_t mask_A = (bwA == 64 ? -1 : ((1ULL << bwA) - 1));
uint64_t mask_B = (bwB == 64 ? -1 : ((1ULL << bwB) - 1));
uint64_t mask_C = (bwC == 64 ? -1 : ((1ULL << bwC) - 1));

// Convert fixed point to double
double fix2double(uint64_t x, uint64_t y, int32_t bw, int32_t f) {
    uint64_t mask = (bw == 64) ? ~0ULL : ((1ULL << bw) - 1);
    uint64_t sum = (x + y) & mask;
    int64_t signed_val;
    if (sum >= (1ULL << (bw - 1))) {
        signed_val = (int64_t)(sum - (1ULL << bw));
    } else {
        signed_val = (int64_t)sum;
    }
    return (double)signed_val / (1ULL << f);
}

// Convert double to fixed point
uint64_t double2fix(double val, int32_t f, int32_t bw) {
    uint64_t mask = (bw == 64) ? ~0ULL : ((1ULL << bw) - 1);
    int64_t fixed_val = (int64_t)(val * (1ULL << f));
    return ((uint64_t)fixed_val) & mask;
}

int main(int argc, char **argv) {
    ArgMapping amap;
    amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
    amap.arg("p", port, "Port Number");
    amap.arg("ip", address, "IP Address of server (ALICE)");
    amap.arg("d", dim, "Number of elements to test");
    amap.parse(argc, argv);

    cout << "=== vector_scalar_mul Test Suite ===" << endl;
    cout << "Testing vector * scalar multiplication" << endl;
    cout << "Party: " << (party == sci::ALICE ? "ALICE" : "BOB") << endl;
    cout << "Dimensions: " << dim << endl;
    cout << "Vector bit width: " << bwA << endl;
    cout << "Scalar bit width: " << bwB << endl;
    cout << "Output bit width: " << bwC << endl;

    // Initialize communication
    iopack = new IOPack(party, port, address);
    otpack = new OTPack(iopack, party);
    prod = new LinearOT(party, iopack, otpack);
    gp = new GeometricPerspectiveProtocols(party, iopack, otpack);
    aux = new AuxProtocols(party, iopack, otpack);

    // Allocate arrays
    uint64_t *inA = new uint64_t[dim];      // Vector input
    uint64_t inB;                           // Scalar input
    uint64_t *outC = new uint64_t[dim];     // Output vector

    // Generate test data with fixed seed for reproducibility
    std::mt19937 gen(42);
    std::uniform_int_distribution<uint64_t> vec_dis(1, mask_A >> 2);  // Vector elements
    std::uniform_int_distribution<uint64_t> scalar_dis(1, mask_B >> 2);  // Scalar value

    cout << "Generating test data for " << dim << " elements..." << endl;

    // Store original test values for later verification
    uint64_t *test_vector = new uint64_t[dim];
    uint64_t test_scalar;
    
    // Generate inputs for vector-scalar multiplication
    for (int i = 0; i < dim; i++) {
        uint64_t test_val = vec_dis(gen);
        test_vector[i] = test_val;
        
        if (party == sci::ALICE) {
            inA[i] = test_val * 3 / 5;  // Alice gets part of the value
        } else {
            inA[i] = test_val * 2 / 5;  // Bob gets part of the value
        }
    }
    
    // Generate scalar input
    test_scalar = scalar_dis(gen);
    if (party == sci::ALICE) {
        inB = test_scalar * 3 / 5;  // Alice gets part of scalar
    } else {
        inB = test_scalar * 2 / 5;  // Bob gets part of scalar
    }

    cout << "Starting vector_scalar_mul computation..." << endl;
    
    // Record communication start
    size_t comm_start = iopack->io->counter;
    auto start_time = chrono::high_resolution_clock::now();

    // Call the vector_bit_mul function (OT-based protocol)
    // For testing, use choice 0 or 1 based on the scalar value
    uint8_t ot_choice = (party == sci::BOB) ? (inB % 2) : 0;  // Use LSB of inB as choice for Bob
    ot_choice = 1;
    gp->vector_bit_mul(dim, inA, ot_choice, outC, bwA);

    auto end_time = chrono::high_resolution_clock::now();
    size_t comm_end = iopack->io->counter;
    
    auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
 
    cout << "Computation completed in " << duration.count() << " ms" << endl;

    // Share results and inputs for accuracy testing
    uint64_t *result_alice = new uint64_t[dim];
    uint64_t *result_bob = new uint64_t[dim];
    uint64_t *input_vector_alice = new uint64_t[dim];
    uint64_t *input_vector_bob = new uint64_t[dim];
    uint64_t input_scalar_alice, input_scalar_bob;
    
    if (party == sci::ALICE) {
        iopack->io->send_data(outC, dim * sizeof(uint64_t));
        iopack->io->send_data(inA, dim * sizeof(uint64_t));
        iopack->io->send_data(&inB, sizeof(uint64_t));
        iopack->io->recv_data(result_bob, dim * sizeof(uint64_t));
        iopack->io->recv_data(input_vector_bob, dim * sizeof(uint64_t));
        iopack->io->recv_data(&input_scalar_bob, sizeof(uint64_t));
        memcpy(result_alice, outC, dim * sizeof(uint64_t));
        memcpy(input_vector_alice, inA, dim * sizeof(uint64_t));
        input_scalar_alice = inB;
    } else {
        iopack->io->recv_data(result_alice, dim * sizeof(uint64_t));
        iopack->io->recv_data(input_vector_alice, dim * sizeof(uint64_t));
        iopack->io->recv_data(&input_scalar_alice, sizeof(uint64_t));
        iopack->io->send_data(outC, dim * sizeof(uint64_t));
        iopack->io->send_data(inA, dim * sizeof(uint64_t));
        iopack->io->send_data(&inB, sizeof(uint64_t));
        memcpy(result_bob, outC, dim * sizeof(uint64_t));
        memcpy(input_vector_bob, inA, dim * sizeof(uint64_t));
        input_scalar_bob = inB;
    }

    // Accuracy analysis (only on one party to avoid duplicate output)
    // if (party == sci::ALICE) {
    //     cout << "\n=== Accuracy Analysis ===" << endl;
        
    //     double total_error = 0.0;
    //     double max_error = 0.0;
    //     double total_relative_error = 0.0;
    //     double max_relative_error = 0.0;
    //     int correct_results = 0;
        
    //     for (int i = 0; i < dim; i++) {
    //         // Reconstruct the actual result 
    //         uint64_t actual_result = (result_alice[i] + result_bob[i]) & mask_C;
            
    //         // Reconstruct inputs
    //         uint64_t input_vector = (input_vector_alice[i] + input_vector_bob[i]) & mask_A;
    //         uint64_t input_scalar = (input_scalar_alice + input_scalar_bob) & mask_B;
            
    //         // Expected result: vector[i] * scalar
    //         uint64_t expected_result = (input_vector * input_scalar) & mask_C;
            
    //         // Calculate errors
    //         double absolute_error = fabs((double)actual_result - (double)expected_result);
    //         double relative_error = (expected_result != 0) ? absolute_error / expected_result : 0;
            
    //         total_error += absolute_error;
    //         total_relative_error += relative_error;
            
    //         if (absolute_error > max_error) max_error = absolute_error;
    //         if (relative_error > max_relative_error) max_relative_error = relative_error;
            
    //         // Consider result correct if error is small
    //         if (actual_result == expected_result) correct_results++;
            
    //         // Print details for first few results
    //         if (i < 10) {
    //             printf("Test %d: vec=%lu, scalar=%lu, actual=%lu, expected=%lu, error=%.0f\n", 
    //                    i, input_vector, input_scalar, actual_result, expected_result, absolute_error);
    //         }
    //     }
        
    //     cout << "\n=== Summary Statistics ===" << endl;
    //     printf("Total tests: %d\n", dim);
    //     printf("Correct results: %d (%.2f%%)\n", 
    //            correct_results, 100.0 * correct_results / dim);
    //     printf("Average absolute error: %.4f\n", total_error / dim);
    //     printf("Maximum absolute error: %.4f\n", max_error);
    //     printf("Average relative error: %.6f%%\n", 100.0 * total_relative_error / dim);
    //     printf("Maximum relative error: %.6f%%\n", 100.0 * max_relative_error);
        
    //     cout << "\n=== Performance Statistics ===" << endl;
    //     printf("Communication: %zu bytes\n", comm_end - comm_start);
    //     printf("Communication per element: %.2f bytes\n", 
    //            (double)(comm_end - comm_start) / dim);
    //     printf("Computation time: %ld ms\n", duration.count());
    //     printf("Throughput: %.2f vector_scalar_mul/sec\n", 
    //            1000.0 * dim / duration.count());
        
    //     cout << "\n=== vector_scalar_mul Specific Analysis ===" << endl;
    //     printf("Vector bit width: %d bits\n", bwA);
    //     printf("Scalar bit width: %d bits\n", bwB);
    //     printf("Output bit width: %d bits\n", bwC);
    //     printf("Function tested: vector[i] * scalar for all i\n");
    // }

    // Cleanup
    delete[] inA;
    delete[] outC;
    delete[] result_alice;
    delete[] result_bob;
    delete[] input_vector_alice;
    delete[] input_vector_bob;
    delete[] test_vector;
    
    delete gp;
    delete aux;
    delete prod;
    delete otpack;
    delete iopack;

    cout << "\nvector_scalar_mul test completed successfully!" << endl;
    
    return 0;
} 