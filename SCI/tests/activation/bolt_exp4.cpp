// Test file for exp4 function
// This file tests the exp4 function which computes exp(x) using polynomial approximation
// Including ULP error analysis, communication measurement, and timing

#include "LinearOT/linear-ot.h"
#include "utils/emp-tool.h"
#include <cstdint>
#include <cstdio>
#include <iostream>

#include "FloatingPoint/floating-point.h"
#include "FloatingPoint/fp-math.h"
#include <limits>
#include <random>
#include "BuildingBlocks/aux-protocols.h"
#include "BuildingBlocks/geometric_perspective_protocols.h"
#include "BuildingBlocks/truncation.h"
#include "Math/math-functions.h"
#include "Millionaire/equality.h"
#include "Millionaire/millionaire.h"
#include "Millionaire/millionaire_with_equality.h"
#include <chrono>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>

using namespace sci;
using namespace std;

// Global parameters
int party, port = 32000;
string address = "127.0.0.1";
IOPack *iopack;
OTPack *otpack;
LinearOT *prod;
AuxProtocols *aux;
FPMath *fpmath;

// Test parameters - configurable
int dim = 1024; 
// int dim = 1048576/4;              // Number of test elements
int32_t ell = 37;           // Input bit width
int32_t scale = 12;         // Input fractional bits
double test_range_min = -10.0;  // Minimum test value
double test_range_max = 0;   // Maximum test value
bool verbose = false;        // Verbose output
double ulp_threshold = 500.0;  // ULP error threshold for reporting high error cases

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
    // Parse command line arguments
    ArgMapping amap;
    amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
    amap.arg("p", port, "Port Number");
    amap.arg("ip", address, "IP Address of server (ALICE)");
    amap.arg("d", dim, "Number of elements to test");
    amap.arg("ell", ell, "Input bit width");
    amap.arg("s", scale, "Input fractional bits");
    amap.arg("min", test_range_min, "Minimum test value");
    amap.arg("max", test_range_max, "Maximum test value");
    amap.arg("v", verbose, "Verbose output");
    amap.arg("ulp_thresh", ulp_threshold, "ULP error threshold for reporting high error cases");
    amap.parse(argc, argv);

    cout << "=== exp4 Test Configuration ===" << endl;
    cout << "Party: " << (party == ALICE ? "ALICE" : "BOB") << endl;
    cout << "Dimensions: " << dim << endl;
    cout << "Input format: " << ell << "." << scale << endl;
    cout << "Test range: [" << test_range_min << ", " << test_range_max << "]" << endl;
    cout << "ULP error threshold: " << ulp_threshold << endl;
    cout << "================================" << endl;

    // Initialize communication and protocols
    iopack = new IOPack(party, port, address);
    otpack = new OTPack(iopack, party);
    prod = new LinearOT(party, iopack, otpack);
    aux = new AuxProtocols(party, iopack, otpack);
    fpmath = new FPMath(party, iopack, otpack);

    // Allocate arrays
    uint64_t *inA = new uint64_t[dim];
    uint64_t *result = new uint64_t[dim];
    
    // Generate random test data
    mt19937 gen(42);  // Fixed seed for reproducibility
    uniform_real_distribution<double> dis(test_range_min, test_range_max);
    
    cout << "Generating test data for " << dim << " elements..." << endl;
    
    // Store original test values for later verification
    double *test_values = new double[dim];
    
    // Generate secret shared inputs (both parties generate same test values)
    for (int i = 0; i < dim; i++) {
        double test_val = dis(gen);
        test_values[i] = test_val;  // Both parties store the same original value
        
        if (party == ALICE) {
            // ALICE gets 70% of the value as her share
            double alice_share = test_val * 0.7;
            inA[i] = double2fix(alice_share, scale, ell);
        } else {
            // BOB gets 30% of the value as his share  
            double bob_share = test_val * 0.3;
            inA[i] = double2fix(bob_share, scale, ell);
        }
    }
    
    // Create FixArray for input
    FixArray input = fpmath->fix->input(ALICE, dim, inA, true, ell, scale);
    
    cout << "Starting exp4 computation..." << endl;
    
    // Measure computation time and communication
    auto start_time = chrono::high_resolution_clock::now();
    size_t comm_start = iopack->get_comm();
    
    // Call exp4 function
    auto result_tuple = fpmath->exp4(input);
    FixArray result_fix = get<0>(result_tuple);
    FixArray l_short_fix = get<1>(result_tuple);
    
    size_t comm_end = iopack->get_comm();
    auto end_time = chrono::high_resolution_clock::now();
    
    // Copy results
    memcpy(result, result_fix.data, dim * sizeof(uint64_t));
    
    // Calculate timing
    auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
    size_t comm_bytes = comm_end - comm_start;
    
    cout << "Computation completed in " << duration.count() << " ms" << endl;
    
    // Exchange data for error analysis
    uint64_t *result_alice = new uint64_t[dim];
    uint64_t *result_bob = new uint64_t[dim];
    uint64_t *input_alice = new uint64_t[dim];
    uint64_t *input_bob = new uint64_t[dim];
    
    if (party == ALICE) {
        iopack->io->send_data(result, dim * sizeof(uint64_t));
        iopack->io->send_data(inA, dim * sizeof(uint64_t));
        iopack->io->recv_data(result_bob, dim * sizeof(uint64_t));
        iopack->io->recv_data(input_bob, dim * sizeof(uint64_t));
        memcpy(result_alice, result, dim * sizeof(uint64_t));
        memcpy(input_alice, inA, dim * sizeof(uint64_t));
        printf("Communication: %zu bytes\n", comm_bytes);
    } else {
        iopack->io->recv_data(result_alice, dim * sizeof(uint64_t));
        iopack->io->recv_data(input_alice, dim * sizeof(uint64_t));
        iopack->io->send_data(result, dim * sizeof(uint64_t));
        iopack->io->send_data(inA, dim * sizeof(uint64_t));
        memcpy(result_bob, result, dim * sizeof(uint64_t));
        memcpy(input_bob, inA, dim * sizeof(uint64_t));
    }
    
    // Accuracy analysis (only on BOB to avoid duplicate output)
    if (party == BOB) {
        cout << "\n=== Accuracy Analysis ===" << endl;
        
        double ulp = 1.0 / (1ULL << scale);
        double total_ulp_error = 0.0;
        double max_ulp_error = 0.0;
        double total_relative_error = 0.0;
        double max_relative_error = 0.0;
        int correct_results = 0;
        int high_ulp_count = 0;
        
        for (int i = 0; i < dim; i++) {
            // Reconstruct actual result from shares
            double actual_result = fix2double(result_alice[i], result_bob[i], ell, scale);
            
            // Use original test value (more accurate than reconstructing from shares)
            double input_val = test_values[i];
            
            // Expected result: exp(input_val)
            double expected_result = exp(input_val);
            
            // Calculate errors
            double absolute_error = fabs(actual_result - expected_result);
            double ulp_error = absolute_error / ulp;
            double relative_error = (expected_result != 0) ? absolute_error / fabs(expected_result) : 0;
            
            total_ulp_error += ulp_error;
            total_relative_error += relative_error;
            
            if (ulp_error > max_ulp_error) max_ulp_error = ulp_error;
            if (relative_error > max_relative_error) max_relative_error = relative_error;
            
            // Consider result correct if ULP error < 10.0
            if (ulp_error < 10.0) correct_results++;
            
            // Count and print high ULP error cases (> threshold)
            if (ulp_error > ulp_threshold) {
                high_ulp_count++;
                printf("HIGH ULP ERROR [%d]: input=%.6f, actual=%.6f, expected=%.6f, ULP_error=%.2f\n", 
                       i, input_val, actual_result, expected_result, ulp_error);
            }
            
            // Print details for verbose mode or first few results
            if (verbose && i < 5) {
                printf("Test %d: input=%.6f, actual=%.6f, expected=%.6f, ULP_error=%.2f\n", 
                       i, input_val, actual_result, expected_result, ulp_error);
            }
        }
        
        cout << "\n=== Summary Statistics ===" << endl;
        printf("Total tests: %d\n", dim);
        printf("Correct results (ULP < 10): %d (%.2f%%)\n", 
               correct_results, 100.0 * correct_results / dim);
        printf("High ULP error cases (ULP > 500): %d (%.2f%%)\n", 
               high_ulp_count, 100.0 * high_ulp_count / dim);
        printf("Average ULP error: %.4f\n", total_ulp_error / dim);
        printf("Maximum ULP error: %.4f\n", max_ulp_error);
        printf("Average relative error: %.6f%%\n", 100.0 * total_relative_error / dim);
        printf("Maximum relative error: %.6f%%\n", 100.0 * max_relative_error);
        
        cout << "\n=== Performance Statistics ===" << endl;
        printf("Communication: %zu bytes\n", comm_bytes);
        printf("Communication per element: %.2f bytes\n", (double)comm_bytes / dim);
        printf("Computation time: %ld ms\n", duration.count());
        printf("Throughput: %.2f exp4/sec\n", 1000.0 * dim / duration.count());
        
        cout << "\n=== exp4 Function Analysis ===" << endl;
        printf("Function tested: exp4(x) - polynomial approximation of exp(x)\n");
        printf("Input range tested: [%.2f, %.2f]\n", test_range_min, test_range_max);
        printf("Secret sharing: ALICE=70%%, BOB=30%% of each input value\n");
    }
    
    // Cleanup
    delete[] inA;
    delete[] result;
    delete[] result_alice;
    delete[] result_bob;
    delete[] input_alice;
    delete[] input_bob;
    delete[] test_values;
    
    delete fpmath;
    delete aux;
    delete prod;
    delete otpack;
    delete iopack;

    cout << "\nexp4 test completed successfully!" << endl;
    
    return 0;
}