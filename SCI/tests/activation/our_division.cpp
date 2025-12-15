#include "utils/emp-tool.h"
#include "BuildingBlocks/geometric_perspective_protocols.h"
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <chrono>
#include <algorithm>

using namespace sci;
using namespace std;

#define ALICE 1
#define BOB 2

int party, port = 32000;
string address = "127.0.0.1";
int test_dim = 1ULL << 18;
int test_bw = 37;
int precision_bits = 12;
uint64_t test_divisor = 7;

uint64_t max_input_val = (1ULL << 37) / 6;
uint64_t B = (1ULL << (test_bw - 1)) * 0.5;  // B = L/4

IOPack *iopack;
OTPack *otpack;
GeometricPerspectiveProtocols *gp;

void test_division_protocol() {
    uint64_t N = 1ULL << test_bw;
    uint64_t mask_bw = (test_bw == 64 ? -1 : ((1ULL << test_bw) - 1));

    uint64_t *input_array = new uint64_t[test_dim];
    uint64_t *output_array = new uint64_t[test_dim];
    uint64_t *total_x = new uint64_t[test_dim];

    srand(time(nullptr));

    for (int i = 0; i < test_dim; i++) {
        uint64_t r1 = rand();
        uint64_t r2 = rand();
        total_x[i] = ((r1 << 16) ^ r2) % max_input_val;
    }

    if (party == ALICE) {
        for (int i = 0; i < test_dim; i++) {
            if (total_x[i] > 0) {
                uint64_t r = ((uint64_t)rand() << 32) | rand();
                input_array[i] = r % total_x[i];
            } else {
                input_array[i] = 0;
            }
            output_array[i] = 0;
        }
        iopack->io->send_data(input_array, test_dim * sizeof(uint64_t));
        iopack->io->send_data(total_x, test_dim * sizeof(uint64_t));
    } else {
        uint64_t *alice_share = new uint64_t[test_dim];
        iopack->io->recv_data(alice_share, test_dim * sizeof(uint64_t));
        iopack->io->recv_data(total_x, test_dim * sizeof(uint64_t));

        for (int i = 0; i < test_dim; i++) {
            input_array[i] = (total_x[i] - alice_share[i]) & mask_bw;
            output_array[i] = 0;
        }
        delete[] alice_share;
    }

    delete[] total_x;

    size_t comm_start = iopack->io->counter;
    auto start_time = std::chrono::high_resolution_clock::now();

    gp->division(test_dim, input_array, output_array, test_divisor, test_bw, B);

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration_us = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

    size_t comm_end = iopack->io->counter;
    size_t comm_bytes = comm_end - comm_start;

    printf("Division Protocol Test (Party %d)\n", party);
    printf("Parameters: n=%d, bw=%d, precision=%d, divisor=%lu\n",
           test_dim, test_bw, precision_bits, (unsigned long)test_divisor);
    printf("Time: %lld us\n", (long long)duration_us.count());
    printf("Communication: %zu bytes\n", comm_bytes);
    printf("Per division: %.2f bytes\n", (double)comm_bytes / test_dim);

    delete[] input_array;
    delete[] output_array;
}

int main(int argc, char **argv) {
    ArgMapping amap;

    amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
    amap.arg("p", port, "Port Number");
    amap.arg("ip", address, "IP Address of server (ALICE)");
    amap.arg("n", test_dim, "Number of test cases");
    amap.arg("bw", test_bw, "Bit width (default: 37)");
    amap.arg("precision", precision_bits, "Precision bits (default: 12)");
    amap.arg("d", test_divisor, "Divisor for testing");
    amap.parse(argc, argv);

    if (test_divisor == 0) {
        printf("Error: divisor cannot be zero\n");
        return -1;
    }

    iopack = new IOPack(party, port, address);
    otpack = new OTPack(iopack, party);
    gp = new GeometricPerspectiveProtocols(party, iopack, otpack);

    test_division_protocol();

    delete gp;
    delete otpack;
    delete iopack;

    return 0;
}
