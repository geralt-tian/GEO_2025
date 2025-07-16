// MW协议测试文件
// 测试MW (Most Significant Wrap) 协议的正确性

#include "utils/emp-tool.h"
#include "BuildingBlocks/geometric_perspective_protocols.h"
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>

using namespace sci;
using namespace std;

// 定义常量
#define ALICE 1
#define BOB 2

int party, port = 32000;
string address = "127.0.0.1";
IOPack *iopack;
OTPack *otpack;
GeometricPerspectiveProtocols *gp;

// MW明文计算函数
void compute_MW_plain(const uint64_t *x0, const uint64_t *x1, uint64_t *MW,
                      size_t len, uint64_t N) {
    for (size_t i = 0; i < len; ++i) {
        uint64_t sum = x0[i] + x1[i];
        if (0 <= sum && sum < N / 2) {
            MW[i] = 0;
        } else if (N / 2 <= sum && sum < 3 * N / 2) {
            MW[i] = 1;
        } else if (3 * N / 2 <= sum && sum < 2 * N) {
            MW[i] = 2;
        } else {
            MW[i] = 3; // 默认无效
        }
    }
}

// 计算绝对值的辅助函数
int64_t abs_int64(int64_t x) {
    return (x < 0) ? -x : x;
}

// 测试基本MW协议的函数
void test_basic_mw() {
    printf("=== Testing Basic MW Protocol ===\n");
    
    // 初始化随机数种子
    srand(time(nullptr));
    
    // 测试参数 - 环大小为2^21
    int test_dim = 10;
    int test_bw = 21;  // l = 21
    uint64_t N = 1ULL << test_bw;  // 2^21 = 2097152
    
    printf("\n--- Testing with l=%d, N=2^%d=%llu ---\n", test_bw, test_bw, N);
    
    // 创建测试数据
    uint64_t *test_input = new uint64_t[test_dim];
    uint64_t *test_output = new uint64_t[test_dim];
    
    // 初始化测试数据 - 根据环大小调整输入值
    for (int i = 0; i < test_dim; i++) {
        // 生成覆盖不同环区间的测试数据
        if (i == 0) test_input[i] = N/8;         // 第一个区间
        else if (i == 1) test_input[i] = N/4;    // 边界
        else if (i == 2) test_input[i] = N/2;    // 中点
        else if (i == 3) test_input[i] = 3*N/4;  // 第三个区间
        else if (i == 4) test_input[i] = 7*N/8;  // 接近边界
        else test_input[i] = (N-1) - (i-5);     // 其他值
        
        test_output[i] = 0;
    }
    
    printf("Test input values (ring): ");
    for (int i = 0; i < test_dim; i++) {
        printf("%llu ", test_input[i]);
    }
    printf("\n");
    
    size_t comm_start = iopack->io->counter;
    
    // 调用MW协议
    gp->mw(test_dim, test_input, test_output, test_bw, 2);
    
    size_t comm_end = iopack->io->counter;
    
    // 恢复明文结果进行验证
    uint64_t *test_output_alice = new uint64_t[test_dim];
    if (party == ALICE) {
        iopack->io->send_data(test_output, test_dim * sizeof(uint64_t));
    } else {
        iopack->io->recv_data(test_output_alice, test_dim * sizeof(uint64_t));
    }
    
    if (party == BOB) {
        printf("MW Protocol Results and Verification:\n");
        int correct_count = 0;
        
        for (int i = 0; i < test_dim; i++) {
            // 恢复明文结果
            uint64_t mpc_result = (test_output_alice[i] + test_output[i]) & 3; // MW输出是2位
            
            // 计算期望的明文结果
            // MW是基于输入值在环中的位置计算的
            uint64_t expected_result;
            uint64_t input_val = test_input[i];
            
            if (input_val < N/4) {
                expected_result = 0;  // 在[0, N/4)区间
            } else if (input_val < 3*N/4) {
                expected_result = 1;  // 在[N/4, 3*N/4)区间  
            } else {
                expected_result = 2;  // 在[3*N/4, N)区间
            }
            
            bool is_correct = (mpc_result == expected_result);
            if (is_correct) correct_count++;
            
            printf("Test[%d]: Input=%llu, Expected=%llu, MPC=%llu %s\n",
                   i, input_val, expected_result, mpc_result,
                   is_correct ? "✓" : "✗");
        }
        
        printf("\n=== MW Test Summary (l=%d) ===\n", test_bw);
        printf("Total tests: %d\n", test_dim);
        printf("Correct results: %d\n", correct_count);
        printf("Accuracy: %.2f%%\n", (double)correct_count / test_dim * 100.0);
        printf("Communication cost: %zu bytes\n", (comm_end - comm_start));
        
        if (correct_count == test_dim) {
            printf("🎉 All MW tests PASSED for l=%d!\n", test_bw);
        } else {
            printf("❌ Some MW tests FAILED for l=%d.\n", test_bw);
        }
    }
    
    // 清理内存
    delete[] test_input;
    delete[] test_output;
    delete[] test_output_alice;
    
    printf("=== Basic MW Test Completed ===\n\n");
}

// 测试带参数B的MW协议
void test_mw_with_B() {
    printf("=== Testing MW with B Protocol ===\n");
    
    // 测试参数 - 环大小为2^21，B = 5/8 * N
    int test_dim = 10;
    int test_bw = 21;  // l = 21
    uint64_t N = 1ULL << test_bw;  // 2^21 = 2097152
    uint64_t B = (5 * N) / 8;  // B = 5/8 * N = 1310720
    
    printf("\n--- Testing MW with B, l=%d, N=2^%d=%llu, B=%llu (5/8*N) ---\n", test_bw, test_bw, N, B);
    
    // 创建测试数据
    uint64_t *test_input = new uint64_t[test_dim];
    uint64_t *test_output = new uint64_t[test_dim];
    
    // 初始化测试数据 - 根据环大小和B值调整输入值
    for (int i = 0; i < test_dim; i++) {
        // 生成围绕B值的测试数据
        if (i == 0) test_input[i] = B/2;        // 小于B
        else if (i == 1) test_input[i] = B-1;   // 刚好小于B
        else if (i == 2) test_input[i] = B;     // 等于B
        else if (i == 3) test_input[i] = B+1;   // 刚好大于B
        else if (i == 4) test_input[i] = B*2;   // 大于B
        else test_input[i] = (B*3) % N;         // 其他值
        
        test_output[i] = 0;
    }
    
    printf("Test input values (ring): ");
    for (int i = 0; i < test_dim; i++) {
        printf("%llu ", test_input[i]);
    }
    printf("\n");
    
    size_t comm_start = iopack->io->counter;
    
    // 调用MW with B协议
    gp->mwwithB(test_dim, B, test_input, test_output, test_bw, 2);
    
    size_t comm_end = iopack->io->counter;
    
    // 恢复明文结果进行验证
    uint64_t *test_output_alice = new uint64_t[test_dim];
    if (party == ALICE) {
        iopack->io->send_data(test_output, test_dim * sizeof(uint64_t));
    } else {
        iopack->io->recv_data(test_output_alice, test_dim * sizeof(uint64_t));
    }
    
    if (party == BOB) {
        printf("MW with B Protocol Results and Verification:\n");
        int correct_count = 0;
        
        for (int i = 0; i < test_dim; i++) {
            // 恢复明文结果
            uint64_t mpc_result = (test_output_alice[i] + test_output[i]) & ((1ULL << test_bw) - 1);
            
            // 计算期望的明文结果
            // MWwithB的逻辑基于输入值与B的关系
            uint64_t expected_result;
            uint64_t input_val = test_input[i];
            
            // 简化的期望结果计算（根据MW with B的具体逻辑调整）
            if (input_val >= B) {
                expected_result = 1;
            } else {
                expected_result = 0;  
            }
            
            bool is_correct = (mpc_result == expected_result);
            if (is_correct) correct_count++;
            
            printf("Test[%d]: Input=%llu, B=%llu, Expected=%llu, MPC=%llu %s\n",
                   i, input_val, B, expected_result, mpc_result,
                   is_correct ? "✓" : "✗");
        }
        
        printf("\n=== MW with B Test Summary (l=%d) ===\n", test_bw);
        printf("Total tests: %d\n", test_dim);
        printf("Correct results: %d\n", correct_count);
        printf("Accuracy: %.2f%%\n", (double)correct_count / test_dim * 100.0);
        printf("Communication cost: %zu bytes\n", (comm_end - comm_start));
        
        if (correct_count == test_dim) {
            printf("🎉 All MW with B tests PASSED for l=%d!\n", test_bw);
        } else {
            printf("❌ Some MW with B tests FAILED for l=%d.\n", test_bw);
        }
    }
    
    // 清理内存
    delete[] test_input;
    delete[] test_output;
    delete[] test_output_alice;
    
    printf("=== MW with B Test Completed ===\n\n");
}

// 测试MW明文计算的正确性
void test_mw_plaintext() {
    printf("=== Testing MW Plaintext Computation ===\n");
    
    // 测试环大小为2^21
    int test_dim = 10;
    int test_bw = 21;  // l = 21
    uint64_t N = 1ULL << test_bw;  // 2^21 = 2097152
    
    printf("\n--- Testing Plaintext MW with l=%d, N=2^%d=%llu ---\n", test_bw, test_bw, N);
    
    uint64_t *x0 = new uint64_t[test_dim];
    uint64_t *x1 = new uint64_t[test_dim];
    uint64_t *MW_result = new uint64_t[test_dim];
    
    // 创建测试数据覆盖不同的sum范围
    for (int i = 0; i < test_dim; i++) {
        // 根据环大小调整测试数据
        x0[i] = (i * N/test_dim) % N;  // 均匀分布在环上
        x1[i] = (i * N/test_dim + N/8) % N;  // 偏移一些
    }
    
    // 计算MW明文结果
    compute_MW_plain(x0, x1, MW_result, test_dim, N);
    
    printf("MW Plaintext Test Results:\n");
    printf("N = %llu, N/2 = %llu, 3*N/2 = %llu\n", N, N/2, 3*N/2);
    
    int correct_count = 0;
    for (int i = 0; i < test_dim; i++) {
        uint64_t sum = x0[i] + x1[i];
        printf("Test[%d]: x0=%llu, x1=%llu, sum=%llu, MW=%llu\n",
               i, x0[i], x1[i], sum, MW_result[i]);
        
        // 验证逻辑
        uint64_t expected;
        if (sum < N/2) expected = 0;
        else if (sum < 3*N/2) expected = 1;
        else if (sum < 2*N) expected = 2;
        else expected = 3;
        
        if (MW_result[i] != expected) {
            printf("  ❌ Error: Expected %llu but got %llu\n", expected, MW_result[i]);
        } else {
            printf("  ✓ Correct\n");
            correct_count++;
        }
    }
    
    printf("=== MW Plaintext Test Summary (l=%d) ===\n", test_bw);
    printf("Correct results: %d/%d\n", correct_count, test_dim);
    printf("Accuracy: %.2f%%\n", (double)correct_count / test_dim * 100.0);
    
    delete[] x0;
    delete[] x1;
    delete[] MW_result;
    
    printf("=== MW Plaintext Test Completed ===\n\n");
}

int main(int argc, char **argv) {
    ArgMapping amap;

    amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
    amap.arg("p", port, "Port Number");
    amap.arg("ip", address, "IP Address of server (ALICE)");
    amap.parse(argc, argv);

    // 初始化通信和协议
    iopack = new IOPack(party, port, "127.0.0.1");
    otpack = new OTPack(iopack, party);
    gp = new GeometricPerspectiveProtocols(party, iopack, otpack);

    // 运行MW协议测试
    printf("Starting MW Protocol Testing...\n");
    printf("Party: %s\n", (party == ALICE) ? "ALICE" : "BOB");
    
    // 先测试明文计算
    if (party == ALICE) {
        test_mw_plaintext();
    }
    
    // 测试基本MW协议
    test_basic_mw();
    
    // 测试带参数B的MW协议
    test_mw_with_B();
    
    printf("================================\n");
    printf("MW Protocol testing completed!\n");

    // 清理资源
    delete gp;
    delete otpack;
    delete iopack;
    
    return 0;
}
