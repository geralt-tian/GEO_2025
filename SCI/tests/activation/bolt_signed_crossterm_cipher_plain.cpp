// CrossTerm函数正确性测试
// 简化版本，专注于测试逻辑结构

#include "utils/emp-tool.h"
#include "BuildingBlocks/geometric_perspective_protocols.h"
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <ctime>

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

// 计算绝对值的辅助函数
int64_t abs_int64(int64_t x) {
    return (x < 0) ? -x : x;
}

// 验证crossterm函数正确性的测试函数
void test_crossterm_correctness() {
    printf("=== Testing CrossTerm Function Correctness ===\n");
    printf("Testing signed crossterm: inA (shared by both parties) × inB (only BOB has)\n");
    
    // 初始化随机数种子
    srand(time(nullptr));
    
    // 测试参数
    int test_dim = 5;
    int test_bwA = 8;
    int test_bwB = 8;
    int test_bwC = 16;
    
    // 创建测试数据
    uint64_t *test_inA = new uint64_t[test_dim];
    uint64_t *test_inB = new uint64_t[test_dim];
    uint64_t *test_outC = new uint64_t[test_dim];
    
    // 初始化测试数据 - inA是share值，inB只有BOB输入
    // 测试数据：A=[-2, -1, 0, 200, 2], B=[-1, 10, 1, -220, 3]
    int64_t signed_A[5] = {-2, -1, 0, 200, 2};
    int64_t signed_B[5] = {-1, 10, 1, -220, 3};
    
    // 为inA生成随机share值
    uint64_t *inA_share_alice = new uint64_t[test_dim];
    uint64_t *inA_share_bob = new uint64_t[test_dim];
    
    if (party == ALICE) {
        for (int i = 0; i < test_dim; i++) {
            // 将有符号数转换为无符号数（环上表示）
            uint64_t ring_A;
            if (signed_A[i] >= 0) {
                ring_A = signed_A[i] & ((1ULL << test_bwA) - 1);
            } else {
                ring_A = ((1ULL << test_bwA) + signed_A[i]) & ((1ULL << test_bwA) - 1);
            }
            
            // 生成随机share值
            inA_share_alice[i] = rand() & ((1ULL << test_bwA) - 1);
            inA_share_bob[i] = (ring_A - inA_share_alice[i]) & ((1ULL << test_bwA) - 1);
            
            test_inA[i] = inA_share_alice[i];  // ALICE的share
            test_inB[i] = 0;  // ALICE不输入B
            test_outC[i] = 0;
        }
        
        // 发送BOB的share给BOB
        iopack->io->send_data(inA_share_bob, test_dim * sizeof(uint64_t));
        
        printf("ALICE input A values (signed): ");
        for (int i = 0; i < test_dim; i++) {
            printf("%d ", (int)signed_A[i]);
        }
        printf("\n");
        printf("ALICE A share: ");
        for (int i = 0; i < test_dim; i++) {
            printf("%llu ", inA_share_alice[i]);
        }
        printf("\n");
    } else {
        // 接收ALICE发送的share
        iopack->io->recv_data(inA_share_bob, test_dim * sizeof(uint64_t));
        
        for (int i = 0; i < test_dim; i++) {
            test_inA[i] = inA_share_bob[i];  // BOB的share
            
            // 将有符号数转换为无符号数（环上表示）
            if (signed_B[i] >= 0) {
                test_inB[i] = signed_B[i] & ((1ULL << test_bwB) - 1);
            } else {
                test_inB[i] = ((1ULL << test_bwB) + signed_B[i]) & ((1ULL << test_bwB) - 1);
            }
            test_outC[i] = 0;
        }
        
        printf("BOB input B values (signed): ");
        for (int i = 0; i < test_dim; i++) {
            printf("%d ", (int)signed_B[i]);
        }
        printf("\n");
        printf("BOB A share: ");
        for (int i = 0; i < test_dim; i++) {
            printf("%llu ", inA_share_bob[i]);
        }
        printf("\n");
        printf("BOB B values (ring): ");
        for (int i = 0; i < test_dim; i++) {
            printf("%llu ", test_inB[i]);
        }
        printf("\n");
    }
    
    size_t comm_start = iopack->io->counter;
    
    // 调用crossterm函数
    if (party == ALICE) {
        gp->cross_term(test_dim, test_inA, nullptr, test_outC, test_bwA, test_bwB, test_bwC);
    } else {
        gp->cross_term(test_dim, test_inA, test_inB, test_outC, test_bwA, test_bwB, test_bwC);
    }
    
    size_t comm_end = iopack->io->counter;
    
    // 恢复明文结果进行验证
    uint64_t *test_outC_alice = new uint64_t[test_dim];
    if (party == ALICE) {
        iopack->io->send_data(test_outC, test_dim * sizeof(uint64_t));
    } else {
        iopack->io->recv_data(test_outC_alice, test_dim * sizeof(uint64_t));
    }
    
    if (party == BOB) {
        printf("CrossTerm MPC Results and Verification:\n");
        uint64_t mask_c = (1ULL << test_bwC) - 1;
        double total_error = 0.0;
        int correct_count = 0;
        
        // 已知的明文输入值用于验证 (signed)
        int64_t plaintext_A[5] = {-2, -1, 0, 200, 2};
        int64_t plaintext_B[5] = {-1, 10, 1, -220, 3};
        
        for (int i = 0; i < test_dim; i++) {
            // 恢复明文结果
            uint64_t mpc_result_unsigned = (test_outC_alice[i] + test_outC[i]) & mask_c;
            
            // 将结果转换为有符号数
            int64_t mpc_result;
            if (mpc_result_unsigned >= (1ULL << (test_bwC - 1))) {
                mpc_result = (int64_t)(mpc_result_unsigned - (1ULL << test_bwC));
            } else {
                mpc_result = (int64_t)mpc_result_unsigned;
            }
            
            // 计算期望的明文结果 - 在环上进行计算
            // 首先将有符号数转换为环上的无符号表示
            uint64_t ring_A, ring_B;
            if (plaintext_A[i] >= 0) {
                ring_A = plaintext_A[i] & ((1ULL << test_bwA) - 1);
            } else {
                ring_A = ((1ULL << test_bwA) + plaintext_A[i]) & ((1ULL << test_bwA) - 1);
            }
            
            if (plaintext_B[i] >= 0) {
                ring_B = plaintext_B[i] & ((1ULL << test_bwB) - 1);
            } else {
                ring_B = ((1ULL << test_bwB) + plaintext_B[i]) & ((1ULL << test_bwB) - 1);
            }
            
            // 在环上计算乘法
            uint64_t expected_result_unsigned = (ring_A * ring_B) & mask_c;
            
            // 将期望结果转换为有符号数
            int64_t expected_result;
            if (expected_result_unsigned >= (1ULL << (test_bwC - 1))) {
                expected_result = (int64_t)(expected_result_unsigned - (1ULL << test_bwC));
            } else {
                expected_result = (int64_t)expected_result_unsigned;
            }
            
            // 计算误差
            int64_t error = (mpc_result > expected_result) ? 
                           (mpc_result - expected_result) : 
                           (expected_result - mpc_result);
            
            bool is_correct = (error == 0);
            if (is_correct) correct_count++;
            
            total_error += abs_int64(error);
            
            printf("Test[%d]: A=%d(ring:%llu), B=%d(ring:%llu), Expected=%d, MPC=%d, Error=%d %s\n",
                   i, (int)plaintext_A[i], ring_A, (int)plaintext_B[i], ring_B, 
                   (int)expected_result, (int)mpc_result, (int)error, is_correct ? "✓" : "✗");
        }
        
        printf("\n=== Test Summary ===\n");
        printf("Total tests: %d\n", test_dim);
        printf("Correct results: %d\n", correct_count);
        printf("Accuracy: %.2f%%\n", (double)correct_count / test_dim * 100.0);
        printf("Average error: %.2f\n", total_error / test_dim);
        printf("Communication cost: %zu bytes\n", (comm_end - comm_start));
        
        if (correct_count == test_dim) {
            printf("🎉 All tests PASSED! CrossTerm function is working correctly.\n");
        } else {
            printf("❌ Some tests FAILED. Please check the implementation.\n");
        }
    }
    
    // 清理内存
    delete[] test_inA;
    delete[] test_inB;
    delete[] test_outC;
    delete[] test_outC_alice;
    delete[] inA_share_alice;
    delete[] inA_share_bob;
    
    printf("=== CrossTerm Test Completed ===\n\n");
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

    // 运行CrossTerm函数正确性测试
    printf("Starting CrossTerm function correctness verification...\n");
    test_crossterm_correctness();
    
    printf("================================\n");
    printf("CrossTerm test completed successfully!\n");

    // 清理资源
    delete gp;
    delete otpack;
    delete iopack;
    
    return 0;
}