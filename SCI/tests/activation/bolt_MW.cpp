// MW协议测试文件 - 重构版
// 专门测试MW (Most Significant Wrap) 协议
// 
// 测试流程：
// 1. 定义真实输入值（包括用户指定的28039, 21038）
// 2. 将真实输入分解为Alice和Bob的秘密分享
// 3. 各方使用自己的share调用gp->mw()协议
// 4. 将MPC结果与基于真实输入的明文MW结果对比
// 5. 验证协议正确性，特别关注bw=15的情况

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

// MW明文计算函数 - 根据MW协议的定义
uint64_t compute_MW_plain(uint64_t input, int bw) {
    uint64_t N = 1ULL << bw;  // 2^bw
    
    // 确保输入在有效范围内
    // input = input % N;
    
    // MW协议的逻辑：检查输入在环中的位置
    // MW检查是否在上半环 [N/2, N)
    if (input <= N/2) {
        return 0;  // 在上半环
    } else if(input <= 3*N/2){
        return 1;
    } else if(input <= 2*N){
        return 2;
    } else{
        return 3;  // 在下半环 [0, N/2)
    }
}

// 测试MW协议的主函数
void test_mw_protocol() {
    printf("=== Testing MW Protocol with Share Inputs ===\n");
    
    // 测试参数
    int test_bw = 15;  // 用户指定的位宽
    uint64_t N = 1ULL << test_bw;  // 2^15 = 32768
    
    printf("\n--- Testing with bw=%d, N=2^%d=%llu ---\n", test_bw, test_bw, N);
    
    // 测试数据维度
    int test_dim = 8;
    uint64_t *input_array = new uint64_t[test_dim];
    uint64_t *output_array = new uint64_t[test_dim];
    
    // 初始化随机数种子
    srand(time(nullptr) + party);
    
    if (party == ALICE) {
        // Alice设置自己的share值
        input_array[0] = 28039;      // 用户指定输入1
        input_array[1] = 21038;      // 用户指定输入2
        input_array[2] = 0;          // 边界值
        input_array[3] = N/4;        // 测试值
        input_array[4] = N/2;        // 中点
        input_array[5] = 3*N/4;      // 测试值
        input_array[6] = N-1;        // 最大值
        input_array[7] = 12345;      // 随机值
        
        printf("Alice shares: ");
        for (int i = 0; i < test_dim; i++) {
            printf("%llu ", input_array[i]);
            output_array[i] = 0;
        }
        printf("\n");
        
        // 发送Alice的shares给Bob用于计算明文期望结果
        iopack->io->send_data(input_array, test_dim * sizeof(uint64_t));
    } else {
        // Bob设置自己的share值
        input_array[0] = 21038;       // 对应Alice的28039
        input_array[1] = 2000;       // 对应Alice的21038  
        input_array[2] = 100;        // 对应Alice的0
        input_array[3] = 500;        // 对应Alice的N/4
        input_array[4] = 1500;       // 对应Alice的N/2
        input_array[5] = 2500;       // 对应Alice的3*N/4
        input_array[6] = 3000;       // 对应Alice的N-1
        input_array[7] = 6789;       // 对应Alice的12345
        
        printf("Bob shares: ");
        for (int i = 0; i < test_dim; i++) {
            printf("%llu ", input_array[i]);
            output_array[i] = 0;
        }
        printf("\n");
    }
    
    // 计算明文期望结果
    uint64_t *alice_shares = new uint64_t[test_dim];
    uint64_t *total_inputs = new uint64_t[test_dim];
    uint64_t *expected_results = new uint64_t[test_dim];
    
    if (party == BOB) {
        // Bob接收Alice的shares
        iopack->io->recv_data(alice_shares, test_dim * sizeof(uint64_t));
        
        // 计算总输入 = Alice_share + Bob_share
        printf("Total inputs (Alice_share + Bob_share): ");
        for (int i = 0; i < test_dim; i++) {
            total_inputs[i] = alice_shares[i] + input_array[i];
            printf("%llu ", total_inputs[i]);
        }
        printf("\n");
        
        // 计算明文期望MW结果
        printf("Expected MW results: ");
        for (int i = 0; i < test_dim; i++) {
            expected_results[i] = compute_MW_plain(total_inputs[i], test_bw);
            printf("%llu ", expected_results[i]);
        }
        printf("\n");
    }
    
    // 记录通信开销
    size_t comm_start = iopack->io->counter;
    
    // 调用MW协议
    gp->mw(test_dim, input_array, output_array, test_bw, 2);
    for (int i = 0; i < test_dim; i++) {
        printf("input_array[%d] = %llu\n", i, input_array[i]);
        printf("output_array[%d] = %llu\n", i, output_array[i]);
    }
    
    size_t comm_end = iopack->io->counter;
    
    // 恢复明文结果进行验证
    uint64_t *output_alice = new uint64_t[test_dim];
    
    if (party == ALICE) {
        iopack->io->send_data(output_array, test_dim * sizeof(uint64_t));
    } else {
        iopack->io->recv_data(output_alice, test_dim * sizeof(uint64_t));
    }
    
    if (party == BOB) {
        printf("\n=== MW Protocol Results and Verification ===\n");
        int correct_count = 0;
        
        for (int i = 0; i < test_dim; i++) {
            // 恢复明文结果 (2位输出)
            uint64_t mpc_result = (output_alice[i] + output_array[i]) & 3;
            
            bool is_correct = (mpc_result == expected_results[i]);
            if (is_correct) correct_count++;
            
            printf("Test[%d]: Total_Input=%llu, Expected=%llu, MPC=%llu %s\n",
                   i, total_inputs[i], expected_results[i], mpc_result,
                   is_correct ? "✓" : "✗");
            
            if (i < 2) {  // 特别关注用户指定的两个输入
                printf("  -> User specified input: %s\n", 
                       is_correct ? "PASSED" : "FAILED");
            }
        }
        
        printf("\n=== Test Summary ===\n");
        printf("Bit width (bw): %d\n", test_bw);
        printf("Ring size (N): %llu\n", N);
        printf("Total tests: %d\n", test_dim);
        printf("Correct results: %d\n", correct_count);
        printf("Accuracy: %.2f%%\n", (double)correct_count / test_dim * 100.0);
        printf("Communication cost: %zu bytes\n", (comm_end - comm_start));
        
        // 特别报告用户指定输入的结果
        printf("\n=== User Specified Inputs Results ===\n");
        for (int i = 0; i < 2 && i < test_dim; i++) {
            uint64_t mpc_result = (output_alice[i] + output_array[i]) & 3;
            printf("Total Input %llu: Expected=%llu, Got=%llu %s\n",
                   total_inputs[i], expected_results[i], mpc_result,
                   (mpc_result == expected_results[i]) ? "✓" : "✗");
        }
        
        if (correct_count == test_dim) {
            printf("\n🎉 All MW tests PASSED!\n");
        } else {
            printf("\n❌ Some MW tests FAILED.\n");
        }
    }
    
    // 清理内存
    delete[] alice_shares;
    delete[] total_inputs;
    delete[] expected_results;
    delete[] input_array;
    delete[] output_array;
    delete[] output_alice;
    
    printf("=== MW Test Completed ===\n\n");
}

int main(int argc, char **argv) {
    ArgMapping amap;

    amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
    amap.arg("p", port, "Port Number");
    amap.arg("ip", address, "IP Address of server (ALICE)");
    amap.parse(argc, argv);

    // 初始化通信和协议
    iopack = new IOPack(party, port, address);
    otpack = new OTPack(iopack, party);
    gp = new GeometricPerspectiveProtocols(party, iopack, otpack);

    printf("Starting MW Protocol Testing...\n");
    printf("Party: %s\n", (party == ALICE) ? "ALICE" : "BOB");
    
    // 运行MW协议测试
    test_mw_protocol();
    
    printf("================================\n");
    printf("MW Protocol testing completed!\n");

    // 清理资源
    delete gp;
    delete otpack;
    delete iopack;
    
    return 0;
}
