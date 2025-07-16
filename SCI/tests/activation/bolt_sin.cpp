#include "LinearOT/linear-ot.h"
#include "utils/emp-tool.h"
#include <cstdint>
#include <cstdio>
#include <iostream>

#include "FloatingPoint/floating-point.h"
#include "FloatingPoint/fp-math.h"
#include <limits>
#include <math.h>
#include <random>
// #include "float_utils.h"
#include "BuildingBlocks/aux-protocols.h"
#include "BuildingBlocks/geometric_perspective_protocols.h"
#include "BuildingBlocks/truncation.h"
#include "Math/math-functions.h"
#include "Millionaire/equality.h"
#include "Millionaire/millionaire.h"
#include "Millionaire/millionaire_with_equality.h"
#include <chrono>
// #include <matplotlibcpp.h>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

using namespace sci;
using namespace std;
// namespace plt = matplotlibcpp;
// using namespace plt;
int party, port = 32000;
string address = "127.0.0.1";
IOPack *iopack;
OTPack *otpack;
LinearOT *prod;
XTProtocol *ext;

int bwL = 22;
uint64_t N = pow(2, bwL);
uint64_t mask_bwL = (bwL == 64 ? -1 : ((1ULL << bwL) - 1));
bool accumulate = true;
bool precomputed_MSBs = false;
MultMode mode = MultMode::None;

uint64_t f = 10;

uint64_t f_MW = 28;
uint64_t bwL_MW = f_MW+2;
uint64_t N_f_MW = pow(2, f_MW);
uint64_t N_MW = pow(2, bwL_MW);
uint64_t mask_N_MW = (bwL_MW == 64 ? -1 : ((1ULL << bwL_MW) - 1));
uint64_t pow_f = pow(2, f);

int bwL_input = 13;
uint64_t mask_bwL_input = (bwL_input == 64 ? -1 : ((1ULL << bwL_input) - 1));
uint64_t N_input = pow(2, bwL_input);
uint64_t f_input = 10;
uint64_t pow_f_input = pow(2, f_input);

Truncation *trunc_oracle;
AuxProtocols *aux;
MillionaireWithEquality *mill_eq;
Equality *eq;
MathFunctions *math;
GeometricPerspectiveProtocols *gp;
// LinearOT *prod;
// int dim = pow(2, 20);
// int dim = 1;
int dim = 16384;
uint64_t acc = 2;
uint64_t init_input = 0;
uint64_t step_size = 1;
uint64_t correct = 1;

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

void ring_exp(const uint64_t *inA, uint64_t *out, size_t len, uint64_t f,
              uint64_t mask_bwL) {
  double pow_f = std::pow(2.0, f);
  for (size_t i = 0; i < len; ++i) {
    double x0_real = static_cast<double>(inA[i]) / pow_f;
    double exp_val = std::exp(x0_real);
    uint64_t exp_fixed =
        static_cast<uint64_t>(std::round(exp_val * pow_f)) & mask_bwL;
    out[i] = exp_fixed;
  }
}

void ring_sin(const uint64_t *inA, uint64_t *out, size_t len, uint64_t f,
              uint64_t mask_bwL) {
  double pow_f = std::pow(2.0, f);
  for (size_t i = 0; i < len; ++i) {
    double x0_real = static_cast<double>(inA[i]) / pow_f;
    double sin_val = std::sin(x0_real);
    uint64_t sin_fixed =
        static_cast<uint64_t>(std::round(sin_val * pow_f)) & mask_bwL;
    out[i] = sin_fixed;
  }
}

void ring_cos(const uint64_t *inA, uint64_t *out, size_t len, uint64_t f,
              uint64_t mask_bwL) {
  double pow_f = std::pow(2.0, f);
  for (size_t i = 0; i < len; ++i) {
    double x0_real = static_cast<double>(inA[i]) / pow_f;
    double cos_val = std::cos(x0_real);
    uint64_t cos_fixed =
        static_cast<uint64_t>(std::round(cos_val * pow_f)) & mask_bwL;
    out[i] = cos_fixed;
  }
}

double fix2double(uint64_t x, uint64_t y, uint64_t L, uint64_t f) {
  uint64_t mask = (L == 64) ? ~0ULL : ((1ULL << L) - 1);
  uint64_t sum = (x + y) & mask;
  int64_t signed_val;
  if (sum >= (1ULL << (L - 1))) {
    signed_val = (int64_t)(sum - (1ULL << L));
  } else {
    signed_val = (int64_t)sum;
  }
  return (double)signed_val / (1ULL << f);
}

int64_t sign_extend_l(uint64_t x, int l) {
  uint64_t mask = (1ULL << l) - 1;
  x &= mask;
  if (x & (1ULL << (l - 1)))
      return (int64_t)(x | (~mask));
  else
      return (int64_t)x;
}

// x0, x1: 分别为x的高、低部分，定点整数（环上的补码数）
// MW: 区间编号（0,1,2），直接用作查表索引
// N: 区间宽度（整数）
// f: x0/x1的定点小数位数
// l_x: x0/x1的实际位宽
// lut_range: LUT输入的总位宽（如13位）
// lut_frac: LUT的小数位数
// sin_lut, cos_lut: 查表数组
// 返回值：定点整数
int64_t sin_float(
    uint64_t x0, uint64_t x1, int MW, int N, int f, int l_x,
    int lut_range, int lut_frac,
    const std::vector<int64_t>& sin_lut, // 预先填好
    const std::vector<int64_t>& cos_lut  // 预先填好
) {
    // 1. 符号扩展
    int64_t x0_signed = sign_extend_l(x0, l_x);
    int64_t x1_signed = sign_extend_l(x1, l_x);
    // MW 直接用作查表索引
    int idx_MW = MW;

    // 2. x0/x1缩放到LUT精度
    // x0_signed, x1_signed 是以f为小数位的定点数
    // LUT查表输入是以lut_frac为小数位的定点数
    // 需要将x0/x1从f位小数缩放到lut_frac位小数
    int64_t x0_lut = x0_signed;
    int64_t x1_lut = x1_signed;
    if (f > lut_frac) {
        x0_lut = x0_signed >> (f - lut_frac);
        x1_lut = x1_signed >> (f - lut_frac);
    } else if (f < lut_frac) {
        x0_lut = x0_signed << (lut_frac - f);
        x1_lut = x1_signed << (lut_frac - f);
    }
    // 限制查表索引在合法范围
    int idx_x0 = (int)(x0_lut + (1LL << (lut_range-1)));
    int idx_x1 = (int)(x1_lut + (1LL << (lut_range-1)));
    idx_x0 = std::max(0, std::min(idx_x0, (1 << lut_range) - 1));
    idx_x1 = std::max(0, std::min(idx_x1, (1 << lut_range) - 1));

    // 3. 查表
    int64_t sin_x0 = sin_lut[idx_x0];
    int64_t cos_x0 = cos_lut[idx_x0];
    int64_t sin_x1 = sin_lut[idx_x1];
    int64_t cos_x1 = cos_lut[idx_x1];
    int64_t MW_cos = cos_lut[idx_MW]; // cos(-MW * N / 2^f)
    int64_t MW_sin = sin_lut[idx_MW]; // sin(-MW * N / 2^f)

    // 4. 计算每一项（用128位防止溢出）
    __int128 term1 = (__int128)sin_x0 * cos_x1 * MW_cos;
    __int128 term2 = (__int128)cos_x0 * sin_x1 * MW_cos;
    __int128 term3 = (__int128)cos_x0 * cos_x1 * MW_sin;
    __int128 term4 = (__int128)sin_x0 * sin_x1 * MW_sin;

    // 5. 累加
    __int128 sum = term1 + term2 + term3 + term4;

    // 6. 缩放回定点（右移2*lut_frac位，因为三次乘法，每次乘法都带lut_frac位小数）
    int shift = 2 * lut_frac;
    int64_t result = (int64_t)(sum >> shift);

    // 7. 如需环映射，掩码（如只保留低l_x位）
    // result = result & ((1ULL << l_x) - 1);

    return result;
}


int main(int argc, char **argv) {
  ArgMapping amap;

  amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
  amap.arg("p", port, "Port Number");
  amap.arg("ip", address, "IP Address of server (ALICE)");
  amap.arg("m", precomputed_MSBs, "MSB_to_Wrap Optimization?");
  amap.arg("a", ::accumulate, "Accumulate?");
  amap.parse(argc, argv);


  iopack = new IOPack(party, port, "127.0.0.1");
  uint64_t *inA = new uint64_t[dim];    // Declare the variable "inA"
  otpack = new OTPack(iopack, party);
  prod = new LinearOT(party, iopack, otpack);
  gp = new GeometricPerspectiveProtocols(party, iopack, otpack);
  aux = new AuxProtocols(party, iopack, otpack);
  uint64_t *inB = new uint64_t[dim];  // Declare the variable "inB"
  uint64_t *outC = new uint64_t[dim]; // Declare the variable "inC"

  bool signed_arithmetic = true;
  uint64_t *MW = new uint64_t[dim];
  uint64_t *MW_sin_lut = new uint64_t[4];
  uint64_t *MW_cos_lut = new uint64_t[4];
  uint64_t *MW_exp = new uint64_t[dim];
  uint64_t *exp_inA = new uint64_t[dim];
  uint64_t *exp_inB = new uint64_t[dim];
  uint64_t *res_exp = new uint64_t[dim];
  for (int i = 0; i < 3; i++) {
    double pow_f_input = std::pow(2.0, f_input);
    // double pow_f = std::pow(2.0, f_MW);
    double angle = -i * static_cast<double>(N_input) / pow_f_input; // 这里的angle单位是弧度
    double sin_val = std::sin(angle);
    double cos_val = std::cos(angle);
    MW_sin_lut[i] = static_cast<uint64_t>(std::round(sin_val * N_f_MW)) & mask_N_MW;
    MW_cos_lut[i] = static_cast<uint64_t>(std::round(cos_val * N_f_MW)) & mask_N_MW;
    // printf("MW_sin_lut[%d]: %llu\n", i, MW_sin_lut[i]);
    // printf("MW_cos_lut[%d]: %llu\n", i, MW_cos_lut[i]);
  }
  MW_sin_lut[3] = 100;
  MW_cos_lut[3] = 100;


  // for (int j = 0; j < 4; j++) {
  //   MW_sin_lut[j] = MW_sin_lut[j] + (1ULL<<(bwL_MW-1));
  //   MW_cos_lut[j] = MW_cos_lut[j] + (1ULL<<(bwL_MW-1));
  // }

  size_t comm_start = iopack->io->counter;
  for (int i = 0; i < dim; i++) {
    inA[i] = (0 + i * 1 +3 +4832) & mask_bwL_input;
    inB[i] = (init_input+2 +4832+i * 1) & mask_bwL_input;
    compute_MW_plain(inA, inB, MW, dim, N_input);
  }
  uint64_t *sin_inA = new uint64_t[dim];
  uint64_t *cos_inA = new uint64_t[dim];
  uint64_t *sin_inB = new uint64_t[dim];
  uint64_t *cos_inB = new uint64_t[dim];
  if (party == ALICE) {
      ring_sin(inA, sin_inA, dim, f_input, mask_bwL_input);
      ring_cos(inA, cos_inA, dim, f_input, mask_bwL_input);
  } else {
    ring_sin(inA, sin_inA, dim, f_input, mask_bwL_input);
      ring_cos(inA, cos_inA, dim, f_input, mask_bwL_input);
    ring_sin(inB, sin_inB, dim, f_input, mask_bwL_input);
    ring_cos(inB, cos_inB, dim, f_input, mask_bwL_input);
  }

  // if (party == ALICE) {
  // for (int i = 0; i < dim; i++) {
  //   printf("inA[%d]: %llu\n", i, inA[i]);
  //   printf("sin_inA[%d]: %llu\n", i, sin_inA[i]);
  //   printf("cos_inA[%d]: %llu\n", i, cos_inA[i]);
  // }
  // } else {
  //   for (int i = 0; i < dim; i++) {
  //     printf("inB[%d]: %llu\n", i, inB[i]);
  //     printf("sin_inB[%d]: %llu\n", i, sin_inB[i]);
  //     printf("cos_inB[%d]: %llu\n", i, cos_inB[i]);
  //   }
  // }

  uint64_t *sin_inA_cos_inB = new uint64_t[dim];
  uint64_t *cos_inA_sin_inB = new uint64_t[dim];
  uint64_t *cos_inA_cos_inB = new uint64_t[dim];
  uint64_t *sin_inA_sin_inB = new uint64_t[dim];

  uint64_t *zero = new uint64_t[dim];
  for (int i = 0; i < dim; i++) {
    zero[i] = 0;
  }

  if (party == ALICE) {
    prod->hadamard_product(dim, sin_inA, zero, sin_inA_cos_inB, bwL+1, bwL+1, bwL+1+bwL+1, true, true, MultMode::None);
    prod->hadamard_product(dim, cos_inA, zero, cos_inA_sin_inB, bwL+1, bwL+1, bwL+1+bwL+1, true, true, MultMode::None);
    prod->hadamard_product(dim, cos_inA, zero, cos_inA_cos_inB, bwL+1, bwL+1, bwL+1+bwL+1, true, true, MultMode::None);
    prod->hadamard_product(dim, sin_inA, zero, sin_inA_sin_inB, bwL+1, bwL+1, bwL+1+bwL+1, true, true, MultMode::None);
  } else {
    prod->hadamard_product(dim, zero, cos_inB, sin_inA_cos_inB, bwL+1, bwL+1, bwL+1+bwL+1, true, true, MultMode::None);
    prod->hadamard_product(dim, zero, sin_inB, cos_inA_sin_inB, bwL+1, bwL+1, bwL+1+bwL+1, true, true, MultMode::None);
    prod->hadamard_product(dim, zero, cos_inB, cos_inA_cos_inB, bwL+1, bwL+1, bwL+1+bwL+1, true, true, MultMode::None);
    prod->hadamard_product(dim, zero, sin_inB, sin_inA_sin_inB, bwL+1, bwL+1, bwL+1+bwL+1, true, true, MultMode::None);
  }

  /////check the cross_term result
  // for (int i = 0; i < dim; i++) {
  //   printf("sin_inA_cos_inB[%d]: %llu\n", i, sin_inA_cos_inB[i]);
  //   printf("cos_inA_sin_inB[%d]: %llu\n", i, cos_inA_sin_inB[i]);
  //   printf("cos_inA_cos_inB[%d]: %llu\n", i, cos_inA_cos_inB[i]);
  //   printf("sin_inA_sin_inB[%d]: %llu\n", i, sin_inA_sin_inB[i]);
  // }

uint8_t *wrap_sin_inA_cos_inB = new uint8_t[dim];
uint8_t *wrap_cos_inA_sin_inB = new uint8_t[dim];
uint8_t *wrap_cos_inA_cos_inB = new uint8_t[dim];
uint8_t *wrap_sin_inA_sin_inB = new uint8_t[dim];



if(party != ALICE) {
  for (int i = 0; i < dim; i++) {
    sin_inA_cos_inB[i] = sin_inA_cos_inB[i] + (1ULL<<(2*bwL+2  -1));
    cos_inA_sin_inB[i] = cos_inA_sin_inB[i] + (1ULL<<(2*bwL+2  -1));
    cos_inA_cos_inB[i] = cos_inA_cos_inB[i] + (1ULL<<(2*bwL+2  -1));
    sin_inA_sin_inB[i] = sin_inA_sin_inB[i] + (1ULL<<(2*bwL+2  -1));
  }
  }


  aux->wrap_computation(sin_inA_cos_inB, wrap_sin_inA_cos_inB, dim, 2*bwL+2);
  aux->wrap_computation(cos_inA_sin_inB, wrap_cos_inA_sin_inB, dim, 2*bwL+2);
  aux->wrap_computation(cos_inA_cos_inB, wrap_cos_inA_cos_inB, dim, 2*bwL+2);
  aux->wrap_computation(sin_inA_sin_inB, wrap_sin_inA_sin_inB, dim, 2*bwL+2);


//   alice将wrap都发给bob,如果两边的wrap都为1，则赋值为0
uint8_t *wrap_sin_inA_cos_inB_recv = new uint8_t[dim];
uint8_t *wrap_cos_inA_sin_inB_recv = new uint8_t[dim];
uint8_t *wrap_cos_inA_cos_inB_recv = new uint8_t[dim];
uint8_t *wrap_sin_inA_sin_inB_recv = new uint8_t[dim];
if (party == ALICE) {
  iopack->io->send_data(wrap_sin_inA_cos_inB, dim * sizeof(uint8_t));
  iopack->io->send_data(wrap_cos_inA_sin_inB, dim * sizeof(uint8_t));
  iopack->io->send_data(wrap_cos_inA_cos_inB, dim * sizeof(uint8_t));
  iopack->io->send_data(wrap_sin_inA_sin_inB, dim * sizeof(uint8_t));
} else {
    iopack->io->recv_data(wrap_sin_inA_cos_inB_recv, dim * sizeof(uint8_t));
    iopack->io->recv_data(wrap_cos_inA_sin_inB_recv, dim * sizeof(uint8_t));
    iopack->io->recv_data(wrap_cos_inA_cos_inB_recv, dim * sizeof(uint8_t));
    iopack->io->recv_data(wrap_sin_inA_sin_inB_recv, dim * sizeof(uint8_t));
    for (int i = 0; i < dim; i++) {
      if (wrap_sin_inA_cos_inB[i] == 1 && wrap_sin_inA_cos_inB_recv[i] == 1) {
        wrap_sin_inA_cos_inB[i] = 0;
        wrap_sin_inA_cos_inB_recv[i] = 0;
      }
      if (wrap_cos_inA_sin_inB[i] == 1 && wrap_cos_inA_sin_inB_recv[i] == 1) {
        wrap_cos_inA_sin_inB[i] = 0;
        wrap_cos_inA_sin_inB_recv[i] = 0;
      }
      if (wrap_cos_inA_cos_inB[i] == 1 && wrap_cos_inA_cos_inB_recv[i] == 1) {
        wrap_cos_inA_cos_inB[i] = 0;
        wrap_cos_inA_cos_inB_recv[i] = 0;
      }
      if (wrap_sin_inA_sin_inB[i] == 1 && wrap_sin_inA_sin_inB_recv[i] == 1) {
        wrap_sin_inA_sin_inB[i] = 0;
        wrap_sin_inA_sin_inB_recv[i] = 0;
      }
    }
}
if (party != ALICE) {
  iopack->io->send_data(wrap_sin_inA_cos_inB_recv, dim * sizeof(uint8_t));
  iopack->io->send_data(wrap_cos_inA_sin_inB_recv, dim * sizeof(uint8_t));
  iopack->io->send_data(wrap_cos_inA_cos_inB_recv, dim * sizeof(uint8_t));
  iopack->io->send_data(wrap_sin_inA_sin_inB_recv, dim * sizeof(uint8_t));
} else {
    iopack->io->recv_data(wrap_sin_inA_cos_inB, dim * sizeof(uint8_t));
    iopack->io->recv_data(wrap_cos_inA_sin_inB, dim * sizeof(uint8_t));
    iopack->io->recv_data(wrap_cos_inA_cos_inB, dim * sizeof(uint8_t));
    iopack->io->recv_data(wrap_sin_inA_sin_inB, dim * sizeof(uint8_t));
}
// for (int i = 0; i < dim; i++) {
//   printf("wrap_sin_inA_cos_inB[%d]: %d\n", i, wrap_sin_inA_cos_inB[i]);
//   printf("wrap_cos_inA_sin_inB[%d]: %d\n", i, wrap_cos_inA_sin_inB[i]);
//   printf("wrap_cos_inA_cos_inB[%d]: %d\n", i, wrap_cos_inA_cos_inB[i]);
//   printf("wrap_sin_inA_sin_inB[%d]: %d\n", i, wrap_sin_inA_sin_inB[i]);
// }
// uint64_t *prod_output = new uint64_t[dim * 4];
// uint64_t *diff_output = new uint64_t[dim * 4];
uint64_t *sin_cos_cos_buffer = new uint64_t[dim * 4];
uint64_t *cos_sin_cos_buffer = new uint64_t[dim * 4];
uint64_t *cos_cos_sin_buffer = new uint64_t[dim * 4];
uint64_t *sin_sin_sin_buffer = new uint64_t[dim * 4];
uint64_t random_numberbob = 0;
uint64_t random_numberalice = 0;



uint64_t mask_MWlut_outC = (1ULL << (bwL + bwL+2)) - 1;
 if (party != ALICE) {
    for (int i = 0; i < dim; ++i) {
        uint64_t shift = f_MW + f * 2 - f * 2;
      for (uint8_t j = 0; j < 4; ++j) {
        ///////////////////////sin_cos_cos_buffer buffer//////////////////////////////////////////////////////////////////
        __uint128_t prod0 = (__uint128_t)MW_cos_lut[j] * (__uint128_t)sin_inA_cos_inB[i];
        // prod_output[i * 4 + j] = prod;
        // 再减去 MW_lut[j] 和 outC[i]
        __uint128_t diff0 = (__uint128_t)prod0 - (__uint128_t)MW_cos_lut[j]*(1ULL<<(2*bwL+2 -1)) 
        -(__uint128_t)sin_inA_cos_inB[i]*(1ULL<<(bwL_MW -1)) + wrap_sin_inA_cos_inB[i] * (1ULL << (2*bwL+2 + bwL_MW -1)) + (1ULL<<(2*bwL+2 + bwL_MW -2)) - wrap_sin_inA_cos_inB[i] * (__uint128_t)MW_cos_lut[j] ;
        // diff_output[i * 4 + j] = diff;
        uint64_t res0 = (uint64_t)(((diff0 + (__uint128_t)(1ULL << (shift - 1))) >> shift) & mask_MWlut_outC);
        sin_cos_cos_buffer[i * 4 + j] = (res0 + random_numberbob) & mask_MWlut_outC;
        // printf("diff: %llu\n", diff);
        // printf("buffer[%d]: %llu\n", i * 4 + j, buffer[i * 4 + j]);
        ///////////////////////cos_sin_cos_buffer buffer//////////////////////////////////////////////////////////////////
        __uint128_t prod1 = (__uint128_t)MW_cos_lut[j] * (__uint128_t)cos_inA_sin_inB[i];
        __uint128_t diff1 = (__uint128_t)prod1 - (__uint128_t)MW_cos_lut[j]*(1ULL<<(2*bwL+2 -1)) 
        -(__uint128_t)cos_inA_sin_inB[i]*(1ULL<<(bwL_MW -1)) + wrap_cos_inA_sin_inB[i] * (1ULL << (2*bwL+2 + bwL_MW -1)) + (1ULL<<(2*bwL+2 + bwL_MW -2)) - wrap_cos_inA_sin_inB[i] * (__uint128_t)MW_cos_lut[j] ;
        uint64_t res1 = (uint64_t)(((diff1 + (__uint128_t)(1ULL << (shift - 1))) >> shift) & mask_MWlut_outC);
        cos_sin_cos_buffer[i * 4 + j] = (res1 + random_numberbob) & mask_MWlut_outC;
        ///////////////////////cos_cos_sin_buffer buffer//////////////////////////////////////////////////////////////////
        __uint128_t prod2 = (__uint128_t)MW_sin_lut[j] * (__uint128_t)cos_inA_cos_inB[i];
        __uint128_t diff2 = (__uint128_t)prod2 - (__uint128_t)MW_sin_lut[j]*(1ULL<<(2*bwL+2 -1)) 
        -(__uint128_t)cos_inA_cos_inB[i]*(1ULL<<(bwL_MW -1)) + wrap_cos_inA_cos_inB[i] * (1ULL << (2*bwL+2 + bwL_MW -1)) + (1ULL<<(2*bwL+2 + bwL_MW -2)) - wrap_cos_inA_cos_inB[i] * (__uint128_t)MW_sin_lut[j] ;
         // 这里假设fc = f*2，和原代码一致
        uint64_t res2 = (uint64_t)(((diff2 + (__uint128_t)(1ULL << (shift - 1))) >> shift) & mask_MWlut_outC);
        cos_cos_sin_buffer[i * 4 + j] = (res2 + random_numberbob) & mask_MWlut_outC;
        ///////////////////////sin_sin_sin_buffer buffer//////////////////////////////////////////////////////////////////
        __uint128_t prod3 = (__uint128_t)MW_sin_lut[j] * (__uint128_t)sin_inA_sin_inB[i];
        __uint128_t diff3 = (__uint128_t)prod3 - (__uint128_t)MW_sin_lut[j]*(1ULL<<(2*bwL+2 -1)) 
        -(__uint128_t)sin_inA_sin_inB[i]*(1ULL<<(bwL_MW -1)) + wrap_sin_inA_sin_inB[i] * (1ULL << (2*bwL+2 + bwL_MW -1)) + (1ULL<<(2*bwL+2 + bwL_MW -2)) - wrap_sin_inA_sin_inB[i] * (__uint128_t)MW_sin_lut[j] ;
        uint64_t res3 = (uint64_t)(((diff3 + (__uint128_t)(1ULL << (shift - 1))) >> shift) & mask_MWlut_outC);
        sin_sin_sin_buffer[i * 4 + j] = (res3 + random_numberbob) & mask_MWlut_outC;
      }
    }
    iopack->io->send_data(sin_cos_cos_buffer, dim * 4 * sizeof(uint64_t));
    iopack->io->send_data(cos_sin_cos_buffer, dim * 4 * sizeof(uint64_t));
    iopack->io->send_data(cos_cos_sin_buffer, dim * 4 * sizeof(uint64_t));
    iopack->io->send_data(sin_sin_sin_buffer, dim * 4 * sizeof(uint64_t));
    // for (int i = 0; i < dim; i++) {
    //   for (int j = 0; j < 4; j++) {
    //   printf("sin_cos_cos_buffer[%d][%d]: %llu\n", i, j, sin_cos_cos_buffer[i * 4 + j]);
    //   printf("cos_sin_cos_buffer[%d][%d]: %llu\n", i, j, cos_sin_cos_buffer[i * 4 + j]);
    //   printf("cos_cos_sin_buffer[%d][%d]: %llu\n", i, j, cos_cos_sin_buffer[i * 4 + j]);
    //   printf("sin_sin_sin_buffer[%d][%d]: %llu\n", i, j, sin_sin_sin_buffer[i * 4 + j]);
    //   }
    // }
  } else {
    iopack->io->recv_data(sin_cos_cos_buffer, dim * 4 * sizeof(uint64_t));
    iopack->io->recv_data(cos_sin_cos_buffer, dim * 4 * sizeof(uint64_t));
    iopack->io->recv_data(cos_cos_sin_buffer, dim * 4 * sizeof(uint64_t));
    iopack->io->recv_data(sin_sin_sin_buffer, dim * 4 * sizeof(uint64_t));
  }
// todo
  // uint64_t *y = new uint64_t[dim];
  uint64_t *sin_cos_cos_y = new uint64_t[dim ];
  uint64_t *cos_sin_cos_y = new uint64_t[dim ];
  uint64_t *cos_cos_sin_y = new uint64_t[dim ];
  uint64_t *sin_sin_sin_y = new uint64_t[dim ];
  if (party == ALICE) {
    uint64_t **sin_cos_cos_spec = new uint64_t *[dim];
    uint64_t **cos_sin_cos_spec = new uint64_t *[dim];
    uint64_t **cos_cos_sin_spec = new uint64_t *[dim];
    uint64_t **sin_sin_sin_spec = new uint64_t *[dim];
    for (int i = 0; i < dim; ++i) {
      sin_cos_cos_spec[i] = new uint64_t[4];
      cos_sin_cos_spec[i] = new uint64_t[4];
      cos_cos_sin_spec[i] = new uint64_t[4];
      sin_sin_sin_spec[i] = new uint64_t[4];
      uint64_t shift = f_MW + f * 2 - f * 2;
      for (uint8_t j = 0; j < 4; ++j) {
        ///////////////////////sin_cos_cos_spec buffer//////////////////////////////////////////////////////////////////
        __uint128_t prod0 = (__uint128_t)MW_cos_lut[j] * (__uint128_t)sin_cos_cos_buffer[i];
        __uint128_t diff0 = (__uint128_t)prod0 - (__uint128_t)sin_inA_cos_inB[i]*(1ULL<<(bwL_MW -1)) + wrap_sin_inA_cos_inB[i] * (1ULL << (2*bwL+2 + bwL_MW -1)) - wrap_sin_inA_cos_inB[i] * (__uint128_t)MW_cos_lut[j] ;
         uint64_t res0 = (uint64_t)(((diff0 + (__uint128_t)(1ULL << (shift - 1))) >> shift) & mask_MWlut_outC);
        sin_cos_cos_spec[i][j] = (res0 + random_numberalice +sin_cos_cos_buffer[i * 4 + j]) & mask_MWlut_outC;
        ///////////////////////cos_sin_cos_spec buffer//////////////////////////////////////////////////////////////////
        __uint128_t prod1 = (__uint128_t)MW_cos_lut[j] * (__uint128_t)cos_sin_cos_buffer[i];
        __uint128_t diff1 = (__uint128_t)prod1 - (__uint128_t)cos_inA_sin_inB[i]*(1ULL<<(bwL_MW -1)) + wrap_cos_inA_sin_inB[i] * (1ULL << (2*bwL+2 + bwL_MW -1)) - wrap_cos_inA_sin_inB[i] * (__uint128_t)MW_cos_lut[j] ;
        uint64_t res1 = (uint64_t)(((diff1 + (__uint128_t)(1ULL << (shift - 1))) >> shift) & mask_MWlut_outC);
        cos_sin_cos_spec[i][j] = (res1 + random_numberalice +cos_sin_cos_buffer[i * 4 + j]) & mask_MWlut_outC;
        ///////////////////////cos_cos_sin_spec buffer//////////////////////////////////////////////////////////////////
        __uint128_t prod2 = (__uint128_t)MW_sin_lut[j] * (__uint128_t)cos_cos_sin_buffer[i];
        __uint128_t diff2 = (__uint128_t)prod2 - (__uint128_t)cos_inA_cos_inB[i]*(1ULL<<(bwL_MW -1)) + wrap_cos_inA_cos_inB[i] * (1ULL << (2*bwL+2 + bwL_MW -1)) - wrap_cos_inA_cos_inB[i] * (__uint128_t)MW_sin_lut[j] ;
        uint64_t res2 = (uint64_t)(((diff2 + (__uint128_t)(1ULL << (shift - 1))) >> shift) & mask_MWlut_outC);
        cos_cos_sin_spec[i][j] = (res2 + random_numberalice +cos_cos_sin_buffer[i * 4 + j]) & mask_MWlut_outC;
        ///////////////////////sin_sin_sin_spec buffer//////////////////////////////////////////////////////////////////
        __uint128_t prod3 = (__uint128_t)MW_sin_lut[j] * (__uint128_t)sin_sin_sin_buffer[i];
        __uint128_t diff3 = (__uint128_t)prod3 - (__uint128_t)sin_inA_sin_inB[i]*(1ULL<<(bwL_MW -1)) + wrap_sin_inA_sin_inB[i] * (1ULL << (2*bwL+2 + bwL_MW -1)) - wrap_sin_inA_sin_inB[i] * (__uint128_t)MW_sin_lut[j] ;
        uint64_t res3 = (uint64_t)(((diff3 + (__uint128_t)(1ULL << (shift - 1))) >> shift) & mask_MWlut_outC);
        sin_sin_sin_spec[i][j] = (res3 + random_numberalice +sin_sin_sin_buffer[i * 4 + j]) & mask_MWlut_outC;
      }
    }
    aux->lookup_table<uint64_t>(sin_cos_cos_spec, nullptr, nullptr, dim, 2,
                                 bwL + bwL+2);
    aux->lookup_table<uint64_t>(cos_sin_cos_spec, nullptr, nullptr, dim, 2,
                                 bwL + bwL+2);
    aux->lookup_table<uint64_t>(cos_cos_sin_spec, nullptr, nullptr, dim, 2,
                                 bwL + bwL+2);
    aux->lookup_table<uint64_t>(sin_sin_sin_spec, nullptr, nullptr, dim, 2,
                                 bwL + bwL+2);

    // for (int i = 0; i < dim; i++) {
    //   for (int j = 0; j < 4; j++) {
    //   printf("sin_cos_cos_spec[%d][%d]: %llu\n", i, j, sin_cos_cos_spec[i][j]);
    //   printf("cos_sin_cos_spec[%d][%d]: %llu\n", i, j, cos_sin_cos_spec[i][j]);
    //   printf("cos_cos_sin_spec[%d][%d]: %llu\n", i, j, cos_cos_sin_spec[i][j]);
    //   printf("sin_sin_sin_spec[%d][%d]: %llu\n", i, j, sin_sin_sin_spec[i][j]);
    //   }
    // }
    for (int i = 0; i < dim; ++i) {
      delete[] sin_cos_cos_spec[i];
      delete[] cos_sin_cos_spec[i];
      delete[] cos_cos_sin_spec[i];
      delete[] sin_sin_sin_spec[i];
    }
    delete[] sin_cos_cos_spec;
    delete[] cos_sin_cos_spec;
    delete[] cos_cos_sin_spec;
    delete[] sin_sin_sin_spec;
  } else  {
    // printf("party is %d\n",party);
    aux->lookup_table<uint64_t>(nullptr, MW, sin_cos_cos_y, dim, 2,bwL + bwL+2);
    aux->lookup_table<uint64_t>(nullptr, MW, cos_sin_cos_y, dim, 2,bwL + bwL+2);
    aux->lookup_table<uint64_t>(nullptr, MW, cos_cos_sin_y, dim, 2,bwL + bwL+2);
    aux->lookup_table<uint64_t>(nullptr, MW, sin_sin_sin_y, dim, 2,bwL + bwL+2);
  }
  if (party != ALICE) {
    for (int i = 0; i < dim; i++) {
      res_exp[i] =
          ((sin_cos_cos_y[i] - random_numberbob) +
          (cos_sin_cos_y[i] - random_numberbob) +
          (cos_cos_sin_y[i] - random_numberbob) +
          (sin_sin_sin_y[i] - random_numberbob) )& mask_MWlut_outC;
    }
  } else {
    for (int i = 0; i < dim; i++) {
      res_exp[i] =
          (- random_numberalice *4) & mask_MWlut_outC; //这里的位宽到底选多少？
    }
  }
  

  ///计算明文值
  double *res_exp_plain = new double[dim];
  double pow_double = pow(2, 10.0);      // 1024.0
  double N_input_double = pow(2, 13.0);  // 8192.0
  if (party != ALICE) {
    for (int i = 0; i < dim; i++) {
      double inA_fp = static_cast<double>(inA[i]) / pow_double;
      double inB_fp = static_cast<double>(inB[i]) / pow_double;
      double MW_angle = -static_cast<double>(MW[i]) * N_input_double / pow_double;

      double res = std::sin(inA_fp) * std::cos(inB_fp) * std::cos(MW_angle)
           + std::cos(inA_fp) * std::sin(inB_fp) * std::cos(MW_angle)
           + std::cos(inA_fp) * std::cos(inB_fp) * std::sin(MW_angle)
           - std::sin(inA_fp) * std::sin(inB_fp) * std::sin(MW_angle);
      res_exp_plain[i] = res;
      // printf("inA[%d]: %d\n", i, inA[i]);
      // printf("inB[%d]: %d\n", i, inB[i]);
      // printf("MW_fp: %f\n", i, MW_angle);
      // printf("N_input: %d\n", N_input_double);
      // printf("pow_f: %d\n", pow_double);
      // printf("std::sin(inA[%d] / pow_f): %f\n", i, std::sin(inA[i] / pow_f));
      // printf("std::cos(inB[%d] / pow_f): %f\n", i, std::cos(inB[i] / pow_f));
      // // printf("N_input_fp: %f\n",  N_input_fp);
      // printf("-MW[%d] * N_input / pow_f: %f\n", i, -MW[i] * N_input / pow_f);
      // printf("DEBUG: -MW[%d] * N_input / pow_f = %f\n", i, -MW[i] * N_input / (double)pow_f);
      // printf("DEBUG: -MW[%d] * N_input_double / pow_double = %f\n", i, -MW[i] * N_input_double / pow_double);
      // printf("std::cos(-MW[%d] * N_input_double / pow_double): %f\n", i, std::cos(MW_angle));
      // printf("std::sin(inA[%d] / pow_f): %f\n", i, std::sin(inA[i] / pow_f));
      // printf("std::sin(inB[%d] / pow_f): %f\n", i, std::sin(inB[i] / pow_f));
      // printf("std::sin(-MW[%d] * N_input / pow_f): %f\n", i, std::sin(-MW[i] * N_input / pow_f));
      // printf("res[%d]: %f\n", i, res);
      // printf("res_exp_plain[%d]: %f\n", i, res_exp_plain[i]);
    }
  }
  //   for (int i = 0; i < dim; i++) {
  //     // int64_t sin_inA_int = sign_extend_l(sin_inA[i], bwL_input);
  //     // int64_t cos_inB_int = sign_extend_l(cos_inB[i], bwL_input);
  //     // int64_t cos_inA_int = sign_extend_l(cos_inA[i], bwL_input);
  //     // int64_t sin_inB_int = sign_extend_l(sin_inB[i], bwL_input);
  //     // uint64_t MW_cos_lut_int = sign_extend_l(MW_cos_lut[MW[i]], bwL_MW);
  //     // uint64_t MW_sin_lut_int = sign_extend_l(MW_sin_lut[MW[i]], bwL_MW);
  //     uint64_t term1 = sin_inA[i] * cos_inB[i] * MW_cos_lut[MW[i]];
  //     uint64_t term2 = cos_inA[i] * sin_inB[i] * MW_cos_lut[MW[i]];
  //     uint64_t term3 = cos_inA[i] * cos_inB[i] * MW_sin_lut[MW[i]];
  //     uint64_t term4 = sin_inA[i] * sin_inB[i] * MW_sin_lut[MW[i]];
    
  //   uint64_t sum = (term1 + term2 + term3 + term4) & ((1ULL << 13*2+30+1) - 1);
  //   sum >>= f_MW; // 先右移
  //   uint64_t sum_ring = (uint64_t)sum & ((1ULL << 26+1  ) - 1);
  //   res_exp[i] = sum_ring;
  //   printf("res_exp[%d]: %llu\n", i, res_exp[i]);
  //   printf("sin_inA[%d]: %llu\n", i, sin_inA[i]);
  //   printf("cos_inB[%d]: %llu\n", i, cos_inB[i]);
  //   printf("cos_inA[%d]: %llu\n", i, cos_inA[i]);
  //   printf("sin_inB[%d]: %llu\n", i, sin_inB[i]);
  //   printf("MW_cos_lut[%d]: %llu\n", i, MW_cos_lut[MW[i]]);
  //   printf("MW_sin_lut[%d]: %llu\n", i, MW_sin_lut[MW[i]]);
  //   // printf("sin_inA_int: %lld\n", sin_inA_int);
  //   // printf("cos_inB_int: %lld\n", cos_inB_int);
  //   // printf("cos_inA_int: %lld\n", cos_inA_int);
  //   // printf("sin_inB_int: %lld\n", sin_inB_int);
  //   // printf("MW_cos_lut_int: %lld\n", MW_cos_lut_int);
  //   // printf("MW_sin_lut_int: %lld\n", MW_sin_lut_int);
  //   printf("term1: %llu\n", term1);
  //   printf("term2: %llu\n", term2);
  //   printf("term3: %llu\n", term3);
  //   printf("term4: %llu\n", term4);
  //   printf("sum: %llu\n", sum);
  //   printf("res_exp[%d]: %llu\n", i, res_exp[i]);
  //   }
  // }

  
  uint64_t *res_exp_alice = new uint64_t[dim];
    if (party == ALICE) {
        iopack->io->send_data(res_exp, dim * sizeof(uint64_t));
    } else {
        iopack->io->recv_data(res_exp_alice, dim * sizeof(uint64_t));
    }
    double *res_sin_plain = new double[dim];
    double *ideal_exp_plain = new double[dim];
    double ulp = 1.0 / (1 << 10); 
    mask_MWlut_outC = (1ULL << (26+1)) - 1;
    for (int i = 0; i < dim; i++) {
        res_sin_plain[i] = static_cast<double>((res_exp_alice[i] + res_exp[i]) & mask_MWlut_outC) / std::pow(2, 2*f);
        ideal_exp_plain[i] = std::sin(fix2double(inA[i], inB[i], bwL_input, f_input));
        // printf("MW: %d\n", MW[i]);
        // printf("inA[%d]: %llu\n", i, inA[i]);
        // printf("inB[%d]: %llu\n", i, inB[i]);
        // printf("res_exp[%d]: %llu\n", i, res_exp[i]);
        printf("res_sin_plain[%d]: %f\n", i, res_exp_plain[i]);
        printf("fix2double(inA[%d], inB[%d], bwL_input, f_input): %f\n", i, i, fix2double(inA[i], inB[i], bwL_input, f_input));
        printf("ideal_exp_plain[%d]: %f\n", i, ideal_exp_plain[i]);
        double ulp_error = fabs(ideal_exp_plain[i] - res_exp_plain[i]) / ulp;
        printf("ULP [%d]: %.2f\n", i, ulp_error);
    }
        double ulp_sum = 0.0;
    double ulp_max = 0.0;
    for (int i = 0; i < dim; i++) {
    double ulp = 1.0 / (1 << 10); // 2的12次精度
    double error = std::abs(res_exp_plain[i] - ideal_exp_plain[i]);
    ulp_sum += error / ulp;
    ulp_max = std::max(ulp_max, error / ulp);
  }
  printf("ULP sum: %f\n", ulp_sum);
  printf("ULP max: %f\n", ulp_max);

  delete prod;
  delete[] inA; // Delete the variable "inA" to avoid memory leaks
  delete[] inB;
}