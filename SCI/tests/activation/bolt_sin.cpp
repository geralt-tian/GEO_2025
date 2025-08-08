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

uint64_t f = 12;

uint64_t f_MW = 28;
uint64_t bwL_MW = f_MW+2;
uint64_t N_f_MW = pow(2, f_MW);
uint64_t N_MW = pow(2, bwL_MW);
uint64_t mask_N_MW = (bwL_MW == 64 ? -1 : ((1ULL << bwL_MW) - 1));
uint64_t pow_f = pow(2, f);

int bwL_input = 15;
uint64_t mask_bwL_input = (bwL_input == 64 ? -1 : ((1ULL << bwL_input) - 1));
uint64_t N_input = pow(2, bwL_input);
uint64_t f_input = 12;
uint64_t pow_f_input = pow(2, f_input);

uint64_t f_output = 12;
uint64_t bwL_output = f_output+10;

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
    __int128 sum = term1 + term2 + term3 - term4;

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
  ext = new XTProtocol(party, iopack, otpack);
  uint64_t *inB = new uint64_t[dim];  // Declare the variable "inB"
  uint64_t *outC = new uint64_t[dim]; // Declare the variable "inC"

  bool signed_arithmetic = true;
  


  // for (int j = 0; j < 4; j++) {
  //   MW_sin_lut[j] = MW_sin_lut[j] + (1ULL<<(bwL_MW-1));
  //   MW_cos_lut[j] = MW_cos_lut[j] + (1ULL<<(bwL_MW-1));
  // }
  uint64_t *MW = new uint64_t[dim];
  size_t comm_start = iopack->io->counter;
  for (int i = 0; i < dim; i++) {
    inA[i] = (0 + i * 20 +3 ) & mask_bwL_input;
    inB[i] = (init_input+2 +i * 11) & mask_bwL_input;
    compute_MW_plain(inA, inB, MW, dim, N_input);
  }

  //step 1: set bitwidth
  uint64_t f_t = 12;
  uint64_t bwL_t = f_t + 2;
  uint64_t mask_bwL_t = (bwL_t == 64 ? -1 : ((1ULL << bwL_t) - 1));
  uint64_t f_T = 30;
  uint64_t bwL_T = f_T + 2;
  uint64_t mask_bwL_T = (bwL_T == 64 ? -1 : ((1ULL << bwL_T) - 1));
  uint64_t N_f_T = pow(2, f_T);

  
  uint64_t *MW_sin_lut = new uint64_t[4];
  uint64_t *MW_cos_lut = new uint64_t[4];
  uint64_t *MW_exp = new uint64_t[dim];
  uint64_t *exp_inA = new uint64_t[dim];
  uint64_t *exp_inB = new uint64_t[dim];
  uint64_t *res_exp = new uint64_t[dim];
  for (int i = 0; i < 3; i++) {
    double pow_f_input = std::pow(2.0, f_input);
    double pow_f = std::pow(2.0, f_T);
    double angle = -i * static_cast<double>(N_input) / pow_f_input; // 这里的angle单位是弧度
    double sin_val = std::sin(angle);
    double cos_val = std::cos(angle);
    printf("sin_val: %f\n", sin_val);
    printf("cos_val: %f\n", cos_val);
    MW_sin_lut[i] = static_cast<uint64_t>(std::round(sin_val * N_f_T)) & mask_bwL_T;
    MW_cos_lut[i] = static_cast<uint64_t>(std::round(cos_val * N_f_T)) & mask_bwL_T;
    printf("MW_sin_lut[%d]: %llu\n", i, MW_sin_lut[i]);
    printf("MW_cos_lut[%d]: %llu\n", i, MW_cos_lut[i]);
  }
  MW_sin_lut[3] = 100;
  MW_cos_lut[3] = 100;
  //step 2: compute sin and cos
  uint64_t *sin_inA = new uint64_t[dim];
  uint64_t *cos_inA = new uint64_t[dim];
  uint64_t *sin_inB = new uint64_t[dim];
  uint64_t *cos_inB = new uint64_t[dim];
  if (party == ALICE) {
      ring_sin(inA, sin_inA, dim, f_t, mask_bwL_t);
      ring_cos(inA, cos_inA, dim, f_t, mask_bwL_t);
  } else {
    ring_sin(inB, sin_inB, dim, f_t, mask_bwL_t);
    ring_cos(inB, cos_inB, dim, f_t, mask_bwL_t);
  }

  //step3: postive  sin cos

  // if (party == ALICE) {
  //   for (int i = 0; i < dim; i++) {
  //     sin_inA[i] = sin_inA[i] + (1ULL<<(f_t - 1));
  //     cos_inA[i] = cos_inA[i] + (1ULL<<(f_t - 1));
  //   }
  // } else {
  //   for (int i = 0; i < dim; i++) {
  //     sin_inB[i] = sin_inB[i] + (1ULL<<(f_t - 1));
  //     cos_inB[i] = cos_inB[i] + (1ULL<<(f_t - 1));
  //   }
  // }





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


  //step 4: compute sin_inA_cos_inB, cos_inA_sin_inB, cos_inA_cos_inB, sin_inA_sin_inB 目前计算的是signed mul，可优化成signed crossterm
  uint64_t *sin_inA_cos_inB = new uint64_t[dim];
  uint64_t *cos_inA_sin_inB = new uint64_t[dim];
  uint64_t *cos_inA_cos_inB = new uint64_t[dim];
  uint64_t *sin_inA_sin_inB = new uint64_t[dim];

  uint64_t *zero = new uint64_t[dim];
  for (int i = 0; i < dim; i++) {
    zero[i] = 0;
  }

  if (party == ALICE) {
    prod->hadamard_product(dim, sin_inA, zero, sin_inA_cos_inB, bwL_t, bwL_t, bwL_t+bwL_t, true, true, MultMode::None);
    prod->hadamard_product(dim, cos_inA, zero, cos_inA_sin_inB, bwL_t, bwL_t, bwL_t+bwL_t, true, true, MultMode::None);
    prod->hadamard_product(dim, cos_inA, zero, cos_inA_cos_inB, bwL_t, bwL_t, bwL_t+bwL_t, true, true, MultMode::None);
    prod->hadamard_product(dim, sin_inA, zero, sin_inA_sin_inB, bwL_t, bwL_t, bwL_t+bwL_t, true, true, MultMode::None);
  } else {
    prod->hadamard_product(dim, zero, cos_inB, sin_inA_cos_inB, bwL_t, bwL_t, bwL_t+bwL_t, true, true, MultMode::None);
    prod->hadamard_product(dim, zero, sin_inB, cos_inA_sin_inB, bwL_t, bwL_t, bwL_t+bwL_t, true, true, MultMode::None);
    prod->hadamard_product(dim, zero, cos_inB, cos_inA_cos_inB, bwL_t, bwL_t, bwL_t+bwL_t, true, true, MultMode::None);
    prod->hadamard_product(dim, zero, sin_inB, sin_inA_sin_inB, bwL_t, bwL_t, bwL_t+bwL_t, true, true, MultMode::None);
  }

  
  for (int i = 0; i < 100; i++) {
    printf("\n");
    if (party == ALICE) {
      printf("sin_inA[%d]: %llu\n", i, sin_inA[i]);
      printf("cos_inA[%d]: %llu\n", i, cos_inA[i]);
      printf("sin_inA_cos_inB[%d]: %llu\n", i, sin_inA_cos_inB[i]);
      printf("cos_inA_sin_inB[%d]: %llu\n", i, cos_inA_sin_inB[i]);
      printf("cos_inA_cos_inB[%d]: %llu\n", i, cos_inA_cos_inB[i]);
      printf("sin_inA_sin_inB[%d]: %llu\n", i, sin_inA_sin_inB[i]);
    } else {
      printf("sin_inB[%d]: %llu\n", i, sin_inB[i]);
      printf("cos_inB[%d]: %llu\n", i, cos_inB[i]);
      printf("sin_inA_cos_inB[%d]: %llu\n", i, sin_inA_cos_inB[i]);
      printf("cos_inA_sin_inB[%d]: %llu\n", i, cos_inA_sin_inB[i]);
      printf("cos_inA_cos_inB[%d]: %llu\n", i, cos_inA_cos_inB[i]);
      printf("sin_inA_sin_inB[%d]: %llu\n", i, sin_inA_sin_inB[i]);
    }
  }

  uint64_t *sc_add_cs_lut_buffer = new uint64_t[dim * 4];
  uint64_t *cc_min_ss_lut_buffer = new uint64_t[dim * 4];
  uint64_t mask_2bwL_t = (2*bwL_t == 64 ? -1 : ((1ULL << (2*bwL_t)) - 1));
  //step 14: compute select table
  for (int i = 0; i < 4; i++) {
    uint64_t *sc_add_cs = new uint64_t[dim];
    uint64_t *cc_min_ss = new uint64_t[dim];
    uint64_t *sc_add_cs_lut = new uint64_t[dim];
    uint64_t *cc_min_ss_lut = new uint64_t[dim];
    for (int j = 0; j < dim; j++) {
      sc_add_cs[j] = (sin_inA_cos_inB[j] + cos_inA_sin_inB[j]) & mask_2bwL_t;
      cc_min_ss[j] = (cos_inA_cos_inB[j] - sin_inA_sin_inB[j]) & mask_2bwL_t;
    }
    uint64_t *MW_cos_lut_extend = new uint64_t[dim];
    uint64_t *MW_sin_lut_extend = new uint64_t[dim];
    for (int j = 0; j < dim; j++) {
      MW_cos_lut_extend[j] = MW_cos_lut[i];
      MW_sin_lut_extend[j] = MW_sin_lut[i];
    }
    if (party == ALICE) {
      prod->hadamard_product(dim, sc_add_cs , MW_cos_lut_extend, sc_add_cs_lut, 2*bwL_t, bwL_T, 2*bwL_t+bwL_T, true, true, MultMode::None);
      prod->hadamard_product(dim, cc_min_ss , MW_sin_lut_extend, cc_min_ss_lut, 2*bwL_t, bwL_T, 2*bwL_t+bwL_T, true, true, MultMode::None);
    } else {
      prod->hadamard_product(dim, sc_add_cs , zero         , sc_add_cs_lut, 2*bwL_t, bwL_T, 2*bwL_t+bwL_T, true, true, MultMode::None);
      prod->hadamard_product(dim, cc_min_ss , zero         , cc_min_ss_lut, 2*bwL_t, bwL_T, 2*bwL_t+bwL_T, true, true, MultMode::None);
    }
    for (int j = 0; j < 100; j++) {
      printf("sc_add_cs[%d][%d]: %llu\n", i, j, sc_add_cs[j]);
      printf("cc_min_ss[%d][%d]: %llu\n", i, j, cc_min_ss[j]);
      printf("MW_cos_lut_extend[%d][%d]: %llu\n", i, j, MW_cos_lut_extend[j]);
      printf("MW_sin_lut_extend[%d][%d]: %llu\n", i, j, MW_sin_lut_extend[j]);
      printf("sc_add_cs_lut[%d][%d]: %llu\n", i, j, sc_add_cs_lut[j]);
      printf("cc_min_ss_lut[%d][%d]: %llu\n", i, j, cc_min_ss_lut[j]);
      
    }
    // 写入：按 i 主序（每个 i 是一个连续的 dim 块）
    for (int j = 0; j < dim; j++) {
      sc_add_cs_lut_buffer[i * dim + j] = sc_add_cs_lut[j];
      cc_min_ss_lut_buffer[i * dim + j] = cc_min_ss_lut[j];
    }
  }
  //step 16: send lut buffer to bob
  uint64_t *sc_add_cs_lut_buffer_recv = new uint64_t[dim * 4];
  uint64_t *cc_min_ss_lut_buffer_recv = new uint64_t[dim * 4];
  uint64_t mask_2bwL_t_bwL_T = (2*bwL_t+bwL_T == 64 ? -1 : ((1ULL << (2*bwL_t+bwL_T)) - 1));
  if (party != ALICE) {
    iopack->io->send_data(sc_add_cs_lut_buffer, dim * 4 * sizeof(uint64_t));
    iopack->io->send_data(cc_min_ss_lut_buffer, dim * 4 * sizeof(uint64_t));
  } else {
    iopack->io->recv_data(sc_add_cs_lut_buffer_recv, dim * 4 * sizeof(uint64_t));
    iopack->io->recv_data(cc_min_ss_lut_buffer_recv, dim * 4 * sizeof(uint64_t));
    for (int i = 0; i < dim; i++) {
      for (int j = 0; j < 4; j++) {
        sc_add_cs_lut_buffer[i * 4 + j] = (sc_add_cs_lut_buffer[i * 4 + j] + sc_add_cs_lut_buffer_recv[i * 4 + j]) & mask_2bwL_t_bwL_T;
        cc_min_ss_lut_buffer[i * 4 + j] = (cc_min_ss_lut_buffer[i * 4 + j] + cc_min_ss_lut_buffer_recv[i * 4 + j]) & mask_2bwL_t_bwL_T;
      }
    }
  }


  // uint64_t *temp0 = new uint64_t[dim];
  // uint64_t *temp1 = new uint64_t[dim];
  uint64_t *T_add_T = new uint64_t[dim];
  if (party == ALICE) {
    uint64_t **sc_cs_lut_spec = new uint64_t *[dim];
    // uint64_t **cc_ss_lut_spec = new uint64_t *[dim];
    for (int i = 0; i < dim; i++) {
      sc_cs_lut_spec[i] = new uint64_t[4];
      // cc_ss_lut_spec[i] = new uint64_t[4];
      for (int j = 0; j < 4; j++) {
        sc_cs_lut_spec[i][j] = (sc_add_cs_lut_buffer[j * dim + i] + cc_min_ss_lut_buffer[j * dim + i]) & mask_2bwL_t_bwL_T;
        // cc_ss_lut_spec[i][j] = ;
      }
    }
    aux->lookup_table<uint64_t>(sc_cs_lut_spec, nullptr, nullptr, dim, 2, 2*bwL_t+bwL_T);
    // aux->lookup_table<uint64_t>(cc_ss_lut_spec, nullptr, nullptr, dim, 2, 2*bwL_t+bwL_T);
  }
  else {
    aux->lookup_table<uint64_t>(nullptr, MW, T_add_T, dim, 2,2*bwL_t+bwL_T);

  }


  uint64_t *T_add_T_reduce = new uint64_t[dim];
  for (int i = 0; i < dim; i++) {
    T_add_T_reduce[i] = T_add_T[i] >> (2*f_t+f_T-f_output);
  }

  ext->s_extend(dim, T_add_T_reduce, res_exp, f_output+5, bwL_output, nullptr);



  
  uint64_t *res_exp_alice = new uint64_t[dim];
    if (party == ALICE) {
        iopack->io->send_data(res_exp, dim * sizeof(uint64_t));
    } else {
        iopack->io->recv_data(res_exp_alice, dim * sizeof(uint64_t));
    }
    double *res_sin_plain = new double[dim];
    double *ideal_exp_plain = new double[dim];
    double ulp = 1.0 / (1 << f_input); 
    uint64_t mask_MWlut_outC = (1ULL << bwL_output) - 1;
    for (int i = 0; i < 1000; i++) {
        // res_sin_plain[i] = static_cast<double>((res_exp_alice[i] + res_exp[i]) & mask_MWlut_outC) / std::pow(2, f_output);
        res_sin_plain[i] = fix2double(res_exp_alice[i], res_exp[i], bwL_output, f_output);
        ideal_exp_plain[i] = std::sin(fix2double(inA[i], inB[i], bwL_input, f_input));
        // printf("MW: %d\n", MW[i]);
        // printf("inA[%d]: %llu\n", i, inA[i]);
        // printf("inB[%d]: %llu\n", i, inB[i]);
        // printf("res_exp[%d]: %llu\n", i, res_exp[i]);
        printf("inA[%u]: %u\n", i, inA[i]);
        printf("inB[%u]: %u\n", i, inB[i]);
        printf("MW: %d\n", MW[i]);
        printf("res_exp_alice[%d]: %llu\n", i, res_exp_alice[i]);
        printf("res_exp[%d]: %llu\n", i, res_exp[i]);
        printf("res_sin_plain[%d]: %.10f\n", i, res_sin_plain[i]);
        printf("fix2double(inA[%d], inB[%d], bwL_input, f_input): %f\n", i, i, fix2double(inA[i], inB[i], bwL_input, f_input));
        printf("ideal_exp_plain[%d]: %.10f\n", i, ideal_exp_plain[i]);
        double ulp_error = fabs(ideal_exp_plain[i] - res_sin_plain[i]) / ulp;
        printf("ULP [%d]: %.6f\n", i, ulp_error);
    }
    for (int i = 0; i < 4; i++) {
        printf("MW_cos_lut[%d]: %llu\n", i, MW_cos_lut[i]);
        printf("MW_sin_lut[%d]: %llu\n", i, MW_sin_lut[i]);
    }
        double ulp_sum = 0.0;
    double ulp_max = 0.0;
    for (int i = 0; i < dim; i++) {
    double ulp = 1.0 / (1 << f_input); // 2的f_input次精度 
    double error = std::abs(res_sin_plain[i] - ideal_exp_plain[i]);
    ulp_sum += error / ulp;
    ulp_max = std::max(ulp_max, error / ulp);
  }
  printf("ULP avg: %f\n", ulp_sum / dim);
  printf("ULP max: %f\n", ulp_max);

  delete prod;
  delete[] inA; // Delete the variable "inA" to avoid memory leaks
  delete[] inB;
}