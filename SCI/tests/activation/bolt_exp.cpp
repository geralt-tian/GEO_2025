#include "LinearOT/linear-ot.h"
#include "utils/emp-tool.h"
#include <cstdint>
#include <cstdio>
#include <iostream>

#include "FloatingPoint/floating-point.h"
#include "FloatingPoint/fp-math.h"
#include <limits>
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
  uint64_t *inzero = new uint64_t[dim]; // Declare the variable "inzero"
  otpack = new OTPack(iopack, party);
  prod = new LinearOT(party, iopack, otpack);
  gp = new GeometricPerspectiveProtocols(party, iopack, otpack);
  aux = new AuxProtocols(party, iopack, otpack);
  uint64_t *inB = new uint64_t[dim];  // Declare the variable "inB"
  uint64_t *outC = new uint64_t[dim]; // Declare the variable "inC"

  bool signed_arithmetic = true;
  uint64_t *MW = new uint64_t[dim];
  uint64_t *MW_lut = new uint64_t[4];
  uint64_t *MW_exp = new uint64_t[dim];
  uint64_t *exp_inA = new uint64_t[dim];
  uint64_t *exp_inB = new uint64_t[dim];
  uint64_t *res_exp = new uint64_t[dim];
  for (int i = 0; i < 3; i++) {
    double pow_f_input = std::pow(2.0, f_input);
    double pow_f = std::pow(2.0, f_MW);
    double exp_val = std::exp(-i * static_cast<double>(N_input) / pow_f_input);
    MW_lut[i] = static_cast<uint64_t>(std::round(exp_val * N_f_MW)) & mask_N_MW;
    // printf("exp_val: %.10f\n", exp_val);
    // printf("MW_lut[%d]: %llu\n", i, MW_lut[i]);
  }
  MW_lut[3] = 100;

  size_t comm_start = iopack->io->counter;
  for (int i = 0; i < dim; i++) {
    // inA[i] = (6873) & mask_bwL_input;
    // inB[i] = (6879) & mask_bwL_input;
    inA[i] = (0 + i * 1 -3) & mask_bwL_input;
    inB[i] = (init_input+3 +i * 1) & mask_bwL_input;
    // inA[i] = (0 + i * 0 + 3) & mask_bwL_input;
    // inB[i] = (init_input + 2) & mask_bwL_input;
    // inA[i] = (0 + i * 0 + 4096+3) & mask_bwL_input;
    // inB[i] = (init_input + 4096+2) & mask_bwL_input;
    // inA[i] = (0 + i * 0 + 8192-3) & mask_bwL_input;
    // inB[i] = (init_input + 8192-2) & mask_bwL_input;
    // inB[i] = (init_input + i * step_size) & mask_bwL;
    compute_MW_plain(inA, inB, MW, dim, N_input);
    // compute_neg_MW_plain(inA, inB, MW_lut, dim, N, f);
    // printf("MW[%d]: %llu\n", i, MW[i]);

    inzero[i] = 0;
  }

  if (party == ALICE) {
    for (int i = 0; i < dim; i++) {
      MW[i] = 0;
      // MW_lut[i] = 0;
      exp_inB[i] = 0;
    }
    ring_exp(inA, exp_inA, dim, f, mask_bwL);
  } else {
    ring_exp(inB, exp_inB, dim, f, mask_bwL);
    for (int i = 0; i < dim; i++) {
      exp_inA[i] = 0;
    }
  }
  // 调用乘法计算inA*inB
  // todo 改成crossterm
  // if (party == ALICE) {
  //   prod->hadamard_product(dim, exp_inA, exp_inB, outC, bwL+1, bwL+1, bwL +1 + bwL +1,
  //                          signed_arithmetic);
  // } else {
  //   prod->hadamard_product(dim, exp_inA, exp_inB, outC, bwL+1, bwL+1, bwL +1+ bwL +1,
  //                          signed_arithmetic);
  // }
  if (party == ALICE) {
    gp->cross_term(dim, exp_inA, nullptr, outC, bwL+1, bwL+1, bwL +1 + bwL +1);
  } else {
    gp->cross_term(dim, nullptr, exp_inB, outC, bwL+1, bwL+1, bwL +1 + bwL +1);
  }

  uint8_t *msb_0 = new uint8_t[dim];
  for (int i = 0; i < dim; i++) {
    msb_0[i] = 0;
  }

  uint8_t *wrap_outC = new uint8_t[dim];

  ///////////////////////////////wrap_outC 两边都为1时赋值为0/////////////////////////////
  aux->MSB_to_Wrap(outC, msb_0, wrap_outC, dim, 2*bwL+2);
  uint8_t *wrap_outC_alice_send = new uint8_t[dim];
  uint8_t *wrap_outC_bob_recv = new uint8_t[dim];
  if(party == ALICE){
        iopack->io->send_data(wrap_outC, dim * sizeof(uint8_t));
  } else {
    iopack->io->recv_data(wrap_outC_bob_recv, dim * sizeof(uint8_t));
  }
  if (party != ALICE) {
    for (int i = 0; i < dim; i++) {
      if (wrap_outC_bob_recv[i] == 1 && wrap_outC[i] == 1) {
        wrap_outC_bob_recv[i] = 0;
        wrap_outC[i] = 0;
      }
    }
  }
  if (party != ALICE) {
    iopack->io->send_data(wrap_outC_bob_recv, dim * sizeof(uint8_t));
  } else {
    iopack->io->recv_data(wrap_outC_alice_send, dim * sizeof(uint8_t));
    for (int i = 0; i < dim; i++) {
        wrap_outC[i] = wrap_outC_alice_send[i];

    }
  }


  uint64_t mask_bwL_wrap = (1ULL << (2*bwL+2)) - 1;
  // for (int i = 0; i < dim; i++) {

  //    outC[i] = (outC[i]  - wrap_outC[i]* (1ULL << (2*bwL))) & mask_bwL_wrap;
  // }
  // printf("wrap_outC: ");
  // for (int i = 0; i < dim; i++) {
  //   printf("wrap_outC[%d]: %d\n", i, wrap_outC[i]);
  // }
  // printf("\n");
  uint64_t random_numberbob = 0;
  uint64_t random_numberalice = 0;
  uint64_t mask_MWlut_outC = (1ULL << (bwL + bwL+2)) - 1;
  uint64_t *buffer = new uint64_t[dim * 4];
  uint64_t *diff_output = new uint64_t[dim * 4];
  uint64_t *prod_output = new uint64_t[dim * 4];
  // uint64_t *res_output = new uint64_t[dim * 4];
  // 在这里把outC*MW_lut算出来。加个随机数发过去。
  if (party != ALICE) {
    for (int i = 0; i < dim; ++i) {
      for (uint8_t j = 0; j < 4; ++j) {
        // 先计算 MW_lut[j] * outC[i]
        __uint128_t prod = (__uint128_t)MW_lut[j] * (__uint128_t)outC[i];
        prod_output[i * 4 + j] = prod;
        // 再减去 MW_lut[j] 和 outC[i]
        __int128_t diff = (__int128_t)prod - wrap_outC[i]*(__int128_t)MW_lut[j] * (1ULL << (2*bwL+2));
        diff_output[i * 4 + j] = diff;
        // 降低精度（右移）
        uint64_t shift = f_MW + f * 2 - f * 2; // 这里假设fc = f*2，和原代码一致
        uint64_t res = (uint64_t)(((diff + (__int128_t)(1ULL << (shift - 1))) >> shift) & mask_MWlut_outC);
        buffer[i * 4 + j] = (res + random_numberbob) & mask_MWlut_outC;
        // printf("diff: %llu\n", diff);
        // printf("buffer[%d]: %llu\n", i * 4 + j, buffer[i * 4 + j]);
      }
    }
    iopack->io->send_data(buffer, dim * 4 * sizeof(uint64_t));
  } else {
    iopack->io->recv_data(buffer, dim * 4 * sizeof(uint64_t));
  }

  uint64_t *y = new uint64_t[dim];
  if (party == ALICE) {
    uint64_t **spec = new uint64_t *[dim];

    for (int i = 0; i < dim; ++i) {
      spec[i] = new uint64_t[4];
      for (uint8_t j = 0; j < 4; ++j) {
        __uint128_t prod = (__uint128_t)MW_lut[j] * (__uint128_t)outC[i];
        prod_output[i * 4 + j] = prod;
        __int128_t diff = (__int128_t)prod - wrap_outC[i]*(__int128_t)MW_lut[j] * (1ULL << (2*bwL+2));
        diff_output[i * 4 + j] = diff;
        // 降低精度（右移）
        uint64_t shift = f_MW + f * 2 - f * 2; // 这里假设fc = f*2，和原代码一致
        uint64_t res = (uint64_t)(((diff + (__int128_t)(1ULL << (shift - 1))) >> shift) & mask_MWlut_outC);
        spec[i][j] = (res + random_numberalice +buffer[i * 4 + j]) & mask_MWlut_outC;
        // printf("diff: %llu\n", diff);
        // printf("MW_lut[%d]: %llu\n", j, MW_lut[j]);
        // printf("outC[%d]: %llu\n", i, outC[i]);
        // printf("wrap_outC[%d]: %d\n", i, wrap_outC[i]);
        // printf("buffer[%d]: %llu\n", i * 4 + j, buffer[i * 4 + j]);
        // printf("res: %llu\n", res);
        // printf("spec[%d][%d]: %llu\n", i, j, spec[i][j]);
      }
    }
    // for (int i = 0; i < dim; i++) {
    //   for (int j = 0; j < 4; j++) {
    //     printf("spec[%d][%d]: %llu\n", i, j, spec[i][j]);
    //   }
    // }
    aux->lookup_table<uint64_t>(spec, nullptr, nullptr, dim, 2,
                                 bwL + bwL+2);
    for (int i = 0; i < dim; ++i)
      delete[] spec[i];
    delete[] spec;
  } else if (party == BOB) {

    // for (int i = 0; i < dim; i++) {
    //   for (int j = 0; j < 4; j++) {
    //     printf("MW[%d]: %llu\n", i, MW[i]);
    //   }
    // }
    aux->lookup_table<uint64_t>(nullptr, MW, y, dim, 2,bwL + bwL+2);
  }
  // uint64_t mask_res_exp = (1ULL << (f+1+bwL+bwL)) - 1;
  if (party != ALICE) {
    for (int i = 0; i < dim; i++) {
      res_exp[i] =
          (y[i] - random_numberbob) & mask_MWlut_outC; //这里的位宽到底选多少？
    }
  } else {
    for (int i = 0; i < dim; i++) {
      res_exp[i] =
          -random_numberalice & mask_MWlut_outC; //这里的位宽到底选多少？
    }
  }

  size_t comm_end = iopack->io->counter;


    //alice将res_exp发送给bob，bob恢复明文后计算误差
    uint64_t *res_exp_alice = new uint64_t[dim];
    if (party == ALICE) {
        iopack->io->send_data(res_exp, dim * sizeof(uint64_t));
    } else {
        iopack->io->recv_data(res_exp_alice, dim * sizeof(uint64_t));
    }
    double *res_exp_plain = new double[dim];
    double *ideal_exp_plain = new double[dim];
    for (int i = 0; i < dim; i++) {
        res_exp_plain[i] = static_cast<double>((res_exp_alice[i] + res_exp[i]) & mask_MWlut_outC) / std::pow(2, 2*f);
        ideal_exp_plain[i] = std::exp(fix2double(inA[i], inB[i], bwL_input, f_input));
    }
        double ulp_sum = 0.0;
    double ulp_max = 0.0;
  for (int i = 0; i < dim; i++) {
    double ulp = 1.0 / (1 << 10); // 2的12次精度

    // printf("inA[%d]: %llu\n", i, inA[i]);
    // printf("inB[%d]: %llu\n", i, inB[i]);
    // printf("exp_inA[%d]: %llu\n", i, exp_inA[i]);
    // printf("exp_inB[%d]: %llu\n", i, exp_inB[i]);
    // printf("outC[%d]: %llu\n", i, outC[i]);
    // printf("MW[%d]: %llu\n", i, MW[i]);
    // printf("MW_lut[%d]: %llu\n", i, MW_lut[i]);
    // printf("MW_exp[%d]: %llu\n", i, MW_exp[i]);
    // printf("y[%d]: %llu\n", i, y[i]);
    // printf("res_exp[%d]: %llu\n", i, res_exp[i]);
    // printf("res_exp_plain[%d]: %f\n", i, res_exp_plain[i]);
    // printf("ideal_exp_plain[%d]: %f\n", i, ideal_exp_plain[i]);
    double ulp_error = fabs(ideal_exp_plain[i] - res_exp_plain[i]) / ulp;
    // if(ulp_error > 100){
    //   printf("inA[%d]: %llu\n", i, inA[i]);
    //   printf("inB[%d]: %llu\n", i, inB[i]);
    //   printf("MW[%d]: %llu\n", i, MW[i]);
    //   printf("outC[%d]: %llu\n", i, outC[i]);
    //   printf("wrap_outC[%d]: %d\n", i, wrap_outC[i]);
    //   for (int j = 0; j < 4; j++) {
    //     printf("MW_lut[%d]: %llu\n", j, MW_lut[j]);
    //     printf("prod_output[%d][%d]: %llu\n", i, j, prod_output[i * 4 + j]);
    //     printf("diff_output[%d][%d]: %llu\n", i, j, diff_output[i * 4 + j]);
    //     printf("buffer[%d][%d]: %llu\n", i, j, buffer[i * 4 + j]);
    //   }
    //   printf("res_exp[%d]: %llu\n", i, res_exp[i]);
    //   printf("res_exp_plain[%d]: %f\n", i, res_exp_plain[i]);
    //   printf("ideal_exp_plain[%d]: %f\n", i, ideal_exp_plain[i]);

    printf("ULP [%d]: %.2f\n", i, ulp_error);
    // }
    ulp_sum += ulp_error;
    if (ulp_error > ulp_max) ulp_max = ulp_error;
    
  }
  printf("--------------------------------\n");
  if(party != ALICE){
  if (dim > 0) {
    printf("Average ULP: %.2f\n", ulp_sum / dim);
    printf("Max ULP: %.2f\n", ulp_max);
    printf("Total Bytes Sent: %zu bytes\n", (comm_end - comm_start)/dim);
    
  }}


  delete prod;
  delete[] inA; // Delete the variable "inA" to avoid memory leaks
  delete[] inB;
}