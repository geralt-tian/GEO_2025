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

// int bwL = 22;
// uint64_t N = pow(2, bwL);
// uint64_t mask_bwL = (bwL == 64 ? -1 : ((1ULL << bwL) - 1));
bool accumulate = true;
bool precomputed_MSBs = false;
MultMode mode = MultMode::None;

// uint64_t f = 12;

uint64_t f_MW = 28;
uint64_t bwL_MW = f_MW+2;
uint64_t N_f_MW = pow(2, f_MW);
uint64_t N_MW = pow(2, bwL_MW);
uint64_t mask_N_MW = (bwL_MW == 64 ? -1 : ((1ULL << bwL_MW) - 1));
// uint64_t pow_f = pow(2, f);

uint64_t f_input = 12;
int bwL_input = f_input + 9;
uint64_t mask_bwL_input = (bwL_input == 64 ? -1 : ((1ULL << bwL_input) - 1));
uint64_t N_input = pow(2, bwL_input);
uint64_t pow_f_input = pow(2, f_input);
uint64_t mask_in_bw = mask_bwL_input;



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
// int dim = 16384;
int dim = 1ULL<<18;

// Declare missing variables  

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

void ring_sin(const uint64_t *inA, uint64_t *out, size_t len, uint64_t f_input,uint64_t f,
              uint64_t mask_bwL) {
  double pow_f = std::pow(2.0, f);
  double pow_f_input = std::pow(2.0, f_input);
  for (size_t i = 0; i < len; ++i) {
    double x0_real = static_cast<double>(inA[i]) / pow_f_input;
    double sin_val = std::sin(x0_real);
    uint64_t sin_fixed =
        static_cast<uint64_t>(std::round(sin_val * pow_f)) & mask_bwL;
    out[i] = sin_fixed;
  }
}

void ring_cos(const uint64_t *inA, uint64_t *out, size_t len, uint64_t f_input,uint64_t f,
              uint64_t mask_bwL) {
  double pow_f = std::pow(2.0, f);
  double pow_f_input = std::pow(2.0, f_input);
  for (size_t i = 0; i < len; ++i) {
    double x0_real = static_cast<double>(inA[i]) / pow_f_input;
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
  uint64_t *MW_plain = new uint64_t[dim];

  uint64_t bwL_input_div = 21;
  uint64_t B = 0.5*2*(1ULL << bwL_input_div)/4;

  uint64_t ddd= 10;
  uint64_t bwL_d = static_cast<uint64_t>(ceil(log2(ddd)));
  size_t comm_start = iopack->io->counter;
  auto start_time = chrono::high_resolution_clock::now();


  uint seed = 10;
  for (int i = 0; i < dim; i++) {
    inA[i] = (rand_r(&seed) & mask_bwL_input) % (B / 2);
    inB[i] = (rand_r(&seed) & mask_bwL_input) % (B / 2); 
  }

  if (B == N_input/4) {
    if (party == ALICE) {
      gp->mw(dim, inA, MW, bwL_input_div, 2);
    } else {
      gp->mw(dim, inB, MW, bwL_input_div, 2);
    }
  } else {
    if (party == ALICE) {
      gp->mwwithB(dim, B, inA, MW, bwL_input_div, 2);
    } else {
      gp->mwwithB(dim, B, inB, MW, bwL_input_div, 2);
    }
  }

  uint64_t *T_add_T = new uint64_t[dim];
  uint64_t *T_add_T_reduce = new uint64_t[dim];
  for (int i = 0; i < dim; i++) {
    T_add_T[i] = MW[i] + MW[i];
    T_add_T_reduce[i] = T_add_T[i] >> (2*f_MW+1);
  }
  uint64_t** sc_cs_lut_spec = new uint64_t*[dim];

  for(int i = 0; i < dim; i++) {
    sc_cs_lut_spec[i] = new uint64_t[4]; // 假设每行有4个元素，根据MW的输出
    
  }
  
  if (party == ALICE) {
    aux->lookup_table<uint64_t>(sc_cs_lut_spec, nullptr, nullptr, dim, 2, bwL_input_div);
    aux->lookup_table<uint64_t>(sc_cs_lut_spec, nullptr, nullptr, dim, 2, bwL_d+1);
    // aux->lookup_table<uint64_t>(cc_ss_lut_spec, nullptr, nullptr, dim, 2, 2*bwL_t+bwL_T);
  }
  else {
    aux->lookup_table<uint64_t>(nullptr, MW, T_add_T, dim, 2,bwL_input_div);
    aux->lookup_table<uint64_t>(nullptr, MW, T_add_T, dim, 2,bwL_d+1);
  }

  uint8_t *drelu = new uint8_t[dim];
  aux->MSB(T_add_T_reduce, drelu, dim, bwL_d+1);

  uint64_t *drelu_input = new uint64_t[dim];
  for (int i = 0; i < dim; i++) {
    drelu_input[i] = (inA[i] + 4 * pow_f_input) & mask_in_bw;
  }

  uint64_t *drelu_input_reduce = new uint64_t[dim];
  aux->B2A((uint8_t*)drelu_input, drelu_input_reduce, dim, bwL_input_div);

  size_t comm_end = iopack->io->counter;
  size_t comm_bytes = comm_end - comm_start;
  auto end_time = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);

  printf("Communication: %zu bytes\n", comm_bytes);
  printf("Time: %ld ms\n", duration.count());
  

  delete prod;
  delete[] inA; // Delete the variable "inA" to avoid memory leaks
  delete[] inB;
}