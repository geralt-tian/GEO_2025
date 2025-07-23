// anonymous authors

#include "geometric_perspective_protocols.h"
#include <cmath>
#include <cstdint>
#include <math.h>
#include <ctime>
#include <cstdlib>
#include <algorithm>


GeometricPerspectiveProtocols::GeometricPerspectiveProtocols(int party, sci::IOPack *iopack, sci::OTPack *otpack) {
    this->party = party;
    this->iopack = iopack;
    this->otpack = otpack;
    this->aux = new AuxProtocols(party, iopack, otpack);
    this->mill = this->aux->mill;
    this->mill_eq = new MillionaireWithEquality(party, iopack, otpack);
    this->eq = new Equality(party, iopack, otpack);
    this->triple_gen = this->mill->triple_gen;
}

GeometricPerspectiveProtocols::~GeometricPerspectiveProtocols() {
    delete this->aux;
    delete this->mill_eq;
    delete this->eq;
}

void GeometricPerspectiveProtocols::new_truncate(int32_t dim, uint64_t *inA, uint64_t *outB,
                          int32_t shift, int32_t bw) {
    if (shift == 0) {
        memcpy(outB, inA, sizeof(uint64_t) * dim);
        return;
    }
    assert((bw - shift) > 0 && "Truncation shouldn't truncate the full bitwidth");
    assert(bw - shift - 1 >= 0);
    assert(inA != outB);

    uint64_t mask_bw = (bw == 64 ? -1 : ((1ULL << bw) - 1));
    uint64_t mask_shift = (shift == 64 ? -1 : ((1ULL << shift) - 1));
    uint64_t mask_upper =
            ((bw - shift) == 64 ? -1 : ((1ULL << (bw - shift)) - 1));

    uint64_t *inA_orig = new uint64_t[dim];
    uint64_t *inA_lower = new uint64_t[dim];
    uint64_t *inA_upper = new uint64_t[dim];
    uint8_t *wrap_lower = new uint8_t[dim];
    uint8_t *wrap_upper = new uint8_t[dim];
    uint8_t *eq_upper = new uint8_t[dim];
    uint8_t *and_upper = new uint8_t[dim];
    for (int i = 0; i < dim; i++) {
        inA_lower[i] = inA[i] & mask_shift;
        inA_upper[i] = (inA[i] >> shift) & mask_upper;
        if (party == sci::BOB) {
            inA_upper[i] = (mask_upper - inA_upper[i]) & mask_upper;
        }
    }
    uint64_t *inA_prime = new uint64_t[dim];
    uint64_t quarter = 1ULL << (bw - 2);
    if (party == sci::ALICE) {
        for (int i = 0; i < dim; i++) {
            inA_prime[i] = (inA[i] - quarter) & mask_bw;
        }
    } else {
        for (int i = 0; i < dim; i++) {
            inA_prime[i] = inA[i];
        }
    }
    uint64_t *r = new uint64_t[dim];
    uint64_t half = 1ULL << (bw - 1);
    for (int i = 0; i < dim; i++) {
        if (inA_prime[i] >= half) {
            r[i] = 1;
        } else {
            r[i] = 0;
        }
    }
    uint64_t *bit_mul = new uint64_t[dim];
    if (party == sci::ALICE) {
        sci::PRG128 prg;
        uint64_t *data0 = new uint64_t[dim];
        prg.random_data(data0, dim * sizeof(uint64_t));
        otpack->iknp_straight->send_cot(data0, r, dim, shift);
        for (int i = 0; i < dim; i++) {
            bit_mul[i] = ((1ULL << shift) - data0[i]) & mask_shift;
        }
        delete[] data0;
    } else { // party == BOB
        bool *choice = new bool[dim];
        for (int i = 0; i < dim; i++) {
            choice[i] = r[i];
        }
        uint64_t *data = new uint64_t[dim];
        otpack->iknp_straight->recv_cot(data, choice, dim, shift);
        for (int i = 0; i < dim; i++) {
            bit_mul[i] = data[i];
        }
        delete[] choice;
    }
    uint64_t *arith_wrap_lower = new uint64_t[dim];
    for (int i = 0; i < dim; i++) {
      arith_wrap_lower[i] = 0;
    }
    this->aux->wrap_computation(inA_lower, wrap_lower, dim, shift);
    this->aux->B2A(wrap_lower, arith_wrap_lower, dim, bw);
    uint64_t offset = 1ULL << (bw - shift - 2);
    for (int i = 0; i < dim; i++) {
        if (party == sci::ALICE) {
            outB[i] = (((inA_prime[i] >> shift) & mask_upper) + offset + arith_wrap_lower[i] -
                       (1ULL << (bw - shift)) * (bit_mul[i] + 1)) &
                      mask_bw;
        } else {
            outB[i] = (((inA_prime[i] >> shift) & mask_upper) + arith_wrap_lower[i] -
                       (1ULL << (bw - shift)) * bit_mul[i]) &
                      mask_bw;
        }
    }
    delete[] inA_orig;
    delete[] inA_lower;
    delete[] inA_upper;
    delete[] wrap_lower;
    delete[] wrap_upper;
    delete[] eq_upper;
    delete[] and_upper;
    delete[] arith_wrap_lower;

    return;
}

void GeometricPerspectiveProtocols::mw(int32_t dim, uint64_t *input, uint64_t *output, int32_t in_bw, int32_t out_bw) {
    uint64_t mask_in = (in_bw == 64 ? -1 : ((1ULL << in_bw) - 1));
    uint64_t mask_out = (out_bw == 64 ? -1 : ((1ULL << out_bw) - 1));
    uint64_t *input_prime = new uint64_t[dim];
    uint64_t N = 1ULL << (in_bw);
    uint64_t quarter = 1ULL << (in_bw - 2);


    if (party == sci::ALICE) {
        for (int i = 0; i < dim; i++) {
            input_prime[i] = static_cast<uint64_t>(floor((N - ((input[i] - quarter) & mask_in))/(N- 2*quarter))) ;
        }
    } else {
        for (int i = 0; i < dim; i++) {
            input_prime[i] = static_cast<uint64_t>(floor(input[i]/(N- 2*quarter))) ;
        }
    }
    uint64_t *r = new uint64_t[dim];
    uint64_t half = 1ULL << (in_bw - 1);
    if (party == sci::ALICE) {
        for (int i = 0; i < dim; i++) {
            if (input_prime[i] == 0) {
                r[i] = 1;
            } else {
                r[i] = 0;
            }
        }
    }
    else {
        for (int i = 0; i < dim; i++) {
            if (input_prime[i] > 0) {
                r[i] = 1;
            } else {
                r[i] = 0;
            }
        }
    }
    uint64_t *mul = new uint64_t[dim];
    bit_mul(dim, r, mul, out_bw);


    for (int i = 0; i < dim; i++) {
        if (party == sci::ALICE) {
            output[i] = mul[i] + 1;
            if (input[i] < quarter) {
                output[i] = (output[i] - 1) & mask_out;
            }
        } else {
            output[i] = mul[i] & mask_out;
        }
    }


    // for (int i = 0; i < 1000; i++){
    //     printf("input[%d]: %llu\n", i, input[i]);
    //     printf("input_prime[%d]: %llu\n", i, input_prime[i]);
    //     printf("mul[%d]: %llu\n", i, mul[i]);
    //     printf("output[%d]: %llu\n", i, output[i]);
    // }
    // printf ("in_bw: %d\n", in_bw);
    // printf ("N- 2*quarter: %llu\n", N- 2*quarter);
}

void GeometricPerspectiveProtocols::mwwithB(int32_t dim,uint64_t B, uint64_t *input, uint64_t *output, int32_t in_bw, int32_t out_bw) {
uint64_t mask_in = (in_bw == 64 ? -1 : ((1ULL << in_bw) - 1));
uint64_t N = pow(2,in_bw);
uint64_t l_star;
uint64_t mask_out = (out_bw == 64 ? -1 : ((1ULL << out_bw) - 1));
uint64_t *compare_input = new uint64_t[dim];
if (B == N/2){
    l_star = in_bw;
}
else{
    l_star = ceil(log2(floor(N/(N-2*B))));
}

uint8_t *M = new uint8_t[dim];
uint8_t *M_eq = new uint8_t[dim];
uint8_t *delta = new uint8_t[dim];
uint64_t *x0_star = new uint64_t[dim];
uint64_t *M_result = new uint64_t[dim];
if (party == sci::ALICE){
    for (int i = 0; i < dim; i++){
   if (input[i] > B){ //这个=貌似不该取
    delta[i] = 1;
   }
   else{
    delta[i] = 0;
    }
    x0_star[i] = (input[i] - B) & mask_in;
   }
}
else {
    for (int i = 0; i < dim; i++){
        delta[i] = 0;
    }
}

if (B <= (3*N/8)){
    // mw(dim, input, output, in_bw, out_bw);
    return;
}
else{ // B > 3N/8
    if (B == (N/2)){
        // compare(N - x0,x1)
        // uint64_t *compare_input = new uint64_t[dim];
        if (party == sci::ALICE){
            for (int i = 0; i < dim; i++){
                compare_input[i] = N - x0_star[i];
            }
            this->mill_eq->compare_with_eq(M,M_eq,compare_input,dim,l_star,false,true);
        }
        else{
            for (int i = 0; i < dim; i++){
                compare_input[i] = input[i];
            }
            this->mill_eq->compare_with_eq(M,M_eq,compare_input,dim,l_star,false,true);
        }
        
        // delete[] compare_input;
    }
    else{

        
        if (party == sci::ALICE){
            for (int i = 0; i < dim; i++){
                compare_input[i] = std::floor((N - x0_star[i])/(N-2*B));
            }
            this->mill_eq->compare_with_eq(M,M_eq,compare_input,dim,l_star,false);
        }
        else{
            for (int i = 0; i < dim; i++){
                compare_input[i] = std::floor((input[i])/(N-2*B));
            }
            this->mill_eq->compare_with_eq(M,M_eq,compare_input,dim,l_star,false);
        }
        
        // delete[] compare_input;
    }
    
    for (int i = 0; i < dim; i++){
        M[i] = M[i] ^ M_eq[i];
    }

    this->aux->B2A(M, M_result, dim, out_bw);

}
for (int i = 0; i < dim; i++){
    output[i] = (M_result[i] + delta[i]) & mask_out;
    // if (party == sci::ALICE){
    //     printf("input[%d]: %llu\n", i, input[i]);
    //     printf("input_star[%d]: %llu\n", i, x0_star[i]);
    //     printf("compare_input[%d]: %llu\n", i, compare_input[i]);
    //     // printf("input_delta[%d]: %llu\n", i, delta[i]);
    // }
    // else{
    //     printf("input[%d]: %llu\n", i, input[i]);
    //     printf("compare_input[%d]: %llu\n", i, compare_input[i]);
    //     // printf("input_star[%d]: %llu\n", i, x0_star[i]);
    //     // printf("input_delta[%d]: %llu\n", i, delta[i]);
    // }
    // printf("M[%d]: %lld\n", i   , M[i]);    
    // printf("M_result[%d]: %llu\n", i, M_result[i]);
    // printf("delta[%d]: %llu\n", i, delta[i]);
    // printf("output[%d]: %llu\n", i, output[i]);
}
}

void GeometricPerspectiveProtocols::mux_3(int32_t dim, uint64_t *inA, uint64_t *inC, uint64_t *out, int32_t bw) {
    uint64_t mask_bw = (bw == 64 ? -1 : ((1ULL << bw) - 1));
    uint8_t *c0 = new uint8_t[dim];
    uint8_t *c1 = new uint8_t[dim];
    for (int i = 0; i < dim; i++) {
        bool *temp = new bool[2];
        sci::int64_to_bool(temp, inC[i], 2);
        c1[i] = temp[0];
        c0[i] = temp[1];
        delete[] temp;
    }
    uint8_t *carry = new uint8_t[dim];
    if (party == sci::ALICE) {
        uint8_t *dummy = new uint8_t[dim];
        for (int i = 0; i < dim; i++) {
            dummy[i] = 0;
        }
        aux->AND(c1, dummy, carry, dim);
    } else {
        uint8_t *dummy = new uint8_t[dim];
        for (int i = 0; i < dim; i++) {
            dummy[i] = 0;
        }
        aux->AND(dummy, c1, carry, dim);
    }
    uint64_t *temp_inA = new uint64_t[dim];
    for (int i = 0; i < dim; i++) {
        c0[i] = c0[i] ^ carry[i];
        temp_inA[i] = (inA[i] * 2) & mask_bw;
    }
    uint64_t *t1 = new uint64_t[dim];
    uint64_t *t2 = new uint64_t[dim];
    aux->multiplexer(c1, inA, t1, dim, bw, bw);
    aux->multiplexer(c0, temp_inA, t2, dim, bw, bw);
    for (int i = 0; i < dim; i++) {
        out[i] = (t1[i] + t2[i]) & mask_bw;
    }
}

void GeometricPerspectiveProtocols::bit_mul(int32_t dim, uint64_t *input, uint64_t *output, int32_t output_bw) {
    uint64_t mask_out = (output_bw == 64 ? -1 : ((1ULL << output_bw) - 1));
    if (party == sci::ALICE) {
        sci::PRG128 prg;
        uint64_t *data0 = new uint64_t[dim];
        prg.random_data(data0, dim * sizeof(uint64_t));
        otpack->iknp_straight->send_cot(data0, input, dim, output_bw);
        for (int i = 0; i < dim; i++) {
            output[i] = ((1ULL << output_bw) - data0[i]) & mask_out;
        }
        delete[] data0;
    } else { // party == BOB
        bool *choice = new bool[dim];
        for (int i = 0; i < dim; i++) {
            choice[i] = input[i];
        }
        uint64_t *data = new uint64_t[dim];
        otpack->iknp_straight->recv_cot(data, choice, dim, output_bw);
        for (int i = 0; i < dim; i++) {
            output[i] = data[i];
        }
        delete[] choice;
    }
}

void GeometricPerspectiveProtocols::cross_term(int32_t dim, uint64_t *inA, uint64_t *inB, uint64_t *outC,
                                               int32_t bwA, int32_t bwB, int32_t bwC) {
    for (int i = 0; i < dim; i++){
        outC[i] = 0;
    }
    uint64_t mask = (1ULL << bwC) - 1;
    if (party == sci::ALICE) {
        sci::PRG128 prg;
        for (int i = 0; i < bwB; i++) {
            auto *data0 = new uint64_t[dim];
            prg.random_data(data0, dim * sizeof(uint64_t));
            for (int j = 0; j < dim; j++) {
                data0[j] = data0[j] & ((1ULL << (bwC - i)) - 1);
            }
            otpack->iknp_straight->send_cot(data0, inA, dim, bwC - i);
            for (int j = 0; j < dim; j++) {
                outC[j] += (-data0[j] * (1ULL << i));
                outC[j] &= mask;
            }
            delete[] data0;
        }
    } else {
        bool choice[bwB][dim];
        for (int i = 0; i < dim; i++) {
            bool *temp = new bool[bwB];
            sci::int64_to_bool(temp, inB[i], bwB);
            for (int j = 0; j < bwB; j++) {
                choice[j][i] = temp[j];
            }
            delete[] temp;
        }
        for (int i = 0; i < bwB; i++) {
            auto *data = new uint64_t[dim];
            bool *c = new bool[dim];
            for (int j = 0; j < dim; j++) {
                bool *temp = new bool[bwB];
                sci::int64_to_bool(temp, inB[j], bwB);
                c[j] = temp[i];
            }
            otpack->iknp_straight->recv_cot(data, c, dim, bwC - i);
            for (int j = 0; j < dim; j++) {
                outC[j] += (data[j] * (1ULL << i));
                outC[j] &= mask;
            }
            delete[] data;
            delete[] c;
        }
    }
}

void GeometricPerspectiveProtocols::cross_term_reverse(int32_t dim, uint64_t *inA, uint64_t *inB, uint64_t *outC,
                                                       int32_t bwA, int32_t bwB, int32_t bwC) {
    uint64_t mask = (1ULL << bwC) - 1;
    if (party == sci::ALICE) {
        bool choice[bwA][dim];
        for (int i = 0; i < dim; i++) {
            bool *temp = new bool[bwA];
            sci::int64_to_bool(temp, inA[i], bwA);
            for (int j = 0; j < bwA; j++) {
                choice[j][i] = temp[j];
            }
            delete[] temp;
        }
        for (int i = 0; i < bwA; i++) {
            auto *data = new uint64_t[dim];
            bool *c = new bool[dim];
            for (int j = 0; j < dim; j++) {
                bool *temp = new bool[bwA];
                sci::int64_to_bool(temp, inA[j], bwA);
                c[j] = temp[i];
            }
            otpack->iknp_reversed->recv_cot(data, c, dim, bwC - i);
            for (int j = 0; j < dim; j++) {
                outC[j] += ((data[j] * (1ULL << i)) & mask);
                outC[j] &= mask;
            }
            delete[] data;
            delete[] c;
        }
    } else {
        sci::PRG128 prg;
        for (int i = 0; i < bwA; i++) {
            auto *data0 = new uint64_t[dim];
            prg.random_data(data0, dim * sizeof(uint64_t));
            otpack->iknp_reversed->send_cot(data0, inB, dim, bwC - i);
            for (int j = 0; j < dim; j++) {
                outC[j] += ((-data0[j] * (1ULL << i)) & mask);
                outC[j] &= mask;
            }
            delete[] data0;
        }
    }
}

void GeometricPerspectiveProtocols::signed_mul(int32_t dim, uint64_t *inA, uint64_t *inB, uint64_t *outC,
                                                    int32_t bwA, int32_t bwB, int32_t bwC) {
    auto *c = new uint64_t[dim];
    auto *d = new uint64_t[dim];
    for (int i = 0; i < dim; i++) {
        c[i] = 0;
        d[i] = 0;
    }
    if (bwA <= bwB) {
        cross_term_reverse(dim, inA, inB, c, bwA, bwB, bwC);
        cross_term(dim, inB, inA, d, bwB, bwA, bwC);
    } else {
        cross_term(dim, inA, inB, c, bwA, bwB, bwC);
        cross_term_reverse(dim, inB, inA, d, bwB, bwA, bwC);
    }
    auto *m_x = new uint64_t[dim];
    auto *m_y = new uint64_t[dim];
    mw(dim, inA, m_x, bwA, 2);
    mw(dim, inB, m_y, bwB, 2);
    auto *g = new uint64_t[dim];
    auto *h = new uint64_t[dim];
    mux_3(dim, inA, m_y, g, bwA);
    mux_3(dim, inB, m_x, h, bwB);
    uint64_t mask = (1ULL << bwC) - 1;
    for (int i = 0; i < dim; i++) {
        outC[i] = (inA[i] * inB[i] + c[i] + d[i] - (g[i] * (1ULL << bwB)) - (h[i] * (1ULL << bwA))) & mask;
    }
    delete[] c;
    delete[] d;
}

void GeometricPerspectiveProtocols::signed_crossterm(int32_t dim, uint64_t *inA, uint64_t *inB, uint64_t *outC,
                                                    int32_t bwA,  int32_t bwB, int32_t bwC){
    uint64_t maskA = (1ULL << bwA) - 1;
    uint64_t maskB = (1ULL << bwB) - 1;
    uint64_t maskC = (1ULL << bwC) - 1;
    // uint64_t *cross_term_output = new uint64_t[dim];
    if (party == sci::ALICE){
        cross_term(dim, inA, inB, outC, bwA, bwB, bwC);
    }
    else{
        cross_term(dim, inA, inB, outC, bwA, bwB, bwC);
    }
    // for (int i = 0; i < dim; i++){
    //     if (party == sci::ALICE){
    //         outC[i] = (cross_term_output[i] - (1ULL << bwB - 1)*inA[i] + (1ULL << bwA +bwB - 2)) & maskC;
    //     }
    //     else{
    //         outC[i] = (cross_term_output[i] - (1ULL << bwA - 1)*inB[i]) & maskC;
    //     }
    // }

    // for (int i = 0; i < dim; i++){
    //     printf("inA[%d]: %llu\n", i, inA[i]);
    //     printf("inB[%d]: %llu\n", i, inB[i]);
    //     // printf("cross_term_output[%d]: %llu\n", i, cross_term_output[i]);
    //     printf("outC[%d]: %llu\n", i, outC[i]);
    // }
}

void GeometricPerspectiveProtocols::signed_cipher_plainc(int32_t dim, uint64_t *inA, uint64_t *inB, uint64_t *outC,
                                                        int32_t bwA, int32_t bwB, int32_t bwC) {
    uint64_t maskA = (1ULL << bwA) - 1;
    uint64_t maskB = (1ULL << bwB) - 1;
    uint64_t maskC = (1ULL << bwC) - 1;
    uint64_t *inA_p = new uint64_t[dim]; //inA'
    uint8_t *wrap_A_p = new uint8_t[dim];
    uint64_t *inB_p = new uint64_t[dim]; //inB'

    if (party == sci::ALICE) {
        for (int i = 0; i < dim; i++) {
            inA_p[i] = inA[i] + (1ULL << bwA - 1) & maskA;
        }
    } else {
        for (int i = 0; i < dim; i++) {
            inA_p[i] = inA[i] & maskA;
        }
    }

    for (int i = 0; i < dim; i++) {
        inB_p[i] = inB[i] + (1ULL << bwB - 1) & maskB;
    }
    // uint64_t inB_p = new uint64_t[dim]; //inB'
    this->aux->wrap_computation(inA_p, wrap_A_p, dim, bwA);

    if (party == sci::ALICE) {
        for (int i = 0; i < dim; i++) {
            outC[i] = (inA_p[i] * inB_p[i]  - (1ULL << bwA ) * wrap_A_p[i] * inB_p[i]  //x0*y'  -2^m * w(x) *y'
            - (1ULL << bwA - 1) * inB_p[i] //-2^(m-1)*y'
            - (1ULL << bwB - 1) * inA_p[i] //-2^(n-1)*x0
            + (1ULL << bwA + bwB - 1) * wrap_A_p[i] //+2^(m+n-1) * w(x)
            + (1ULL << bwA + bwB - 2)) & maskC; // + 2^(m+n-2)
        }
    }
    else {
        for (int i = 0; i < dim; i++) {
            outC[i] = (inA_p[i] * inB_p[i]  - (1ULL << bwA ) * wrap_A_p[i] * inB_p[i]  //x1*y'  -2^m * w(x) *y'
            - (1ULL << bwB - 1) * inA_p[i] //-2^(n-1)*x1
            + (1ULL << bwA + bwB - 1) * wrap_A_p[i]) & maskC; //+2^(m+n-1) * w(x)
        }
    }   
                                                        
}


void GeometricPerspectiveProtocols::sign_extension(int32_t dim, uint64_t *in, uint64_t *out, int32_t in_bw, int32_t out_bw) {
    if (in_bw == out_bw) {
      memcpy(out, in, sizeof(uint64_t) * dim);
      return;
    }
    assert(in_bw < out_bw);
    uint64_t *c = new uint64_t[dim];
    for (int i = 0; i < dim; i++) {
        c[i] = 0;
    }
    mw(dim, in, c, in_bw, out_bw - in_bw);
    uint64_t M = 1ULL << in_bw;
    uint64_t N = 1ULL << out_bw;
    uint64_t mask = N - 1;
    for (int i = 0; i < dim; i++) {
        out[i] = (in[i] + c[i] * (N - M)) & mask;
    }
}

void GeometricPerspectiveProtocols::truncate_with_one_bit_error(int32_t dim, uint64_t *inA, uint64_t *outB,
                                                                int32_t shift, int32_t bw) {
    if (shift == 0) {
        memcpy(outB, inA, sizeof(uint64_t) * dim);
        return;
    }
    assert((bw - shift) > 0 && "Truncation shouldn't truncate the full bitwidth");
    assert(inA != outB);
    uint64_t mask_bw = (bw == 64 ? -1 : ((1ULL << bw) - 1));
    uint64_t mask_shift = (shift == 64 ? -1 : ((1ULL << shift) - 1));
    uint64_t mask_upper =
            ((bw - shift) == 64 ? -1 : ((1ULL << (bw - shift)) - 1));

    uint64_t *inA_orig = new uint64_t[dim];
    uint64_t *inA_lower = new uint64_t[dim];
    uint64_t *inA_upper = new uint64_t[dim];
    uint8_t *wrap_lower = new uint8_t[dim];
    uint8_t *wrap_upper = new uint8_t[dim];
    uint8_t *eq_upper = new uint8_t[dim];
    uint8_t *and_upper = new uint8_t[dim];
    for (int i = 0; i < dim; i++) {
        inA_lower[i] = inA[i] & mask_shift;
        inA_upper[i] = (inA[i] >> shift) & mask_upper;
        if (party == sci::BOB) {
            inA_upper[i] = (mask_upper - inA_upper[i]) & mask_upper;
        }
    }
    uint64_t *inA_prime = new uint64_t[dim];
    uint64_t quarter = 1ULL << (bw - 2);
    if (party == sci::ALICE) {
        // x - L/4
        for (int i = 0; i < dim; i++) {
            inA_prime[i] = (inA[i] - quarter) & mask_bw;
        }
    } else {
        for (int i = 0; i < dim; i++) {
            inA_prime[i] = inA[i];
        }
    }
    uint64_t *r = new uint64_t[dim];
    uint64_t half = 1ULL << (bw - 1);
    for (int i = 0; i < dim; i++) {
        if (inA_prime[i] >= half) {
            r[i] = 1;
        } else {
            r[i] = 0;
        }
    }
    uint64_t *bit_mul = new uint64_t[dim];
    if (party == sci::ALICE) {
        uint64_t mask_shift = (1ULL << shift) - 1;
        sci::PRG128 prg;
        uint64_t *data0 = new uint64_t[dim];
        prg.random_data(data0, dim * sizeof(uint64_t));
        otpack->iknp_straight->send_cot(data0, r, dim, shift);
        for (int i = 0; i < dim; i++) {
            bit_mul[i] = ((1ULL << shift) - data0[i]) & mask_shift;
        }
        delete[] data0;
    } else { // party == BOB
        bool *choice = new bool[dim];
        for (int i = 0; i < dim; i++) {
            choice[i] = r[i];
        }
        uint64_t *data = new uint64_t[dim];
        otpack->iknp_straight->recv_cot(data, choice, dim, shift);
        for (int i = 0; i < dim; i++) {
            bit_mul[i] = data[i];
        }
        delete[] choice;
    }
    uint64_t offset = 1ULL << (bw - shift - 2);
    for (int i = 0; i < dim; i++) {
        if (party == sci::ALICE) {
            outB[i] = (((inA_prime[i] >> shift) & mask_upper) + offset -
                       (1ULL << (bw - shift)) * (bit_mul[i] + 1)) &
                      mask_bw;
        } else {
            outB[i] = (((inA_prime[i] >> shift) & mask_upper) -
                       (1ULL << (bw - shift)) * bit_mul[i]) &
                      mask_bw;
        }
    }
    delete[] inA_orig;
    delete[] inA_lower;
    delete[] inA_upper;
    delete[] wrap_lower;
    delete[] wrap_upper;
    delete[] eq_upper;
    delete[] and_upper;
    return;
}

void GeometricPerspectiveProtocols::sirnn_unsigned_mul(int32_t dim, uint64_t *inA, uint64_t *inB, uint64_t *outC,
                                                       int32_t bwA, int32_t bwB, int32_t bwC) {
    auto *c = new uint64_t[dim]();
    auto *d = new uint64_t[dim];
    for (int i = 0; i < dim; i++) {
        c[i] = 0;
        d[i] = 0;
    }
    if (bwA <= bwB) {
        cross_term_reverse(dim, inA, inB, c, bwA, bwB, bwC);
        cross_term(dim, inB, inA, d, bwB, bwA, bwC);
    } else {
        cross_term(dim, inA, inB, c, bwA, bwB, bwC);
        cross_term_reverse(dim, inB, inA, d, bwB, bwA, bwC);
    }
    auto *wx = new uint8_t[dim];
    auto *wy = new uint8_t[dim];
    aux->wrap_computation(inA, wx, dim, bwA);
    aux->wrap_computation(inB, wy, dim, bwB);
    auto *h = new uint64_t[dim];
    auto *g = new uint64_t[dim];
    aux->multiplexer(wx, inB, h, dim, bwB, bwB);
    aux->multiplexer(wy, inA, g, dim, bwA, bwA);
    uint64_t mask = (1ULL << bwC) - 1;
    for (int i = 0; i < dim; i++) {
        outC[i] = (inA[i] * inB[i] + c[i] + d[i] - (g[i] * (1ULL << bwB)) - (h[i] * (1ULL << bwA))) & mask;
    }
    delete[] c;
    delete[] d;
}

void GeometricPerspectiveProtocols::msb0_truncation(int32_t dim, uint64_t *inA, uint64_t *outB, int32_t shift, int32_t bw) {
    if (shift == 0) {
        memcpy(outB, inA, sizeof(uint64_t) * dim);
        return;
    }
    assert((bw - shift) > 0 && "Truncation shouldn't truncate the full bitwidth");
    assert(inA != outB);
    uint64_t mask_bw = (bw == 64 ? -1 : ((1ULL << bw) - 1));
    uint64_t mask_shift = (shift == 64 ? -1 : ((1ULL << shift) - 1));
    uint64_t mask_upper =
            ((bw - shift) == 64 ? -1 : ((1ULL << (bw - shift)) - 1));
    uint64_t *inA_orig = new uint64_t[dim];
    uint64_t *inA_lower = new uint64_t[dim];
    uint64_t *inA_upper = new uint64_t[dim];
    uint8_t *wrap_lower = new uint8_t[dim];
    uint8_t *wrap_upper = new uint8_t[dim];
    uint8_t *eq_upper = new uint8_t[dim];
    uint8_t *and_upper = new uint8_t[dim];
    for (int i = 0; i < dim; i++) {
        inA_lower[i] = inA[i] & mask_shift;
        inA_upper[i] = (inA[i] >> shift) & mask_upper;
    }
    uint8_t *r = new uint8_t[dim];
    uint64_t half = 1ULL << (bw - 1);
    for (int i = 0; i < dim; i++) {
        if (inA[i] >= half) {
            r[i] = 1;
        } else {
            r[i] = 0;
        }
    }
    uint8_t *and_result = new uint8_t[dim];
    uint8_t *dummy = new uint8_t[dim];
    for (int i = 0; i < dim; i++) {
        dummy[i] = 0;
    }
    if (party == sci::ALICE) {
        this->aux->AND(r, dummy, and_result, dim);
    } else {
        this->aux->AND(dummy, r, and_result, dim);
    }
    for (int i = 0; i < dim; i++) {
        and_result[i] = and_result[i] ^ r[i];
    }
    uint64_t* mw = new uint64_t[dim];
    this->aux->B2A(and_result, mw, dim, shift);
    uint64_t *arith_wrap_lower = new uint64_t[dim];
    for (int i = 0; i < dim; i++) {
      arith_wrap_lower[i] = 0;
    }
//    this->aux->wrap_computation(inA_lower, wrap_lower, dim, shift);
//    this->aux->B2A(wrap_lower, arith_wrap_lower, dim, bw);
    for (int i = 0; i < dim; i++) {
        if (party == sci::ALICE) {
            outB[i] = (((inA[i] >> shift) & mask_upper) + arith_wrap_lower[i] - (1ULL << (bw - shift)) * (mw[i])) & mask_bw;
        } else {
            outB[i] = (((inA[i] >> shift) & mask_upper) + arith_wrap_lower[i] - (1ULL << (bw - shift)) * mw[i]) & mask_bw;
        }
    }
    delete[] inA_orig;
    delete[] inA_lower;
    delete[] inA_upper;
    delete[] wrap_lower;
    delete[] wrap_upper;
    delete[] eq_upper;
    delete[] and_upper;
    delete[] arith_wrap_lower;
    return;
}

void GeometricPerspectiveProtocols::new_msb0_truncation(int32_t dim, uint64_t *inA, uint64_t *outB, int32_t shift, int32_t bw) {
    if (shift == 0) {
        memcpy(outB, inA, sizeof(uint64_t) * dim);
        return;
    }
    assert((bw - shift) > 0 && "Truncation shouldn't truncate the full bitwidth");
    assert(inA != outB);
    uint64_t mask_bw = (bw == 64 ? -1 : ((1ULL << bw) - 1));
    uint64_t mask_shift = (shift == 64 ? -1 : ((1ULL << shift) - 1));
    uint64_t mask_upper =
            ((bw - shift) == 64 ? -1 : ((1ULL << (bw - shift)) - 1));

    uint64_t *inA_orig = new uint64_t[dim];
    uint64_t *inA_lower = new uint64_t[dim];
    uint64_t *inA_upper = new uint64_t[dim];
    uint8_t *wrap_lower = new uint8_t[dim];
    uint8_t *wrap_upper = new uint8_t[dim];
    uint8_t *eq_upper = new uint8_t[dim];
    uint8_t *and_upper = new uint8_t[dim];
    for (int i = 0; i < dim; i++) {
        inA_lower[i] = inA[i] & mask_shift;
        inA_upper[i] = (inA[i] >> shift) & mask_upper;
        if (party == sci::BOB) {
            inA_upper[i] = (mask_upper - inA_upper[i]) & mask_upper;
        }
    }
    uint64_t *inA_prime = new uint64_t[dim];
    uint64_t quarter = 1ULL << (bw - 2);
    if (party == sci::ALICE) {
        // x - L/4
        for (int i = 0; i < dim; i++) {
            inA_prime[i] = (inA[i] - quarter) & mask_bw;
        }
    } else {
        for (int i = 0; i < dim; i++) {
            inA_prime[i] = inA[i];
        }
    }
    uint64_t *r = new uint64_t[dim];
    uint64_t half = 1ULL << (bw - 1);
    for (int i = 0; i < dim; i++) {
        if (inA_prime[i] >= half) {
            r[i] = 1;
        } else {
            r[i] = 0;
        }
    }
    uint64_t *bit_mul = new uint64_t[dim];
    if (party == sci::ALICE) {
        uint64_t mask_shift = (1ULL << shift) - 1;
        sci::PRG128 prg;
        uint64_t *data0 = new uint64_t[dim];
        prg.random_data(data0, dim * sizeof(uint64_t));
        otpack->iknp_straight->send_cot(data0, r, dim, shift);
        for (int i = 0; i < dim; i++) {
            bit_mul[i] = ((1ULL << shift) - data0[i]) & mask_shift;
        }
        delete[] data0;
    } else { // party == BOB
        bool *choice = new bool[dim];
        for (int i = 0; i < dim; i++) {
            choice[i] = r[i];
        }
        uint64_t *data = new uint64_t[dim];
        otpack->iknp_straight->recv_cot(data, choice, dim, shift);
        for (int i = 0; i < dim; i++) {
            bit_mul[i] = data[i];
        }
        delete[] choice;
    }
    uint64_t offset = 1ULL << (bw - shift - 2);
    for (int i = 0; i < dim; i++) {
        if (party == sci::ALICE) {
            outB[i] = (((inA_prime[i] >> shift) & mask_upper) + offset -
                       (1ULL << (bw - shift)) * (bit_mul[i] + 1)) &
                      mask_bw;
        } else {
            outB[i] = (((inA_prime[i] >> shift) & mask_upper) -
                       (1ULL << (bw - shift)) * bit_mul[i]) &
                      mask_bw;
        }
    }
    delete[] inA_orig;
    delete[] inA_lower;
    delete[] inA_upper;
    delete[] wrap_lower;
    delete[] wrap_upper;
    delete[] eq_upper;
    delete[] and_upper;
    return;
}

void GeometricPerspectiveProtocols::ring_exp(const uint64_t *inA, uint64_t *out, size_t len, uint64_t f,
    uint64_t bwL,uint64_t in_f) {
uint64_t mask_bwL = (bwL == 64 ? -1 : ((1ULL << bwL) - 1));
double pow_f_input = std::pow(2.0, in_f);
double pow_f = std::pow(2.0, f);
for (size_t i = 0; i < len; ++i) {
double x0_real = static_cast<double>(inA[i]) / pow_f_input;
double exp_val = std::exp(x0_real);
uint64_t exp_fixed =
static_cast<uint64_t>(std::round(exp_val * pow_f)) & mask_bwL;
out[i] = exp_fixed;
}
}


void GeometricPerspectiveProtocols::exp(int32_t dim, uint64_t *inA, uint64_t *result, int32_t in_bw,int32_t in_f, int32_t out_bw, int32_t out_f,
     int32_t localexp_bw, int32_t localexp_f, int32_t locallut_bw, int32_t locallut_f) {
    uint64_t *MW = new uint64_t[dim];
    uint64_t *exp_inA = new uint64_t[dim];
    uint64_t *MW_lut = new uint64_t[4];
    uint64_t N_input = 1ULL << in_bw;
    uint64_t MW_B = N_input/2;
    uint64_t mask_locallut_bw = (locallut_bw == 64 ? -1 : ((1ULL << locallut_bw) - 1));
    uint64_t mask_exp_bw = (1ULL << (2*localexp_bw)) - 1;
    uint64_t *outB = new uint64_t[dim];


    ///////////////////////计算MW_lut///////////////////////
    for (int i = 0; i < 3; i++) {
        double pow_f_input = std::pow(2.0, in_f);
        double pow_f = std::pow(2.0, locallut_f);
        double exp_val = std::exp(-i * static_cast<double>(N_input) / pow_f_input); // 可以直接在这里乘exp(-4)
        MW_lut[i] = static_cast<uint64_t>(std::round(exp_val * pow(2, locallut_f))) & mask_locallut_bw;
      }
      MW_lut[3] = 100;

      ///////////////////////计算MW///////////////////////
      mwwithB(dim, MW_B, inA, MW, in_bw, 2); //compute MW

      if (party == sci::ALICE){  //send MW to bob. Bob get whole MW. Alice get 0.
        iopack->io->send_data(MW, dim * sizeof(uint64_t));
        for (int i = 0; i < dim; i++) {
            MW[i] = 0;
        }
      } else {
        uint64_t *recv_MW = new uint64_t[dim];
        iopack->io->recv_data(recv_MW, dim * sizeof(uint64_t));
        for (int i = 0; i < dim; i++) {
          MW[i] = (MW[i] + recv_MW[i]) & ((1ULL << 2) - 1);
        }
        delete[] recv_MW;
      }


    ///////////////////////计算exp_inA///////////////////////
    ring_exp(inA, exp_inA, dim, localexp_f, localexp_bw, in_f);
    uint64_t *zero = new uint64_t[dim];
    for (int i = 0; i < dim; i++) {
        zero[i] = 0;
    }

    ///////////////////////计算cross term///////////////////////
    if (party == sci::ALICE) {
        cross_term(dim, exp_inA, zero, outB, localexp_bw, localexp_bw, localexp_bw + localexp_bw);
      } else {
        cross_term(dim, zero, exp_inA, outB, localexp_bw, localexp_bw, localexp_bw + localexp_bw);
      }
    
      ///////////////////////计算密文*明文///////////////////////
      uint8_t *msb_0 = new uint8_t[dim];
      for (int i = 0; i < dim; i++) {
        msb_0[i] = 0;
      }
    
      uint8_t *wrap_outB = new uint8_t[dim];
    
      ///////////////////////////////wrap_outB 两边都为1时赋值为0/////////////////////////////
      this->aux->MSB_to_Wrap(outB, msb_0, wrap_outB, dim, 2*localexp_bw);
      uint8_t *wrap_outB_alice_send = new uint8_t[dim];
      uint8_t *wrap_outB_bob_recv = new uint8_t[dim];
      if(party == sci::ALICE){
            iopack->io->send_data(wrap_outB, dim * sizeof(uint8_t));
      } else {
        iopack->io->recv_data(wrap_outB_bob_recv, dim * sizeof(uint8_t));
      }
      if (party != sci::ALICE) {
        for (int i = 0; i < dim; i++) {
          if (wrap_outB_bob_recv[i] == 1 && wrap_outB[i] == 1) {
            wrap_outB_bob_recv[i] = 0;
            wrap_outB[i] = 0;
          }
        }
      }
      if (party != sci::ALICE) {
        iopack->io->send_data(wrap_outB_bob_recv, dim * sizeof(uint8_t));
      } else {
        iopack->io->recv_data(wrap_outB_alice_send, dim * sizeof(uint8_t));
        for (int i = 0; i < dim; i++) {
            wrap_outB[i] = wrap_outB_alice_send[i];
    
        }
      }

      uint64_t *buffer = new uint64_t[dim * 4];
      uint64_t random_numberbob = 0;
      uint64_t random_numberalice = 0;
      uint64_t mask_res_exp = (1ULL << (2*localexp_bw - 2*localexp_f +2 +out_f)) - 1;

      if (party != sci::ALICE) {
        for (int i = 0; i < dim; ++i) {
          for (uint8_t j = 0; j < 4; ++j) {
            // 先计算 MW_lut[j] * outB[i]
            __uint128_t prod = (__uint128_t)MW_lut[j] * (__uint128_t)outB[i];
            // 再减去 MW_lut[j] 和 outB[i]
            __int128_t diff = (__int128_t)prod - wrap_outB[i]*(__int128_t)MW_lut[j] * (1ULL << (2*localexp_bw));
            // 降低精度（右移）
            uint64_t shift = locallut_f + localexp_f * 2 - out_f; 
            uint64_t res = (uint64_t)(((diff + (__int128_t)(1ULL << (shift - 1))) >> shift) & mask_res_exp);
            buffer[i * 4 + j] = (res + random_numberbob) & mask_res_exp;
          }
        }
        iopack->io->send_data(buffer, dim * 4 * sizeof(uint64_t));
      } else {
        iopack->io->recv_data(buffer, dim * 4 * sizeof(uint64_t));
      }

      uint64_t *y = new uint64_t[dim];
      if (party == sci::ALICE) {
        uint64_t **spec = new uint64_t *[dim];
    
        for (int i = 0; i < dim; ++i) {
          spec[i] = new uint64_t[4];
          for (uint8_t j = 0; j < 4; ++j) {
            __uint128_t prod = (__uint128_t)MW_lut[j] * (__uint128_t)outB[i];
            __int128_t diff = (__int128_t)prod - wrap_outB[i]*(__int128_t)MW_lut[j] * (1ULL << (2*localexp_bw));
            // 降低精度（右移）
            uint64_t shift = locallut_f + localexp_f * 2 - out_f; // 这里假设fc = f*2，和原代码一致
            uint64_t res = (uint64_t)(((diff + (__int128_t)(1ULL << (shift - 1))) >> shift) & mask_res_exp);
            spec[i][j] = (res + random_numberalice +buffer[i * 4 + j]) & mask_res_exp;
          }
        }
      ///////////////////////计算密文*明文///////////////////////
      this->aux->lookup_table<uint64_t>(spec, nullptr, nullptr, dim, 2,
        2*localexp_bw - 2*localexp_f +2 +out_f);
        for (int i = 0; i < dim; ++i)
        delete[] spec[i];
        delete[] spec;
        } else if (party == sci::BOB) {
        this->aux->lookup_table<uint64_t>(nullptr, MW, y, dim, 2,2*localexp_bw - 2*localexp_f +2 +out_f);
        }


        uint64_t mask_result = (1ULL << (out_bw)) - 1;
        if (party != sci::ALICE) {
            for (int i = 0; i < dim; i++) {
                result[i] =
                  (y[i] - random_numberbob) & mask_result; //这里的位宽到底选多少？
            }
          } else {
            for (int i = 0; i < dim; i++) {
                result[i] =
                  -random_numberalice & mask_result; //这里的位宽到底选多少？
            }
          } 

          delete[] MW;
          delete[] exp_inA;
          delete[] MW_lut;
          delete[] outB;
          delete[] zero;
          delete[] msb_0;
          delete[] wrap_outB;
          delete[] wrap_outB_alice_send;
          delete[] wrap_outB_bob_recv;
          delete[] buffer;
          delete[] y;
}



void GeometricPerspectiveProtocols::exp_nag4(int32_t dim, uint64_t *inA, uint64_t *result,uint64_t *MW, int32_t in_bw,int32_t in_f, int32_t out_bw, int32_t out_f,
    int32_t localexp_bw, int32_t localexp_f, int32_t locallut_bw, int32_t locallut_f) {
//    uint64_t *MW = new uint64_t[dim];
   uint64_t *exp_inA = new uint64_t[dim];
   uint64_t *MW_lut = new uint64_t[4];
   uint64_t N_input = 1ULL << in_bw;
   uint64_t MW_B = N_input/2;
   uint64_t mask_locallut_bw = (locallut_bw == 64 ? -1 : ((1ULL << locallut_bw) - 1));
   uint64_t mask_exp_bw = (1ULL << (2*localexp_bw)) - 1;
   uint64_t *outB = new uint64_t[dim];


   ///////////////////////计算MW_lut///////////////////////
   for (int i = 0; i < 3; i++) {
       double pow_f_input = std::pow(2.0, in_f);
       double pow_f = std::pow(2.0, locallut_f);
       double exp_val = std::exp(-i * static_cast<double>(N_input) / pow_f_input) * std::exp(-4); // 可以直接在这里乘exp(-4)
       MW_lut[i] = static_cast<uint64_t>(std::round(exp_val * pow(2, locallut_f))) & mask_locallut_bw;
     }
     MW_lut[3] = 100;

     ///////////////////////计算MW///////////////////////

    

    size_t comm_start = iopack->io->counter;

    //  mw(dim, inA, MW, MW_bw, 2); //compute MW B=N/4
     
    //  mwwithB(dim, MW_B, inA, MW, in_bw, 2); //compute MW
     size_t comm_end = iopack->io->counter;
    //  printf("MW communication cost: %zu\n", comm_end - comm_start);


     if (party == sci::ALICE){  //send MW to bob. Bob get whole MW. Alice get 0.
       iopack->io->send_data(MW, dim * sizeof(uint64_t));
       for (int i = 0; i < dim; i++) {
           MW[i] = 0;
       }
     } else {
       uint64_t *recv_MW = new uint64_t[dim];
       iopack->io->recv_data(recv_MW, dim * sizeof(uint64_t));
       for (int i = 0; i < dim; i++) {
         MW[i] = (MW[i] + recv_MW[i]) & ((1ULL << 2) - 1);
       }
       delete[] recv_MW;
     }


   ///////////////////////计算exp_inA///////////////////////
   ring_exp(inA, exp_inA, dim, localexp_f, localexp_bw, in_f);
   uint64_t *zero = new uint64_t[dim];
   for (int i = 0; i < dim; i++) {
       zero[i] = 0;
   }

   ///////////////////////计算cross term///////////////////////
   if (party == sci::ALICE) {
       cross_term(dim, exp_inA, zero, outB, localexp_bw, localexp_bw, localexp_bw + localexp_bw);
     } else {
       cross_term(dim, zero, exp_inA, outB, localexp_bw, localexp_bw, localexp_bw + localexp_bw);
     }
   
     ///////////////////////计算密文*明文///////////////////////
     uint8_t *msb_0 = new uint8_t[dim];
     for (int i = 0; i < dim; i++) {
       msb_0[i] = 0;
     }
   
     uint8_t *wrap_outB = new uint8_t[dim];
   
     ///////////////////////////////wrap_outB 两边都为1时赋值为0/////////////////////////////
     this->aux->MSB_to_Wrap(outB, msb_0, wrap_outB, dim, 2*localexp_bw);
     uint8_t *wrap_outB_alice_send = new uint8_t[dim];
     uint8_t *wrap_outB_bob_recv = new uint8_t[dim];
     if(party == sci::ALICE){
           iopack->io->send_data(wrap_outB, dim * sizeof(uint8_t));
     } else {
       iopack->io->recv_data(wrap_outB_bob_recv, dim * sizeof(uint8_t));
     }
     if (party != sci::ALICE) {
       for (int i = 0; i < dim; i++) {
         if (wrap_outB_bob_recv[i] == 1 && wrap_outB[i] == 1) {
           wrap_outB_bob_recv[i] = 0;
           wrap_outB[i] = 0;
         }
       }
     }
     if (party != sci::ALICE) {
       iopack->io->send_data(wrap_outB_bob_recv, dim * sizeof(uint8_t));
     } else {
       iopack->io->recv_data(wrap_outB_alice_send, dim * sizeof(uint8_t));
       for (int i = 0; i < dim; i++) {
           wrap_outB[i] = wrap_outB_alice_send[i];
   
       }
     }

     uint64_t *buffer = new uint64_t[dim * 4];
     uint64_t random_numberbob = 0;
     uint64_t random_numberalice = 0;
     uint64_t mask_res_exp = (1ULL << (2*localexp_bw - 2*localexp_f +2 +out_f)) - 1;

     if (party != sci::ALICE) {
       for (int i = 0; i < dim; ++i) {
         for (uint8_t j = 0; j < 4; ++j) {
           // 先计算 MW_lut[j] * outB[i]
           __uint128_t prod = (__uint128_t)MW_lut[j] * (__uint128_t)outB[i];
           // 再减去 MW_lut[j] 和 outB[i]
           __int128_t diff = (__int128_t)prod - wrap_outB[i]*(__int128_t)MW_lut[j] * (1ULL << (2*localexp_bw));
           // 降低精度（右移）
           uint64_t shift = locallut_f + localexp_f * 2 - out_f; 
           uint64_t res = (uint64_t)(((diff + (__int128_t)(1ULL << (shift - 1))) >> shift) & mask_res_exp);
           buffer[i * 4 + j] = (res + random_numberbob) & mask_res_exp;
         }
       }
       iopack->io->send_data(buffer, dim * 4 * sizeof(uint64_t));
     } else {
       iopack->io->recv_data(buffer, dim * 4 * sizeof(uint64_t));
     }

     uint64_t *y = new uint64_t[dim];
     if (party == sci::ALICE) {
       uint64_t **spec = new uint64_t *[dim];
   
       for (int i = 0; i < dim; ++i) {
         spec[i] = new uint64_t[4];
         for (uint8_t j = 0; j < 4; ++j) {
           __uint128_t prod = (__uint128_t)MW_lut[j] * (__uint128_t)outB[i];
           __int128_t diff = (__int128_t)prod - wrap_outB[i]*(__int128_t)MW_lut[j] * (1ULL << (2*localexp_bw));
           // 降低精度（右移）
           uint64_t shift = locallut_f + localexp_f * 2 - out_f; // 这里假设fc = f*2，和原代码一致
           uint64_t res = (uint64_t)(((diff + (__int128_t)(1ULL << (shift - 1))) >> shift) & mask_res_exp);
           spec[i][j] = (res + random_numberalice +buffer[i * 4 + j]) & mask_res_exp;
         }
       }
     ///////////////////////计算密文*明文 end///////////////////////
     this->aux->lookup_table<uint64_t>(spec, nullptr, nullptr, dim, 2,
       2*localexp_bw - 2*localexp_f +2 +out_f);
       for (int i = 0; i < dim; ++i)
       delete[] spec[i];
       delete[] spec;
       } else if (party == sci::BOB) {
       this->aux->lookup_table<uint64_t>(nullptr, MW, y, dim, 2,2*localexp_bw - 2*localexp_f +2 +out_f);
       }


       uint64_t mask_result = (1ULL << (out_bw)) - 1;
       if (party != sci::ALICE) {
           for (int i = 0; i < dim; i++) {
               result[i] =
                 (y[i] - random_numberbob) & mask_result; //这里的位宽到底选多少？
           }
         } else {
           for (int i = 0; i < dim; i++) {
               result[i] =
                 -random_numberalice & mask_result; //这里的位宽到底选多少？
           }
         } 

         delete[] MW;
         delete[] exp_inA;
         delete[] MW_lut;
         delete[] outB;
         delete[] zero;
         delete[] msb_0;
         delete[] wrap_outB;
         delete[] wrap_outB_alice_send;
         delete[] wrap_outB_bob_recv;
         delete[] buffer;
         delete[] y;
}

void GeometricPerspectiveProtocols::exp_nagx(int32_t dim, uint64_t *inA, uint64_t *result, int32_t in_bw,int32_t in_f, 
    int32_t localexp_bw, int32_t localexp_f, int32_t locallut_bw, int32_t locallut_f) //compute exp( - inA), where inA >0
    {
        uint8_t *drelu = new uint8_t[dim];
        uint64_t *drelu_input = new uint64_t[dim];
        uint64_t mask_in_bw = (in_bw == 64 ? -1 : ((1ULL << in_bw) - 1));
        uint64_t mask_in_fplus3 = (in_f + 3 == 64 ? -1 : ((1ULL << (in_f + 3)) - 1));
        uint64_t pow_f_input = std::pow(2, in_f);
        uint64_t *inA_expinput = new uint64_t[dim];
        uint64_t *exp_result = new uint64_t[dim];
        uint64_t *MW = new uint64_t[dim];


        ///////////////////////计算drelu_input///////////////////////
        if (party == sci::ALICE) {
        for (int i = 0; i < dim; i++) {
            drelu_input[i] = (-inA[i] + 4*pow_f_input) & mask_in_bw;
          }
        } else {
          for (int i = 0; i < dim; i++) {
            drelu_input[i] = (-inA[i] + 4*pow_f_input) & mask_in_bw;
          }
        }





        aux->MSB(drelu_input, drelu, dim, in_bw);

          if(party == sci::ALICE){
          for (int i = 0; i < dim; i++) {
            drelu[i] = 1- drelu[i]; 
          }
          }

          ///////////// 
          if (party == sci::ALICE) {
          for (int i = 0; i < dim; i++) {
            inA_expinput[i] = (-inA[i] + 2*pow_f_input) & mask_in_fplus3;
          }
        } else {
          for (int i = 0; i < dim; i++) {
            inA_expinput[i] = (-inA[i] + 2*pow_f_input) & mask_in_fplus3;
          }
        }

        ///////////////计算exp( - inA) * exp(-4)///////////////////////
        size_t comm_start_exp = iopack->get_comm();
        

        uint64_t *MW_input = new uint64_t[dim];
        if (party == sci::ALICE) {
            for (int i = 0; i < dim; i++) {
                MW_input[i] = (-inA[i] + 2*pow_f_input) & mask_in_bw;
            }
          } else {
            for (int i = 0; i < dim; i++) {
                MW_input[i] = (-inA[i] + 2*pow_f_input) & mask_in_bw;
            }
          }
  
        mw(dim, MW_input, MW, in_bw, 2);
          exp_nag4(dim, inA_expinput,exp_result, MW , in_f + 3, in_f, in_f + 3, in_f, localexp_bw, localexp_f, locallut_bw, locallut_f);
          size_t comm_end_exp = iopack->get_comm();
        std::cout << "exp nag4 comm: " << comm_end_exp - comm_start_exp << std::endl;
        
          this->aux->multiplexer(drelu, exp_result, result, dim, in_bw, in_bw);

        // for (int i = 0; i < 1000; i++) {
        //     printf("inA[%d]: %llu\n", i, inA[i]);
        //     printf("drelu_input[%d]: %d\n", i, drelu_input[i]);
        //     printf("drelu[%d]: %d\n", i, drelu[i]);
        //     printf("exp_result[%d]: %llu\n", i, exp_result[i]);
        //     printf("result[%d]: %llu\n", i, result[i]);
        //     printf("\n");
        // }
        // printf("mask_in_bw: %llu\n", mask_in_bw);
    }




    // vector_scalar multiplication bwA>=bwB
    void GeometricPerspectiveProtocols::vector_bit_mul(int32_t dim, uint64_t *inA, uint8_t ot_choice, uint64_t *outC,
    int32_t bwA) {
        uint64_t mask_A = (1ULL << bwA) - 1;
        // uint64_t mask_C = (1ULL << bwC) - 1;

        // 计算紧密打包后的总block数量 (在所有parties中都需要)
        int32_t values_per_block = 128 / bwA;
        int32_t total_blocks = (dim + values_per_block - 1) / values_per_block;

        // Alice端：生成数据和密钥
        sci::block128* encoded_blocks = nullptr;
        sci::block128* encoded_blocks_f = nullptr;
        
        // 密钥只由Alice生成，Bob不知道
        // Bob永远不会直接访问这些密钥，只能通过OT获得其中一个
        unsigned char key_bytes[16];
        unsigned char key_bytes_f[16];
        
        if (party == sci::ALICE) {
            // 只有Alice生成密钥
            unsigned char alice_key_0[16] = {0,9,0,0,0,0,0,0,0,0,0,0,0,0,9,9};
            unsigned char alice_key_1[16] = {1,1,1,'a',0,0,1,0,0,0,0,5,0,0,4,0};
            memcpy(key_bytes, alice_key_0, 16);
            memcpy(key_bytes_f, alice_key_1, 16);
            
            uint64_t *A_random = new uint64_t[dim];
            uint64_t *random_number = new uint64_t[dim];
            
            // 设置固定种子确保每次运行都生成相同的随机数序列
            srand(1);
            for (int i = 0; i < dim; i++) {
                random_number[i] = (rand() % (mask_A + 1));
                outC[i] = (-random_number[i]) & mask_A;
            }
            //将inA和random_number相加
            for (int i = 0; i < dim; i++) {
                A_random[i] = (inA[i] + random_number[i]) & mask_A;
            }
            for (int i = 0; i < 100; i++) {
                printf("random_number[%d]: %llu\n", i, random_number[i]); 
                printf("A_random[%d]: %llu\n", i, A_random[i]);
            }
            
            
            //AES key generation
            sci::AESNI_KEY encrypt_key, encrypt_key_f;
            sci::AESNI_set_encrypt_key(&encrypt_key, key_bytes, 16);
            sci::AESNI_set_encrypt_key(&encrypt_key_f, key_bytes_f, 16);

            // 编码和加密
            encoded_blocks = encode_A_random_to_blocks(A_random, dim, bwA);
            encoded_blocks_f = encode_A_random_to_blocks(random_number, dim, bwA);
            
            // 打印调试信息
            uint64_t block_data[2];
            _mm_storeu_si128((__m128i*)block_data, encoded_blocks[0]);
            printf("encoded_blocks[0]: 0x%016llx%016llx (contains %d values)\n", 
                   block_data[1], block_data[0], values_per_block);

            sci::AESNI_ecb_encrypt_blks(encoded_blocks, total_blocks, &encrypt_key);
            sci::AESNI_ecb_encrypt_blks(encoded_blocks_f, total_blocks, &encrypt_key_f);

            _mm_storeu_si128((__m128i*)block_data, encoded_blocks[0]);
            printf("encrypt encoded_blocks[0]: 0x%016llx%016llx\n", block_data[1], block_data[0]);
            
            delete[] A_random;
            delete[] random_number;
        }

        // Alice发送加密数据，Bob接收
        sci::block128* received_encrypt_blocks = new sci::block128[total_blocks];
        sci::block128* received_encrypt_blocks_f = new sci::block128[total_blocks];
        
        if (party == sci::ALICE) {
            iopack->io->send_data(encoded_blocks, total_blocks * sizeof(sci::block128));
            iopack->io->send_data(encoded_blocks_f, total_blocks * sizeof(sci::block128));
        } else {
            iopack->io->recv_data(received_encrypt_blocks, total_blocks * sizeof(sci::block128));
            iopack->io->recv_data(received_encrypt_blocks_f, total_blocks * sizeof(sci::block128));
        }

        // OT传输密钥 - 使用多组OT数组提高效率
        const int num_ots = 8;  // 可以根据需要调整OT数量
        sci::block128 key_block_0, key_block_1;
        sci::block128 received_key_block[num_ots];
        bool ot_choice_bool[num_ots];
        
        if (party == sci::ALICE) {
            // Alice: 将16字节密钥转换为block128
            key_block_0 = _mm_loadu_si128((__m128i*)key_bytes_f);
            key_block_1 = _mm_loadu_si128((__m128i*)key_bytes);
            
            // 准备多组OT数据
            sci::block128 ot_data0[num_ots];
            sci::block128 ot_data1[num_ots];
            
            for (int i = 0; i < num_ots; i++) {
                ot_data0[i] = key_block_0;  // 选择0的密钥
                ot_data1[i] = key_block_1;  // 选择1的密钥
            }
            
            // 执行多组OT
            otpack->iknp_straight->send(ot_data0, ot_data1, num_ots);
            printf("Alice sent %d OTs with 128-bit keys\n", num_ots);
            for (int i = 0; i < num_ots; i++) {
                uint64_t data0_parts[2], data1_parts[2];
                _mm_storeu_si128((__m128i*)data0_parts, ot_data0[i]);
                _mm_storeu_si128((__m128i*)data1_parts, ot_data1[i]);
                printf("ot_data0[%d]: 0x%016llx%016llx\n", i, data0_parts[1], data0_parts[0]);
                printf("ot_data1[%d]: 0x%016llx%016llx\n", i, data1_parts[1], data1_parts[0]);
            }
        } else {
            // Bob: 准备多组选择位
            for (int i = 0; i < num_ots; i++) {
                ot_choice_bool[i] = (ot_choice == 1);  // 所有OT使用相同选择位
            }
            
            // 接收多组OT结果
            otpack->iknp_straight->recv(received_key_block, ot_choice_bool, num_ots);
            printf("Bob received %d OTs, using first result\n", num_ots);
            for (int i = 0; i < num_ots; i++) {
                uint64_t received_parts[2];
                _mm_storeu_si128((__m128i*)received_parts, received_key_block[i]);
                printf("received_key_block[%d]: 0x%016llx%016llx\n", i, received_parts[1], received_parts[0]);
            }
        }

        // 只有Bob进行解密
        if (party == sci::BOB) {
            // 将接收到的block128密钥转换为AES可用格式
            unsigned char received_key[16];
            _mm_storeu_si128((__m128i*)received_key, received_key_block[0]);
            
            sci::AESNI_KEY received_decrypt_key;
            sci::AESNI_set_decrypt_key(&received_decrypt_key, received_key, 16);
            
            if (ot_choice == 1) {
                sci::AESNI_ecb_decrypt_blks(received_encrypt_blocks, total_blocks, &received_decrypt_key);
                uint64_t* decoded_data = decode_blocks_to_A_random(received_encrypt_blocks, dim, bwA);
                // decoded_data 现在包含 A_random = inA + random_number
                memcpy(outC, decoded_data, dim * sizeof(uint64_t));
                delete[] decoded_data;
            } else if (ot_choice == 0) {
                sci::AESNI_ecb_decrypt_blks(received_encrypt_blocks_f, total_blocks, &received_decrypt_key);
                uint64_t* decoded_data = decode_blocks_to_A_random(received_encrypt_blocks_f, dim, bwA);
                // decoded_data 现在包含 random_number
                memcpy(outC, decoded_data, dim * sizeof(uint64_t));
                delete[] decoded_data;
            }
            for (int i = 0; i < 100; i++) {
                printf("outC[%d]: %llu \n", i, outC[i]); 
            }
            // printf("outC[0]: %llu \n",outC[0]); 
        }

        // 清理内存
        if (party == sci::ALICE) {
            delete[] encoded_blocks;
            delete[] encoded_blocks_f;
        }
        delete[] received_encrypt_blocks;
        delete[] received_encrypt_blocks_f;
    }



    // 🔥 新增：将A_random数组转换为block128数组（选项1：紧密打包多个值）
    sci::block128* GeometricPerspectiveProtocols::encode_A_random_to_blocks(
        const uint64_t* A_random, int32_t dim, int32_t bwA) {
        
        // 计算每个block128能放多少个值
        int32_t values_per_block = 128 / bwA;  // 整数除法
        int32_t total_blocks = (dim + values_per_block - 1) / values_per_block;  // 向上取整
        
        // 分配block128数组
        sci::block128* block_array = new sci::block128[total_blocks];
        
        // 处理每个block
        for (int32_t block_idx = 0; block_idx < total_blocks; block_idx++) {
            // 初始化当前block为0
            uint64_t block_data_low = 0;
            uint64_t block_data_high = 0;
            
            // 在当前block中打包多个值
            for (int32_t i = 0; i < values_per_block; i++) {
                int32_t value_idx = block_idx * values_per_block + i;
                
                if (value_idx < dim) {
                    uint64_t value = A_random[value_idx];
                    int32_t bit_offset = i * bwA;
                    
                    if (bit_offset < 64) {
                        // 值放在低64位
                        if (bit_offset + bwA <= 64) {
                            // 完全在低64位内
                            block_data_low |= (value << bit_offset);
                        } else {
                            // 跨越低64位和高64位
                            int32_t bits_in_low = 64 - bit_offset;
                            int32_t bits_in_high = bwA - bits_in_low;
                            block_data_low |= (value << bit_offset);
                            block_data_high |= (value >> bits_in_low);
                        }
                    } else {
                        // 值完全在高64位
                        block_data_high |= (value << (bit_offset - 64));
                    }
                }
            }
            
            // 创建block128
            block_array[block_idx] = sci::makeBlock128(block_data_high, block_data_low);
        }
        
        return block_array;
    }
    
    // 🔥 新增：从block128数组解码回A_random数组（选项1：紧密打包多个值）
    uint64_t* GeometricPerspectiveProtocols::decode_blocks_to_A_random(
        const sci::block128* block_array, int32_t dim, int32_t bwA) {
        
        // 计算每个block128能放多少个值
        int32_t values_per_block = 128 / bwA;  // 整数除法
        int32_t total_blocks = (dim + values_per_block - 1) / values_per_block;  // 向上取整
        
        // 分配A_random数组
        uint64_t* A_random = new uint64_t[dim];
        
        // 从每个block中提取多个值
        for (int32_t block_idx = 0; block_idx < total_blocks; block_idx++) {
            // 提取当前block的数据
            uint64_t block_data[2];
            _mm_storeu_si128((__m128i*)block_data, block_array[block_idx]);
            uint64_t block_data_low = block_data[0];
            uint64_t block_data_high = block_data[1];
            
            // 从当前block中提取多个值
            for (int32_t i = 0; i < values_per_block; i++) {
                int32_t value_idx = block_idx * values_per_block + i;
                
                if (value_idx < dim) {
                    int32_t bit_offset = i * bwA;
                    uint64_t value = 0;
                    
                    if (bit_offset < 64) {
                        // 值从低64位开始
                        if (bit_offset + bwA <= 64) {
                            // 完全在低64位内
                            uint64_t mask = (bwA == 64) ? 0xFFFFFFFFFFFFFFFFULL : ((1ULL << bwA) - 1);
                            value = (block_data_low >> bit_offset) & mask;
                        } else {
                            // 跨越低64位和高64位
                            int32_t bits_in_low = 64 - bit_offset;
                            int32_t bits_in_high = bwA - bits_in_low;
                            uint64_t low_part = (block_data_low >> bit_offset);
                            uint64_t high_mask = (bits_in_high == 64) ? 0xFFFFFFFFFFFFFFFFULL : ((1ULL << bits_in_high) - 1);
                            uint64_t high_part = (block_data_high & high_mask) << bits_in_low;
                            value = low_part | high_part;
                        }
                    } else {
                        // 值完全在高64位
                        uint64_t mask = (bwA == 64) ? 0xFFFFFFFFFFFFFFFFULL : ((1ULL << bwA) - 1);
                        value = (block_data_high >> (bit_offset - 64)) & mask;
                    }
                    
                    A_random[value_idx] = value;
                }
            }
        }
        
        return A_random;
    }
    




// 矩阵是m*n大小，vector是m长。矩阵的i行都需要与vector中的i项做crossterm，但是为了降低轮数，mn矩阵和m*1的vector要一起做crossterm。
// crossterm中inA和inB不是share值，但outC是share值（输出类似cot）。
    void GeometricPerspectiveProtocols::matrix_vector_crossterm(int32_t m,int32_t n, uint64_t *inA, uint64_t *inB, uint64_t *outC,
        int32_t bwA, int32_t bwB) {
            int32_t bwC = bwA + bwB;

    //step 0 alice生成一组m*n的随机数 random。再将随机数与inA相加得到inA_random。



    //step 1 将random 和inA_random encode成block128格式得到random_encode128和inA_random_encode128(因为AES加密需要使用block128格式)


    //step 2 生成aes加密密钥key_random和key_inA_random
    //step 2.1 将random_encode128和inA_random_encode128使用aes加密成密文enc_random和enc_inA_random

    //step 3 将enc_random和enc_inA_random发送给bob

    //step 4 alice调用otpack->iknp_straight->send，输入m个（key_random，key_inA_random）pair
    //step 4.1 bob调用otpack->iknp_straight->recv，输入m个choice bit

    //step 5 bob根据输入的choice bit，使用拿到的密钥解密拿到的密文enc_random和enc_inA_random，得到明文random和inA_random

    //step 6 将random和inA_random decode成uint64_t格式得到random_decode和inA_random_decode

    //step 7 alice和bob执行outC +=  参考sirnn crossterm


    // 循环上面步骤，bob更改选择bit




    
    }
// 矩阵是m*n大小，vector是m长。矩阵的i行都需要与vector中的i项做乘法，但是为了降低轮数，mn矩阵和m*1的vector要一起做乘法。
// inA和inB都是share值，outC是也是share值
// 写mul之前要先写matrix_vectorcrossterm和crossterm reverse
    void GeometricPerspectiveProtocols::matrix_vector_mul(int32_t m,int32_t n, uint64_t *inA, uint64_t *inB, uint64_t *outC,
    int32_t bwA, int32_t bwB, int32_t bwC) {


    }

