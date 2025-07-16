// anonymous authors

#include "geometric_perspective_protocols.h"
#include <cmath>
#include <math.h>

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
          exp_nag4(dim, inA_expinput,exp_result, MW , in_f + 3, in_f, in_f + 3, in_f, localexp_bw, localexp_f, locallut_bw, locallut_f);

        
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