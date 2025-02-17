#ifndef NEWTON_H_INCLUDED
#define NEWTON_H_INCLUDED

#include <cstdint>
#include <bit>
#include <algorithm>
#include "../ImprovedApproach/LookupApproach.h"
#include "../ImprovedApproach/CombinedApproach.h"


// Inverse formula from   
// https://marc-b-reynolds.github.io/math/2017/09/18/ModInverse.html#fnref:hurchalla
// https://arxiv.org/abs/2204.04342
static inline uint64_t mod_inverse_h(uint64_t a)
{
  uint64_t x = (3*a)^2; 
  uint64_t y  = 1 - a*x;
  x = x*(1+y); y *= y;
  x = x*(1+y); y *= y;
  x = x*(1+y); y *= y;
  x = x*(1+y);
  return x;
}

// Compute the k = L*2^vth root of n.
// Assume v > 0, thus n = 1 (mod 2^{v+2})
// For n = 64, answer should have no more than 1+64/k digits.
static const uint64_t EulerEvenRoot(uint64_t n, uint64_t dop, uint64_t v, uint64_t L_inv) {
    uint64_t x_two_val = std::countr_zero(n >> v); // (Really (n-1) >> v)
    
    // Each square gains a net one bit of accuracy. 
    int addl_squares = dop+2- 2*x_two_val - v;
    addl_squares = std::max(addl_squares, 0);
    
    // This could perhaps be done more efficiently by using the adic_pow function 
    // modified to only compute the exact bits of precision needed.
    for(int i = 0; i < addl_squares; ++i) {
            n = n*n;
    }

    uint64_t x = L_inv * (n>>(v+addl_squares)); // (Think (n-1) >> v+addl_powers)
    // Now calculate exp(x) to dop digits.
    return twoadic_exp_precision(x, dop+1);
}

// Compute kth root of n. 
static const uint64_t twoAdicKthRoot(uint64_t n, uint64_t k) {
    // Digits of precision; a kth root will never have more than this many bits
    uint64_t dop = 1+(63/k);

    uint64_t k_two_val = std::countr_zero(k);
    uint64_t k_odd_inv = mod_inverse_h(k >> k_two_val);

    uint64_t n_two_val = std::countr_zero(n);
    uint64_t sol_two_val = (n_two_val >> k_two_val)*k_odd_inv;

    uint64_t n_odd = n >> n_two_val;

    if (k_two_val == 0) {
        // k is odd; implicitly handles n = 3 (mod 4). 
        // Note this could be further optimized to use dop, as the 
        // combined_pow function assumes all 64 bits of precision.  
        return combined_pow(n_odd, k_odd_inv) << sol_two_val;
    }

    // The commented out function does the more straightforward "Root Equation"
    // approach of computing the log, dividing, then computing the exp. But it 
    // does not exploit information about digits of precision, which partially 
    // explains the slower running time compared to the EulerEvenRoot function. 
    //
    // auto rv = twoadic_exp(k_odd_inv*(twoadic_log(n_odd)>>k_two_val));
    auto rv = EulerEvenRoot(n_odd, dop, k_two_val, k_odd_inv);
    
    uint64_t sign_check = 1ull << dop;
    if (rv & sign_check) {
        rv = -rv;
    }
    return(rv & (sign_check-1))<<sol_two_val;
}




// Compute kth root of n. 
static const uint64_t NewtonKthRoot(uint64_t n, uint64_t k) {

    uint64_t y = 1ull << ((64-std::countl_zero(n))/k); 
    uint64_t y_prev = 0;
    
    while (y_prev != y) {
        y_prev = y; 
        y = ((k-1)*y + n/trad_pow(y, k-1))/k;
    }

    return y;
}

#endif // NEWTON_H_INCLUDED