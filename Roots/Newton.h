#ifndef NEWTON_H_INCLUDED
#define NEWTON_H_INCLUDED

#include <cstdint>
#include <bit>
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
// 
static const uint64_t twoPowerRoot(uint64_t n, uint64_t k, uint64_t v, uint64_t L_inv) {
    // bool debug = false;
    // if (n == 81 && k == 4){
    //     debug = true;
    //     // 81 = 1010001 = 1 + 5*2^4
    //     // v = 2
    // }
    
    uint64_t digits_of_precision = 1+(63/k);

    uint64_t x_two_val = std::countr_zero(n >> v); // (Really (n-1) >> v)
    // x = 1 + 4 * 5*2^2. x_two_val = 2.
    // 

    int g = digits_of_precision+2- 2*x_two_val - v;
    g = std::max(g, 0);
    for(int i = 0; i < g; ++i) {
            n = n*n;
    }
    

    uint64_t x = L_inv * (n>>(v+g)); // (Really (n-1) >> v+g)
    // Now calculate exp(x) to dop digits.
    auto rv = twoadic_exp_precision(x, digits_of_precision+1);

    // if (debug) {
    //     std::string s = "";
    //     auto tmp = rv;
    //     while (tmp!=0) {
    //         if(tmp %2 == 0) {
    //             s = "0" + s;
    //         } else {
    //             s = "1" + s;
    //         }
    //         tmp/=2;
    //     }
    //     std::cout << s << std::endl;
    // }

    auto sign_check = 1ull << (digits_of_precision);
 
    if (rv & sign_check) {
        rv = -rv;
    }

    return rv & (sign_check-1);


}

// Compute kth root of n. 
static const uint64_t twoAdicKthRoot(uint64_t n, uint64_t k) {

    uint64_t k_two_val = std::countr_zero(k);

    uint64_t k_odd_inv = mod_inverse_h(k >> k_two_val);

    uint64_t n_two_val = std::countr_zero(n);
    uint64_t sol_two_val = (n_two_val >> k_two_val)*k_odd_inv;

    uint64_t n_odd = n >> n_two_val;

    if (k_two_val == 0) {
        // k is odd; implicitly handles n = 3 (mod 4). 
        return combined_pow(n_odd, k_odd_inv) << sol_two_val;
    } else {
        
        return twoPowerRoot(n_odd, k, k_two_val, k_odd_inv) << sol_two_val;
    }

    
    uint64_t final_mod = 1ull << ((64/k)+1);

    auto rv = (twoadic_exp(k_odd_inv*(twoadic_log(n_odd)>>k_two_val)));
    
    if (rv & final_mod) {
        rv = -rv;
    }
    return(rv & (final_mod-1))<<sol_two_val;

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