#ifndef LOOKUPAPPROACH_H_INCLUDED
#define LOOKUPAPPROACH_H_INCLUDED

// Private:

// Some magic numbers you might come across:
// 3^-1 = 12297829382473034411ull
// 5^-1 = 14757395258967641293ull


// This array satisfies logarray[i] == log_2(2^i + 1) whenever both sides
// make sense. Note that log_2(2) and log_2(3) are undefined.
static const uint64_t logarray[] = {0,0,6713115954038056572ull,
    6165135171829223912ull, 6071351533721251728ull, 5630105577977904672ull,
    2053435297058519104ull, 13957893008916471936ull, 13340215882559750400ull,
    1841679027358532096ull};


// Input assumption: x%1024 == 0
//
// Computes the 2-adic exponential of x via the sixth degree Taylor
// approximation of exp_2, i.e. 1+x+x^2/2!+x^3/3!+x^4/4!+x^5/5!+x^6/6!
//
// Equivalently, (((x^2/6+x)*1/5+1)*x^2/2+x)*x^2/6+x^2/2+x+1
uint64_t twoadic_exp_1024(uint64_t x) {
    const uint64_t x2over2 = x*(x>>1);
    const uint64_t x2over6 = x2over2*12297829382473034411ull; // 1/3

    return
        (((x2over6+x)*14757395258967641293ull+1)*(x2over2>>1)+x)*x2over6
            + x2over2 + x + 1;
}

uint64_t twoadic_exp_terms(uint64_t x, uint64_t t) {
    const uint64_t x2over2 = x*(x>>1);
    if (t == 2) {
        return 1 + x + x2over2;
    }
    const uint64_t x2over6 = x2over2*12297829382473034411ull;
    if (t == 3) {
        return 1 + x + x2over2 + x2over6*x;
    }
    return ((x2over2>>1)+x)*x2over6+x2over2+x+1;
}

// Input assumption: x%1024 == 1
//
// Computes the 2-adic logarithm of x via the sixth degree Taylor
// approximation of log_2, i.e. x-x^2/2+x^3/3-x^4/4+x^5/5-x^6/6
//
// Equivalently, (((((-1/3*(x/2)+1/5)*x^2)-x/4+1/3)*x^2 - x/2 +1)*x
uint64_t twoadic_log_1024(uint64_t x) {
    x--;
    const uint64_t x2 = x*x;

    return
        (((((-(12297829382473034411ull))*(x>>1) + 14757395258967641293ull)*x2)
           - (x>>2) + 12297829382473034411ull)*x2 - (x>>1) + 1)*x;
}

// Input assumption: x%4 == 1
//
// Computes the 2-adic logarithm of x
uint64_t twoadic_log(uint64_t x) {
    uint64_t err = 0;

    // Reduce to the case x%1024 == 1
    for(int i = 2; i < 10; ++i) {
        if (x & (1ull<<i)) {
            x += (x<<i);
            err += logarray[i]; // log((2^i +1))
        }
    }
    return twoadic_log_1024(x)-err;
}

// Input assumption: x%4 == 0
//
// Computes the 2-adic exponential of x
uint64_t twoadic_exp(uint64_t x) {
    uint64_t rv = 1;

    // Reduce to the case x%1024 == 0
    for(int i = 2; i < 10; ++i) {
        if (x & (1 << i)) {
            // Multiply rv by 2^i + 1, subtract off the log of 2^i + 1.
            rv += (rv << i);
            x -= logarray[i]; // log((2^i +1))
        }
    }

    return twoadic_exp_1024(x)*rv;
}

// Input assumption: x%4 == 0
//
// Computes the 2-adic exponential of x to m digits
uint64_t twoadic_exp_precision(uint64_t x, int m) {
    uint64_t rv = 1;
    // The Taylor series computes T terms. That means the T+1th term must have valuation greater than or equal to m.
    // Which means T+1 * (l-1)+1 >= m. Full precision is T = 6, l = 10 (works because 7th term of taylor series actually 
    // has valuation 70-4 = 66 >= 64).
    //
    // For T =4, l = 7, the 5th term of the taylor series has val >= 32.
    // For T = 3, l = 5 , the 4rd term of the taylor series has val >= 16
    // For T = 2, l = 3, the 3rd term of the taylor series has val >= 8, 4th term has val >= 8.
    // These choices are heuristic and could be optimized. 
    int l, t;
    if (m >= 33) {
        l = 10;
    } else if (m >= 16) {
        l = 7;
        t = 4;
    } else if (m >= 8) {
        l = 5;
        t = 3;
    } else {
        l = 3;
        t = 2;
    }


    // Reduce to the case x%1024 == 0
    for(int i = 2; i < l; ++i) {
        if (x & (1 << i)) {
            // Multiply rv by 2^i + 1, subtract off the log of 2^i + 1.
            rv += (rv << i);
            x -= logarray[i]; // log((2^i +1))
        }
    }
    if (l == 10) {
        return twoadic_exp_1024(x)*rv;
    }
    return twoadic_exp_terms(x, t)*rv;
}



// Input assumption: a%2 == 1
uint64_t adic_pow_odd(uint64_t a, uint64_t b) {
    if (a&2) {
        // Handle 3(mod 4) by negating
        a = -a;
        if (b % 2 != 0) {
            return -twoadic_exp(b*twoadic_log(a));
        }
    }
    return twoadic_exp(b*twoadic_log(a));
}

// Public:
uint64_t adic_pow(uint64_t a, uint64_t b) {
    // Special cases to avoid messiness with ctz builtin
    // Note the convention 0^0 == 1
    if (b == 0) {
        return 1;
    } else if (a == 0) {
        return 0;
    }
    // Handle powers of two directly
    uint64_t natz = __builtin_ctzll(a);
    auto nrtz = natz*b;
    return (nrtz >= 64) ? 0 : ((adic_pow_odd(a>>natz, b))<<nrtz);
}

uint64_t trad_pow(uint64_t a, uint64_t b) {
    uint64_t rv = 1;
    if(b&1) {
        rv*=a;
    }
    b >>= 1;
    while(b != 0) {
        a*=a;
        if(b&1) {
            rv*=a;
        }
        b >>= 1;
    }
    return rv;
}


#endif // LOOKUPAPPROACH_H_INCLUDED
