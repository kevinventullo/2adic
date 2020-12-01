#ifndef COMBINEDAPPROACH_H_INCLUDED
#define COMBINEDAPPROACH_H_INCLUDED

#include "LookupApproach.h"

uint64_t combined_pow(uint64_t a, uint64_t b) {
    uint64_t rv = 1;
    for(size_t i = 0; i < 8; ++i) {
        if(b&1) {
            rv*=a;
        }
        b >>= 1;
        a*=a;
    }
    if (b == 0) {
        return rv;
    } else if (a == 0) {
        return 0;
    } else {
        return rv*twoadic_exp_1024(b*twoadic_log_1024(a));
    }
}

#endif // COMBINEDAPPROACH_H_INCLUDED
