#include <iostream>
#include "NaiveApproach.h"


int main() {
    uint64_t num_errors = 0;

    for(uint64_t i = 1; i < 4000; i+=4) {
        for (uint64_t j = 1; j < 4000; j+=4) {
            if (i*j != two_exp(two_log(i) + two_log(j))) {
                num_errors++;
            }
        }
    }

    for(uint64_t i = 0; i < 4000; i+=4) {
        for (uint64_t j = 0; j < 4000; j+=4) {
            if (i+j != two_log(two_exp(i)*two_exp(j))) {
                num_errors++;
            }
        }
    }

    std::cout << num_errors << " errors" << std::endl;

    return 0;
}
