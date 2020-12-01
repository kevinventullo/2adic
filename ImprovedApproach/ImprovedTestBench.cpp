#include <chrono>
#include <iostream>
#include <random>
#include "LookupApproach.h"
#include "CombinedApproach.h"

void CorrectnessTest() {
    uint64_t num_errors = 0;
    const auto tsize = 3000;

    for(uint64_t i = 0; i < tsize; i++) {
        for (uint64_t j = 0; j < tsize; j++) {
            if (trad_pow(i,j) != adic_pow(i,j) || adic_pow(i,j) != combined_pow(i,j)) {
                num_errors++;
            }
        }
    }

    std::cout << num_errors << " errors among base, exp < " << tsize << std::endl;
    std::cout << std::endl;
}

void Benchmark() {
    const uint64_t num_samples = 1000000;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_int_distribution<uint64_t> distribution(0, -1ull);

    std::vector<uint64_t> as, bs;
    for(uint64_t i = 0; i < num_samples; ++i) {
        auto temp = distribution(generator);
        // Force the base to be odd, since even bases almost always exponentiate to zero.
        as.push_back(temp|1);
        temp = distribution(generator);
        bs.push_back(temp);
    }
    std::vector<uint64_t> results;
    results.reserve(4*num_samples);
    auto s = std::chrono::high_resolution_clock::now();

    for(unsigned i = 0; i < num_samples; ++i) {
        results.push_back(trad_pow(as[i],bs[i]));
    }
    auto m = std::chrono::high_resolution_clock::now();

    for(unsigned i = 0; i < num_samples; ++i) {
        results.push_back(adic_pow(as[i], bs[i]));
    }
    auto m2 = std::chrono::high_resolution_clock::now();

    for(unsigned i = 0; i < num_samples; ++i) {
        results.push_back(combined_pow(as[i], bs[i]));
    }
    auto m3 = std::chrono::high_resolution_clock::now();

    for(unsigned i = 0; i < num_samples; ++i) {
        results.push_back(as[i]);
    }

    auto e = std::chrono::high_resolution_clock::now();

    auto first_time = m - s;
    auto second_time(m2 - m);
    auto third_time(m3-m2);
    auto fourth_time(e-m3);

    std::cout << "Standard: " << std::chrono::duration_cast<std::chrono::milliseconds>(first_time).count() << "ms" << std::endl;
    std::cout << "Adic: " << std::chrono::duration_cast<std::chrono::milliseconds>(second_time).count() << "ms" << std::endl;
    std::cout << "Combined: " << std::chrono::duration_cast<std::chrono::milliseconds>(third_time).count() << "ms" << std::endl;
    std::cout << "Baseline: " << std::chrono::duration_cast<std::chrono::milliseconds>(fourth_time).count() << "ms" << std::endl;
    std::cout << std::endl;
    int num_wrong = 0;

    for(unsigned i = 0; i < num_samples; ++i) {
        if (results[i] != results[i+num_samples]) {
            num_wrong++;
        }
    }

    std::cout << num_wrong << " wrong " << std::endl;
}


int main() {
    CorrectnessTest();
    Benchmark();
    return 0;
}
