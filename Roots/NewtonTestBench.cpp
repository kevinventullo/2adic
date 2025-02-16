#include <chrono>
#include <iostream>
#include <random>
#include "Newton.h"



bool CorrectnessTest(int r) {
    std::vector<std::pair<uint64_t, uint64_t>> root_pairs;
    
    for(uint64_t i = 1; i < std::min((1ull << 64/r), 10000000ull); i++){

        uint64_t pow = trad_pow(i, r);
        root_pairs.emplace_back(pow, i);
    }
    for(uint64_t i = (1ull << 64/r)-1; i > (1ull << 64/r)-1000000; i--){

        uint64_t pow = trad_pow(i, r);
        root_pairs.emplace_back(pow, i);
    }

    for(auto p : root_pairs) {
        auto kr = twoAdicKthRoot(p.first, r);
        if (kr != p.second) {
            std::cout << "Adic computed " << r << "th root of " << p.first << " as " << kr << " not " << p.second << std::endl;
            return false;
        }
        kr = NewtonKthRoot(p.first, r);
        if (kr != p.second) {
            std::cout << "Newton computed " << r << "th root of " << p.first << " as " << kr << " not " << p.second << std::endl;
            return false;
        }
        
    }

    return true;
    

}


void Benchmark(int k) {
    std::vector<std::pair<uint64_t, uint64_t>> root_pairs;
    


    for(uint64_t i = 1; i < 1000000; i+=2){
        uint64_t pow = trad_pow(i%(1ull<<64/k), k);
        root_pairs.emplace_back(pow, i%(1ull<<64/k));
    }

    std::vector<uint64_t> results;
    results.reserve(2*root_pairs.size());

    auto s = std::chrono::high_resolution_clock::now();

    for(auto p : root_pairs) {
        results.push_back(twoAdicKthRoot(p.first, k));
    }
    auto m = std::chrono::high_resolution_clock::now();

    for(auto p : root_pairs) {
        results.push_back(NewtonKthRoot(p.first, k));
    }
    auto m2 = std::chrono::high_resolution_clock::now();

    auto first_time = m - s;
    auto second_time(m2 - m);

    std::cout << "Adic: " << std::chrono::duration_cast<std::chrono::milliseconds>(first_time).count() << "ms" << std::endl;
    std::cout << "Newton: " << std::chrono::duration_cast<std::chrono::milliseconds>(second_time).count() << "ms" << std::endl;
    std::cout << std::endl;

    int num_wrong = 0;
    for(unsigned i = 0; i < results.size()/2; ++i) {
        if (results[i] != results[i+results.size()/2]) {
            num_wrong++;
        }
    }

    if (num_wrong!=0) {
        std::cout << "Num wrong = " << num_wrong << std::endl;
    }
}


int main() {
    for (int k = 2; k < 13; ++k) {
        if (!CorrectnessTest(k)) return -1;
    }
    for (int k = 2; k < 13; ++k) {
        std::cout << "Benchmark k = " << k << ":" << std::endl;
        Benchmark(k);
    }
    return 0;
}