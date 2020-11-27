#ifndef NAIVEAPPROACH_H_INCLUDED
#define NAIVEAPPROACH_H_INCLUDED

// Assume x == 1(mod 4)
// Need to compute Taylor series up to x^32/32
uint64_t two_log(uint64_t x) {
    x--; // Taylor expansion is in terms of x-1


    uint64_t rv = x;
    uint64_t curr = x; // Keeps track of the power of x we're working with,

    // -x^2/2
    curr *= (x>>1);
    rv -= curr;
    curr <<= 1;

    // +x^3/3
    curr*=x;
    rv += curr*12297829382473034411ull; // 3^{-1}

    // -x^4/4
    curr*= (x>>2);
    rv -= curr;
    curr <<= 2;

    curr*=x;
    rv += curr*14757395258967641293ull; // 5^{-1}

    curr*= (x>>1);
    rv -= curr*12297829382473034411ull; // 3^{-1}
    curr <<= 1;

    curr*=x>>1;
    rv += curr*2*7905747460161236407ull; // 7^{-1}

    curr*= (x>>2);
    rv -= curr;
    curr <<= 3;

    curr*=x;
    rv += curr*10248191152060862009ull; // 9^{-1}

    curr*= (x>>1);
    rv -= curr*14757395258967641293ull; // 5^{-1}
    curr <<= 1;

    curr*=x;
    rv += curr*3353953467947191203ull; // 11^{-1}

    curr*= (x>>2);
    rv -= curr*12297829382473034411ull; // 3^{-1}
    curr <<= 2;

    curr*=x;
    rv += curr*5675921253449092805ull; // 13^{-1}

    curr*= (x>>1);
    rv -= curr*7905747460161236407ull; // 7^{-1}
    curr <<= 1;

    curr*=(x>>2);
    rv += curr*4*17216961135462248175ull; // 15^{-1}

    curr*= (x>>2);
    rv -= curr;
    curr <<= 4;

    curr*=x;
    rv += curr*17361641481138401521ull; // 17^{-1}

    curr*= (x>>1);
    rv -= curr*10248191152060862009ull; // 9^{-1}
    curr <<= 1;

    curr*=x;
    rv += curr*9708812670373448219ull; // 19^{-1}

    curr *= (x>>2);
    rv -= curr*14757395258967641293ull; // 5^{-1}
    curr <<= 2;

    curr *= x;
    rv += curr*14933078535860113213ull; // 21^{-1}

    curr *= (x>>1);
    rv -= curr*3353953467947191203ull; // 11^{-1}
    curr <<= 1;

    curr *= (x>>1);
    rv += curr*2*15238614669586151335ull; // 23^{-1}

    curr *= (x>>2);
    rv -= curr*12297829382473034411ull; // 3^{-1}
    curr <<= 3;

    curr *= x;
    rv += curr*10330176681277348905ull; // 25^{-1}

    curr *= (x>>1);
    rv -= curr*5675921253449092805ull; // 13^{-1}
    curr <<= 1;

    curr *= x;
    rv += curr*9564978408590137875ull; // 27^{-1}

    curr *= (x>>2);
    rv -= curr*7905747460161236407ull; // 7^{-1}
    curr <<= 2;

    curr *= x;
    rv += curr*3816567739388183093ull; // 29^{-1}

    curr *= (x>>1);
    rv -= curr*17216961135462248175ull; // 15^{-1}
    // Curr = x^30/2

    curr *= (x>>2);
    rv += curr*8*17256631552825064415ull; // 31^{-1}
    // Curr = x^31/8

    curr *= (x >> 2);
    rv -= curr;

    return rv;
}

// Assume x = 0 (mod 4)
// Need to compute Taylor series up to x^58/58!
uint64_t two_exp(uint64_t x) {

    uint64_t rv = 1+x;
    uint64_t curr = x*(x >> 1);

    rv+=curr;

    curr*=12297829382473034411ull; // 3^{-1}
    curr*=x;
    rv+=curr;

    curr >>= 2; // 4^{-1}
    curr*=x;
    rv+=curr;

    curr*=14757395258967641293ull; // 5^{-1}
    curr*=x;
    rv+=curr;

    curr >>= 1; // 2^{-1}
    curr*=12297829382473034411ull; // 3^{-1}
    curr*=x;
    rv+=curr;

    curr*=7905747460161236407ull; // 7^{-1}
    curr*=(x>>1);
    rv+=(curr<<1);
    // curr == x^7/7!*2

    curr*=(x>>2);
    rv+=curr;
    // curr == x^8/8!

    curr*=10248191152060862009ull; // 9^{-1}
    curr*=x;
    rv+=curr;

    curr >>= 1; // 2^{-1}
    curr*=14757395258967641293ull; // 5^{-1}
    curr*=x;
    rv+=curr;

    curr*=3353953467947191203ull; // 11^{-1}
    curr*=x;
    rv+=curr;

    curr >>= 2; // 4^{-1}
    curr*=12297829382473034411ull; // 3^{-1}
    curr*=x;
    rv+=curr;

    curr*=5675921253449092805ull; // 13^{-1}
    curr*=x;
    rv+=curr;

    curr >>= 1;
    curr*=7905747460161236407ull; // 7^{-1}
    curr*=x;
    rv+=curr;

    curr*=17216961135462248175ull; // 15^{-1}
    curr*=(x>>2);
    rv+=(curr << 2);
    // Curr = x^15/15!*4

    curr*=(x>>2);
    rv+=curr;
    // Curr = x^16/16!

    curr*=17361641481138401521ull; // 17^{-1}
    curr*=x;
    rv+=curr;

    curr >>= 1; // 2^{-1}
    curr*=10248191152060862009ull; // 9^{-1}
    curr*=x;
    rv+=curr;

    curr*=9708812670373448219ull; // 19^{-1}
    curr*=x;
    rv+=curr;

    curr >>= 2; // 4^{-1}
    curr*=14757395258967641293ull; // 5^{-1}
    curr*=x;
    rv+=curr;


    curr*=14933078535860113213ull; // 21^{-1}
    curr*=x;
    rv+=curr;

    curr >>= 1;
    curr*= 3353953467947191203ull; // 11^{-1}
    curr*=x;
    rv+=curr;

    curr*=15238614669586151335ull; // 23^{-1}
    curr*=(x>>1);
    rv+=(curr<<1);
    // curr == x^23/23!*2

    curr *= 12297829382473034411ull; // 3^{-1}
    curr *= (x>>2);
    rv+=curr;

    curr *= 10330176681277348905ull; // 25^{-1}
    curr *= x;
    rv+= curr;

    curr >>= 1;
    curr *= 5675921253449092805ull; // 13^{-1}
    curr *= x;
    rv += curr;

    curr *= 9564978408590137875ull; // 27^{-1}
    curr *= x;
    rv += curr;

    curr >>= 2;
    curr *= 7905747460161236407ull; // 7^{-1}
    curr *= x;
    rv += curr;

    curr *= 3816567739388183093ull; // 29^{-1}
    curr *= x;
    rv += curr;

    curr *= 17216961135462248175ull; // 15^{-1}
    curr *= (x>>2);
    rv += (curr<<1);
    // curr == x^30/2*30!

    curr *= 17256631552825064415ull; // 31^{-1}
    curr *= (x>>2);
    rv += (curr<<3);
    // curr == x^31/8*31!

    curr *= (x>>2);
    rv += curr;

    curr *= 1117984489315730401ull; // 33^{-1}
    curr *= x;
    rv += curr;

    curr >>= 1;
    curr *= 17361641481138401521ull; // 17^{-1}
    curr *= x;
    rv += curr;

    curr *= 12649195936257978251ull; // 35^{-1}
    curr *= x;
    rv += curr;


    curr >>= 2;
    curr *= 10248191152060862009ull; // 9^{-1}
    curr *= x;
    rv += curr;

    curr *= 1495681951922396077ull; // 37^{-1}
    curr *= x;
    rv += curr;

    curr >>= 1;
    curr *= 9708812670373448219ull; // 19^{-1}
    curr *= x;
    rv += curr;

    curr *= 8040888442386214807ull; // 39^{-1}
    curr *= (x>>1);
    rv += (curr<<1);
    // curr == x^39/39!*2

    curr *= 14757395258967641293ull; // 5^{-1}
    curr *= (x>>2);
    rv+=curr;

    curr *= 10348173504763894809ull; // 41^{-1}
    curr *= x;
    rv += curr;

    curr >>= 1;
    curr *= 14933078535860113213ull; // 21^{-1}
    curr *= x;
    rv += curr;

    curr *= 9437869060967677571ull; // 43^{-1}
    curr *= x;
    rv += curr;

    curr >>= 2;
    curr *= 3353953467947191203ull; // 11^{-1}
    curr *= x;
    rv += curr;

    curr *= 5738987045154082725ull; // 45^{-1}
    curr *= x;
    rv += curr;

    curr >>= 1;
    curr *= 15238614669586151335ull; // 23^{-1}
    curr *= x;
    rv += curr;

    curr *= 5887258746928580303ull; // 47^{-1}
    curr *= (x>>2);
    rv += (curr<<2);
    // curr == x^47/4*47!

    curr *= 12297829382473034411ull; // 3^{-1}
    curr *= (x>>2);
    rv += curr;

    curr *= 9035139954469984465ull; // 49^{-1}
    curr *= x;
    rv += curr;

    curr >>= 1;
    curr *= 10330176681277348905ull; // 25^{-1}
    curr *= x;
    rv += curr;

    curr *= 18085043209519168251ull; // 51^{-1}
    curr *= x;
    rv += curr;

    curr >>= 2;
    curr *= 5675921253449092805ull; // 13^{-1}
    curr *= x;
    rv += curr;

    curr *= 2436362424829563421ull; // 53^{-1}
    curr *= x;
    rv += curr;

    curr >>= 1;
    curr *= 9564978408590137875ull; // 27^{-1}
    curr *= x;
    rv += curr;

    curr *= 8049488323073258887ull; // 55^{-1}
    curr *= (x>>1);
    rv += (curr<<1);
    // curr == x^55/2*55!

    curr *= 7905747460161236407ull; // 7^{-1}
    curr *= (x>>2);
    rv += curr;

    curr *= 9385185581360999945ull; // 57^{-1}
    curr *= x;
    rv += curr;

    curr >>= 1;
    curr *= 3816567739388183093ull; // 29^{-1}
    curr *= x;
    rv += curr;

    return rv;
}

#endif // NAIVEAPPROACH_H_INCLUDED
