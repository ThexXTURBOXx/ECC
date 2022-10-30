#ifndef ECC_UTILS_H
#define ECC_UTILS_H

#include "CoCoA/library.H"

using namespace std;

namespace CoCoA {

    template<class T>
    RingElem getOr(vector <T> vec, size_t i, T defaultVal);

    int parseNum(char c);

    char toChar(long i);

    RingElem toPolynomial(const string &str, long k, ConstRefRingElem x);

    string toString(ConstRefRingElem p, long n, ConstRefRingElem x);

}

#endif //ECC_UTILS_H
