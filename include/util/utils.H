#ifndef ECC_UTILS_H
#define ECC_UTILS_H

#include "CoCoA/library.H"

using namespace std;

namespace CoCoA {

    template<class T>
    RingElem getOr(const vector <T> &vec, size_t i, T defaultVal);

    int parseNum(char c);

    char toChar(long i);

    RingElem toPolynomial(const string &str, ConstRefRingElem x);

    matrix toMatrix(const string &str, const ring &R);

    string toString(ConstRefRingElem p, long n, ConstRefRingElem x);

    string toString(const matrix &m);

    template<class T>
    vector <vector<T>> tuples(const vector <T> &set, long tupleSize);

}

#endif //ECC_UTILS_H
