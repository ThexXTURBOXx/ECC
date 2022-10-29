#include "CoCoA/library.H"

using namespace std;

namespace CoCoA {

    int parseNum(const char c) {
        if (c >= '0' && c <= '9')
            return c - '0';
        if (c >= 'A' && c <= 'Z')
            return c - 'A' + 10;
        if (c >= 'a' && c <= 'z')
            return c - 'a' + 10;
        throw invalid_argument("Invalid character!");
    }

    char toChar(const long i) {
        if (i >= 0 && i <= 9)
            return (char) ('0' + i);
        if (i >= 10 && i <= 35)
            return (char) ('A' + (i - 10));
        throw invalid_argument("Invalid number!");
    }

    RingElem toPolynomial(const string &str, const long k, ConstRefRingElem x) {
        RingElem poly = zero(owner(x));
        for (long i = 0; i < k; ++i) {
            poly += parseNum(str[i]) * power(x, k - i - 1);
        }
        return poly;
    }

    template<class T>
    RingElem getOr(vector <T> vec, size_t i, T defaultVal) {
        if (i >= vec.size()) return defaultVal;
        return vec[i];
    }

    string toString(ConstRefRingElem p, const long n, ConstRefRingElem x) {
        const vector<RingElem> coeffVec = CoeffVecWRT(p, x);
        const RingElem z = zero(owner(p));
        string str;
        long buf;
        for (long i = n - 1; i >= 0; --i) {
            if (!IsConvertible(buf, getOr(coeffVec, i, z))) {
                throw invalid_argument("Invalid coefficient!");
            }
            str += toChar(buf);
        }
        return str;
    }

}