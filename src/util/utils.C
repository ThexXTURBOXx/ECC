#include "utils.H"
#include "CoCoA/library.H"

using namespace std;

namespace CoCoA {

    template<class T>
    RingElem getOr(vector<T> vec, size_t i, T defaultVal) {
        if (i >= vec.size()) return defaultVal;
        return vec[i];
    }

    int parseNum(const char c) {
        if (c >= '0' && c <= '9')
            return c - '0';
        if (c >= 'A' && c <= 'Z')
            return c - 'A' + 10;
        if (c >= 'a' && c <= 'z')
            return c - 'a' + 10;
        CoCoA_THROW_ERROR("Invalid character", __func__);
        return 0; // Shut up compiler warnings
    }

    char toChar(const long i) {
        if (i >= 0 && i <= 9)
            return (char) ('0' + i);
        if (i >= 10 && i <= 35)
            return (char) ('A' + (i - 10));
        CoCoA_THROW_ERROR("Invalid number", __func__);
        return ' '; // Shut up compiler warnings
    }

    RingElem toPolynomial(const string &str, const long k, ConstRefRingElem x) {
        RingElem poly = zero(owner(x));
        for (long i = 0; i < k; ++i) {
            poly += parseNum(str[i]) * power(x, k - i - 1);
        }
        return poly;
    }

    string toString(ConstRefRingElem p, const long n, ConstRefRingElem x) {
        const vector<RingElem> coeffVec = CoeffVecWRT(p, x);
        const RingElem z = zero(owner(p));
        string str;
        long buf;
        for (long i = n - 1; i >= 0; --i) {
            if (!IsConvertible(buf, getOr(coeffVec, i, z))) {
                CoCoA_THROW_ERROR("Invalid coefficient!", __func__);
            }
            str += toChar(buf);
        }
        return str;
    }

}