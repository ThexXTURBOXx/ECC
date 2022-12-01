#include <numeric>
#include "util/utils.H"
#include "CoCoA/library.H"

using namespace std;

namespace CoCoA {

  template<class T>
  RingElem getOr(const vector<T> &vec, size_t i, T defaultVal) {
    if (i >= vec.size())
      return defaultVal;
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

  RingElem toPolynomial(const string &str, ConstRefRingElem x) {
    const long k = (long) str.size();
    RingElem poly = zero(owner(x));
    for (long i = 0; i < k; ++i) {
      poly += parseNum(str[i]) * power(x, k - i - 1);
    }
    return poly;
  }

  matrix toMatrix(const string &str, const ring &R) {
    const long k = (long) str.size();
    matrix m = NewDenseMat(R, 1, k);
    for (long i = 0; i < k; ++i) {
      SetEntry(m, 0, i, parseNum(str[i]));
    }
    return m;
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

  string toString(const matrix &m) {
    const long n = NumCols(m);
    string str;
    long buf;
    for (long i = 0; i < n; ++i) {
      if (!IsConvertible(buf, m(0, i))) {
        CoCoA_THROW_ERROR("Invalid coefficient!", __func__);
      }
      str += toChar(buf);
    }
    return str;
  }

  template<class T>
  vector<vector<T>> tuples(const vector<T> &set, const long tupleSize) {
    vector<vector<T>> result;

    const long maxValue = SmallPower(set.size(), tupleSize);
    for (unsigned long counter = 0; counter < maxValue; ++counter) {
      vector<T> tuple(tupleSize);

      unsigned long currentValue = counter;
      for (long i = 0; i < tupleSize; ++i) {
        unsigned long digit = currentValue % set.size();
        tuple[tupleSize - i - 1] = set[digit];
        currentValue /= set.size();
      }

      result.push_back(tuple);
    }

    return result;
  }

  // Explicitly generate template function to avoid linker errors
  template vector<vector<long>> tuples<long>(const vector<long> &, long);

  matrix e(const ring &R, const long i, const RingElem &b, const long n) {
    matrix m = NewDenseMat(ZeroMat(R, 1, n));
    SetEntry(m, 0, i, b);
    return m;
  }

  template<class T>
  vector<T> cycShift(const vector<T> &vec, const long s) {
    const long len = (long) vec.size();
    vector<T> ret(len);
    for (long i = 0; i < len; ++i) {
      ret[(i + s + len) % len] = vec[i];
    }
    return ret;
  }

  matrix revCirculantMatrix(const vector<RingElem> &firstRow) {
    const long len = (long) firstRow.size();
    vector<vector<RingElem>> m(len);
    m[0] = firstRow;
    for (long i = 1; i < len; ++i) {
      m[i] = cycShift(firstRow, -i);
    }
    return NewDenseMat(owner(firstRow[0]), m);
  }

  long wt(const vector<long> &v) {
    return accumulate(v.cbegin(), v.cend(), 0L,
                      [](const long a, const long b) {
                        return a + sign(b);
                      });
  }

  long wt(const ConstMatrixView &m) {
    long cols = NumCols(m);
    long rows = NumRows(m);
    long ret = 0;
    for (long i = 0; i < rows; ++i) {
      for (long j = 0; j < cols; ++j) {
        if (!IsZero(m(i, j))) {
          ++ret;
        }
      }
    }
    return ret;
  }

  vector<long> ChienSearch(ConstRefRingElem f, ConstRefRingElem a, const long qn, ConstRefRingElem x) {
    const ring &R = owner(f);
    const RingElem z = zero(R);
    const long n = deg(f);
    vector<long> rootPowers = {};
    vector<RingElem> b = CoeffVecWRT(f, x);
    for (long i = 0; i < qn - 1; ++i) {
      if (IsZero(accumulate(b.cbegin(), b.cend(), z)))
        rootPowers.push_back(i);
      for (long j = 0; j <= n; ++j)
        b[j] *= power(a, j);
    }
    return rootPowers;
  }

  long ChienSearchSingleRoot(ConstRefRingElem f, ConstRefRingElem a, const long qn, ConstRefRingElem x) {
    const ring &R = owner(f);
    const RingElem z = zero(R);
    const long n = deg(f);
    vector<RingElem> b = CoeffVecWRT(f, x);
    for (long i = 0; i < qn - 1; ++i) {
      if (IsZero(accumulate(b.cbegin(), b.cend(), z)))
        return i;
      for (long j = 0; j <= n; ++j)
        b[j] *= power(a, j);
    }
    CoCoA_THROW_ERROR("Polynomial does not have a root", __func__);
    return 0;
  }

  RingElem getUniPoly(const vector<RingElem> &G, const long indetIndex, ConstRefRingElem fallback) {
    for (auto &g: G) {
      if (!IsConstant(g) && UnivariateIndetIndex(g)==indetIndex)
        return g;
    }
    return fallback;
  }

}
