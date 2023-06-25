#include <numeric>
#include "util/utils.H"
#include "CoCoA/library.H"

using namespace std;

namespace CoCoA {
  namespace ECC {

    template<class T>
    T getOr(const vector<T> &vec, const size_t i, T defaultVal) {
      if (i >= vec.size())
        return defaultVal;
      return vec[i];
    }

    /**
     * Explicitly instantiate template function to avoid linker errors.
     * @return {@link getOr}
     */
    template RingElem getOr<RingElem>(const vector<RingElem> &, const size_t, RingElem);

    static int parseNum(const char c) {
      if (c >= '0' && c <= '9')
        return c - '0';
      if (c >= 'A' && c <= 'Z')
        return c - 'A' + 10;
      if (c >= 'a' && c <= 'z')
        return c - 'a' + 10;
      CoCoA_THROW_ERROR("Invalid character", __func__);
    }

    char toChar(const long i) {
      if (i >= 0 && i <= 9)
        return (char) ('0' + i);
      if (i >= 10 && i <= 35)
        return (char) ('A' + (i - 10));
      CoCoA_THROW_ERROR("Invalid number", __func__);
    }

    RingElem toPolynomial(const string &str, ConstRefRingElem x) {
      const long k = (long) str.size();
      RingElem poly = zero(owner(x));
      for (long i = 0; i < k; ++i)
        poly += parseNum(str[i]) * power(x, k - i - 1);
      return poly;
    }

    RingElem toPolynomial(const matrix &mat, ConstRefRingElem x) {
      const long k = NumCols(mat);
      RingElem poly = zero(owner(x));
      for (long i = 0; i < k; ++i)
        poly += mat(0, i) * power(x, k - i - 1);
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

    matrix toMatrix(ConstRefRingElem p, long n, ConstRefRingElem x) {
      const vector<RingElem> coeffVec = CoeffVecWRT(p, x);
      const RingElem z = zero(owner(p));
      matrix m = NewDenseMat(owner(p), 1, n);
      long buf;
      for (long i = n - 1; i >= 0; --i) {
        if (!IsConvertible(buf, getOr(coeffVec, i, z)))
          CoCoA_THROW_ERROR("Invalid coefficient!", __func__);
        SetEntry(m, 0, n - (i + 1), buf);
      }
      return m;
    }

    string toString(ConstRefRingElem p, const long n, ConstRefRingElem x) {
      const vector<RingElem> coeffVec = CoeffVecWRT(p, x);
      const RingElem z = zero(owner(p));
      string str;
      long buf;
      for (long i = n - 1; i >= 0; --i) {
        if (!IsConvertible(buf, getOr(coeffVec, i, z)))
          CoCoA_THROW_ERROR("Invalid coefficient!", __func__);
        str += toChar(buf);
      }
      return str;
    }

    string toString(const matrix &m) {
      const long n = NumCols(m);
      string str;
      long buf;
      for (long i = 0; i < n; ++i) {
        if (!IsConvertible(buf, m(0, i)))
          CoCoA_THROW_ERROR("Invalid coefficient!", __func__);
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

    /**
     * Explicitly instantiate template function to avoid linker errors.
     * @return {@link tuples}
     */
    template vector<vector<long>> tuples<long>(const vector<long> &, const long);

    /**
     * Explicitly instantiate template function to avoid linker errors.
     * @return {@link tuples}
     */
    template vector<vector<RingElem>> tuples<RingElem>(const vector<RingElem> &, const long);

    namespace { /* anonymous */
      /**
       * Recursive function to calculate all subsets of a given size of a given set.
       * @tparam T The generic type of the vector
       * @param arr The set of which the subsets should be calculated
       * @param size The size of the set
       * @param left The number of elements that should still be added to the subset
       * @param index The index of the current element
       * @param l The current subset
       * @param bl The buffer for the set of all subsets
       * @see This is an improved version of `examples/ex-MVT-Simplicial.C` in CoCoALib
       */
      template<class T>
      void subsetsInternal(const vector<T> &arr, const int size, const long left, // NOLINT(misc-no-recursion)
                           const int index, vector<T> &l, vector<vector<T>> &bl) {
        if (left == 0) {
          bl.push_back(l);
          return;
        }
        for (int i = index; i < size; i++) {
          l.push_back(arr[i]);
          subsetsInternal(arr, size, left - 1, i + 1, l, bl);
          l.pop_back();
        }
      }
    }

    template<class T>
    void subsets(const vector<T> &set, const long setSize, vector<vector<T>> &ret) {
      vector<T> buf;
      return subsetsInternal(set, set.size(), setSize, 0, buf, ret);
    }

    /**
     * Explicitly instantiate template function to avoid linker errors.
     */
    template void subsets<long>(const vector<long> &, const long, vector<vector<long>> &);

    /**
     * Explicitly instantiate template function to avoid linker errors.
     */
    template void subsets<RingElem>(const vector<RingElem> &, const long, vector<vector<RingElem>> &);

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
    }

    RingElem getUniPoly(const vector<RingElem> &G, const long IndetIndex, ConstRefRingElem fallback) {
      for (auto &g: G) {
        if (!IsConstant(g) && UnivariateIndetIndex(g) == IndetIndex)
          return g;
      }
      return fallback;
    }

    long binom(const long n, const long k) {
      BigInt b = binomial(n, k);
      long ret;
      if (!IsConvertible(ret, b))
        CoCoA_THROW_ERROR(ERR::ArgTooBig, "binom");
      return ret;
    }

    long GetLength(const matrix &m) {
      long c = m->myNumCols();
      long r = m->myNumRows();
      return max(c, r);
    }

    // As of CoCoALib 0.99815, these have been incorporated into CoCoALib itself
    /*
     * Checks if f is a primitive polynomial in ZZ/(p), skipping sanity checks.
     * @param f The polynomial to check
     * @return Whether f is primitive
     */
    /*
    bool IsPrimitivePoly_NoArgChecks(ConstRefRingElem f) {
      if (IsZero(ConstantCoeff(f))) return false;
      if (!IsOne(LC(f))) return false;
      if (!IsIrred(f)) return false;

      const ring &Px = owner(f);
      const long IndetIndex = UnivariateIndetIndex(f);

      // We check if f is an n-th primitive polynomial in ZZ/(p)
      const long n = deg(f);
      const BigInt p = characteristic(CoeffRing(Px));
      const BigInt M = power(p, n) - 1;
      const vector<BigInt> fac = factor(M).myFactors();
      return none_of(fac.begin(), fac.end(), [f, Px, IndetIndex, M](const BigInt &m) {
        return IsOne(NR(IndetPower(Px, IndetIndex, M / m), {f}));
      });
    }

    bool IsPrimitivePoly(ConstRefRingElem f) {
      const char *const FnName = "IsPrimitivePoly";

      if (IsZero(f)) return false;
      const ring &Px = owner(f);
      if (!IsSparsePolyRing(Px))
        CoCoA_THROW_ERROR(ERR::NotSparsePolyRing, FnName);
      if (UnivariateIndetIndex(f) < 0)
        CoCoA_THROW_ERROR(ERR::NotUnivariate, FnName);
      const ring &P = CoeffRing(Px);
      if (!IsFiniteField(P) || !IsQuotientRing(P) || !IsZZ(BaseRing(P))) // TODO: Does that only allow ZZ/(p)?
        CoCoA_THROW_ERROR(ERR::BadRing, FnName);

      return IsPrimitivePoly_NoArgChecks(f);
    }
    */

    RingElem BruteForcePrimPoly(const ring &Px, const long n, const long IndetIndex) {
      const char *const FnName = "BruteForcePrimPoly";

      if (n < 1)
        CoCoA_THROW_ERROR(ERR::BadArg, FnName);
      if (!IsSparsePolyRing(Px))
        CoCoA_THROW_ERROR(ERR::NotSparsePolyRing, FnName);
      const ring &P = CoeffRing(Px);
      if (!IsFiniteField(P) || !IsQuotientRing(P) || !IsZZ(BaseRing(P))) // TODO: Does that only allow ZZ/(p)?
        CoCoA_THROW_ERROR(ERR::BadRing, FnName);

      const BigInt p = characteristic(P);
      if (n == 1) {
        RingElem f = indet(Px, IndetIndex) + 1;
        for (long i = 1; i < p; ++i) {
          if (IsPrimitivePoly(f))
            return f;
          f += 1;
        }
      }

      // Randomised version
      RandomSeqBigInt seqNonZero(1, p - 1);
      RandomSeqBigInt seqDef(0, p - 1);
      RingElem f = IndetPower(Px, IndetIndex, n) + *seqNonZero;
      ++seqNonZero;
      while (true) {
        for (long i = 1; i < n; ++i) {
          f += *seqDef * IndetPower(Px, IndetIndex, i);
          ++seqDef;
        }
        if (IsPrimitivePoly(f))
          return f;
        f = IndetPower(Px, IndetIndex, n) + *seqNonZero;
        ++seqNonZero;
      }
    }

    bool isMultiple(const vector<RingElem> &a, const vector<RingElem> &b, const long c) {
      for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] != c * b[i])
          return false;
      }
      return true;
    }

    long divide(const vector<RingElem> &a, const vector<RingElem> &b, const long q) {
      for (long i = 1; i <= q; ++i) {
        if (isMultiple(a, b, i))
          return i;
      }
      return 0;
    }

  }
}
