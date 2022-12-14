#include <numeric>
#include "CoCoA/library.H"
#include "ecc/rm.H"
#include "util/utils.H"
#include "types/linear.H"

using namespace std;

namespace CoCoA {

  namespace { /* anonymous */
    vector<RingElem> constructVector(ConstRefRingElem zeroR, ConstRefRingElem oneR, const long size,
                                     const long m, const long i) {
      vector<RingElem> ret(size, zeroR);
      const long stepJ = SmallPower(2, m - i);
      const long stepK = SmallPower(2, m - i - 1);
      for (long j = 0; j < size; j += stepJ) {
        for (long k = 0; k < stepK; ++k) {
          ret[j + k] = oneR;
        }
      }
      return ret;
    }

    vector<vector<RingElem>>
    generateAllRows(const ring &R, const long m, const vector<long> &S) { // NOLINT(misc-no-recursion)
      const long size = SmallPower(2, m);
      if (S.empty())
        return vector<vector<RingElem>>{vector<RingElem>(size, one(R))};

      const RingElem &zeroR = zero(R);
      const RingElem &oneR = one(R);
      vector<vector<RingElem>> Srest = generateAllRows(R, m, vector<long>(S.begin() + 1, S.end()));
      vector<RingElem> xiRow = constructVector(zeroR, oneR, size, m, S[0]);
      vector<RingElem> notXiRow;
      notXiRow.reserve(size);
      std::transform(xiRow.begin(), xiRow.end(), std::back_inserter(notXiRow), [oneR](const auto &e) {
        return oneR - e;
      });

      size_t SrestSize = Srest.size();
      vector<vector<RingElem>> ret(2 * SrestSize);

      for (size_t i = 0; i < SrestSize; ++i) {
        vector<RingElem> row = Srest[i];
        vector<RingElem> temp(size);
        for (long j = 0; j < size; ++j) {
          temp[j] = xiRow[j] * row[j];
        }
        ret[i] = temp;
      }
      for (size_t i = 0; i < SrestSize; ++i) {
        vector<RingElem> row = Srest[i];
        vector<RingElem> temp(size);
        for (long j = 0; j < size; ++j) {
          temp[j] = notXiRow[j] * row[j];
        }
        ret[SrestSize + i] = temp;
      }

      return ret;
    }
  }

  long RM::genK(const long r, const long m) {
    BigInt dim = BigInt(0);
    for (long i = 0; i <= r; ++i)
      dim += binomial(m, i);

    long ret;
    if (!IsConvertible(ret, dim))
      CoCoA_THROW_ERROR(ERR::ArgTooBig, "RM ctor");
    return ret;
  }

  vector<vector<RingElem>> RM::genXrows(const ring &R, long m) {
    vector<vector<RingElem>> ret(m);
    const RingElem &zeroR = zero(R);
    const RingElem &oneR = one(R);
    const long size = SmallPower(2, m);
    for (long i = 0; i < m; ++i)
      ret[i] = constructVector(zeroR, oneR, size, m, i);
    return ret;
  }

  matrix RM::genG(const ring &R, const long r, const long m, const vector<vector<RingElem>> &xrows) {
    vector<long> elems(m);
    iota(begin(elems), end(elems), 0);

    long cols = SmallPower(2, m);
    vector<vector<long>> Ss;
    ConstMatrixView G = ZeroMat(R, 0, cols);
    for (long s = 0; s <= r; ++s) {
      Ss.clear();
      subsets(elems, s, Ss);

      for (const auto &S: Ss) {
        vector<RingElem> ret(cols, one(R));
        for (const auto &i: S) {
          vector<RingElem> row = xrows[i];
          for (long j = 0; j < cols; ++j)
            ret[j] *= row[j];
        }
        G = ConcatVer(G, NewDenseMat(RowMat(ret)));
      }
    }
    return NewDenseMat(G);
  }

  vector<vector<vector<RingElem>>> RM::genVotingRows(const ring &R, const long r, const long m) {
    vector<long> elems(m);
    iota(begin(elems), end(elems), 0);

    vector<vector<long>> Ss;
    vector<vector<vector<RingElem>>> ret;
    for (long s = 0; s <= r; ++s) {
      Ss.clear();
      subsets(elems, s, Ss);

      for (const auto &S: Ss) {
        vector<long> is;
        for (long i = 0; i < m; ++i)
          if (find(S.cbegin(), S.cend(), i)==S.cend())
            is.push_back(i); // All the i not contained in S
        ret.push_back(generateAllRows(R, m, is));
      }
    }
    return ret;
  }

  vector<long> RM::genRibd(long r, long m) {
    vector<long> ribd(r + 1, 0);
    for (long i = 1; i <= r; ++i)
      ribd[i] = ribd[i - 1] + binom(m, i);
    return ribd;
  }

  matrix encodeRM(const RM &rm, const matrix &w) {
    return linEncode(rm.G, w);
  }

  matrix decodeRM(const RM &rm, const matrix &w) {
    return w;
  }

}
