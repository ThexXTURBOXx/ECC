#include <numeric>
#include "CoCoA/library.H"
#include "ecc/rm.H"
#include "util/utils.H"
#include "types/linear.H"

using namespace std;

namespace CoCoA {
  namespace ECC {
    namespace { /* anonymous */
      /**
       * Calculates the row which corresponds to the variable `x_i` of the code.
       * @param zeroR The zero element of the ring
       * @param oneR The one element of the ring
       * @param size The size of the row
       * @param m The variety of the code
       * @param i The index of the variable
       * @return The row which corresponds to the variable `x_i` of the code
       */
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

      /**
       * Recursively calculates all rows over the monomials in the given set.
       * @param R The ring over which the code is defined
       * @param m The variety of the code
       * @param S The set of monomials
       * @param off The offset of the current monomial
       * @return The list of all rows over the monomials in the given set
       */
      vector<vector<RingElem>> generateAllRows(const ring& R, const long m, const vector<long>& S, const long off) {
        // NOLINT(misc-no-recursion)
        const long size = SmallPower(2, m);
        if (off == S.size())
          return {vector<RingElem>(size, one(R))};

        vector<vector<RingElem>> Srest = generateAllRows(R, m, S, off + 1);
        vector<RingElem> xiRow = constructVector(zero(R), one(R), size, m, S[off]);

        size_t SrestSize = Srest.size();
        vector<vector<RingElem>> ret(2 * SrestSize);

        for (size_t i = 0; i < SrestSize; ++i) {
          vector<RingElem> temp(size);
          vector<RingElem> temp2(size);
          for (long j = 0; j < size; ++j) {
            temp[j] = xiRow[j] * Srest[i][j];
            temp2[j] = (xiRow[j] + 1) * Srest[i][j];
          }
          ret[i] = temp;
          ret[SrestSize + i] = temp2;
        }

        return ret;
      }

      /**
       * Calculates all rows over the monomials in the given set.
       * @param R The ring over which the code is defined
       * @param m The variety of the code
       * @param S The set of monomials
       * @return The list of all rows over the monomials in the given set
       */
      vector<vector<RingElem>> generateAllRows(const ring& R, const long m, const vector<long>& S) {
        return generateAllRows(R, m, S, 0);
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

    vector<vector<RingElem>> RM::genXrows(const ring& R, const long m) {
      vector<vector<RingElem>> ret(m);
      const RingElem& zeroR = zero(R);
      const RingElem& oneR = one(R);
      const long size = SmallPower(2, m);
      for (long i = 0; i < m; ++i)
        ret[i] = constructVector(zeroR, oneR, size, m, i);
      return ret;
    }

    matrix RM::genG(const ring& R, const long r, const long m, const vector<vector<RingElem>>& xrows) {
      vector<long> elems(m);
      iota(begin(elems), end(elems), 0);

      long cols = SmallPower(2, m);
      vector<vector<long>> Ss;
      ConstMatrixView G = ZeroMat(R, 0, cols);
      for (long s = 0; s <= r; ++s) {
        Ss.clear();
        subsets(elems, s, Ss);

        for (const auto& S : Ss) {
          vector<RingElem> ret(cols, one(R));
          for (const auto& i : S)
            for (long j = 0; j < cols; ++j)
              ret[j] *= xrows[i][j];
          G = ConcatVer(G, NewDenseMat(RowMat(ret)));
        }
      }
      return NewDenseMat(G);
    }

    vector<vector<vector<RingElem>>> RM::genVotingRows(const ring& R, const long r, const long m) {
      vector<long> elems(m);
      iota(begin(elems), end(elems), 0);

      vector<vector<long>> Ss;
      vector<vector<vector<RingElem>>> ret;
      for (long s = 0; s <= r; ++s) {
        Ss.clear();
        subsets(elems, s, Ss);

        for (const auto& S : Ss) {
          vector<long> is;
          for (long i = 0; i < m; ++i)
            if (find(S.cbegin(), S.cend(), i) == S.cend())
              is.push_back(i); // All the i not contained in S
          ret.push_back(generateAllRows(R, m, is));
        }
      }
      return ret;
    }

    vector<long> RM::genRibd(const long r, const long m) {
      vector<long> ribd(r + 1, 0);
      for (long i = 1; i <= r; ++i)
        ribd[i] = ribd[i - 1] + binom(m, i);
      return ribd;
    }

    matrix encodeRM(const RM& rm, const matrix& w) {
      return linEncode(rm.G, w);
    }

    matrix decodeRM(const RM& rm, matrix w) {
      matrix word = NewDenseMat(rm.R, 1, rm.k);
      RingElem dotProductTemp = zero(rm.R);
      for (long degree = rm.r; degree >= 0; --degree) {
        long upperR = rm.ribd[degree];
        long lowerR = degree == 0 ? 0 : rm.ribd[degree - 1] + 1;

        for (long pos = lowerR; pos <= upperR; ++pos) {
          long ones = 0;
          long zeros = 0;

          for (const auto& vrow : rm.votingRows[pos]) {
            dotProductTemp = zero(rm.R);
            for (long i = 0; i < rm.n; ++i) {
              dotProductTemp += w(0, i) * vrow[i];
            }
            if (IsZero(dotProductTemp))
              ++zeros;
            else
              ++ones;
          }

          if (ones == zeros)
            CoCoA_THROW_ERROR("Cannot decode!", __func__);

          SetEntry(word, 0, pos, zeros > ones ? 0 : 1);
        }

        for (long i = 0; i < rm.n; ++i) {
          dotProductTemp = zero(rm.R);
          for (long j = lowerR; j <= upperR; ++j) {
            dotProductTemp += word(0, j) * rm.G(j, i);
          }
          SetEntry(w, 0, i, w(0, i) + dotProductTemp);
        }
      }
      return word;
    }
  }
}
