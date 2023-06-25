#ifndef ECC_RM_H
#define ECC_RM_H

#include "CoCoA/library.H"

using namespace std;

namespace CoCoA {
  namespace ECC {

    /**
     * Holds the parameters of a Reed-Muller code.
     * @see Muller, D. E. (1954). Application of Boolean algebra to switching circuit design and to error detection
     * @see Reed, I. S. (1953). A class of multiple-error-correcting codes and the decoding scheme
     * @see Raaphorst, S. (2003). Reed-muller codes
     */
    struct RM {
      /**
       * The ring over which the code is defined.
       */
      const ring R;

      /**
       * The degree of the code.
       */
      const long r;

      /**
       * The variety of the code.
       */
      const long m;

      /**
       * The length of the code.
       */
      const long n;

      /**
       * The dimension of the code.
       */
      const long k;

      /**
       * The minimum Hamming distance of the code.
       */
      const long d;

      /**
       * The rows which correspond each to the variable `x_i` of the code.
       */
      const vector<vector<RingElem>> xrows;

      /**
       * The generator matrix of the code.
       */
      const matrix G;

      /**
       * The list of all monomials not in each row.
       */
      const vector<vector<vector<RingElem>>> votingRows;

      /**
       * The list of indices of the rows corresponding to monomials of degree `i`.
       */
      const vector<long> ribd;

      /**
       * Constructs a {@link RM} code.
       * @param r The degree of the code
       * @param m The variety of the code
       */
      RM(const long r, const long m) : RM(NewZZmod(2), r, m) {
      }

      /**
       * Constructs a {@link RM} code.
       * @param R The ring over which the code is defined
       * @param r The degree of the code
       * @param m The variety of the code
       */
      RM(const ring &R, const long r, const long m) : R(R), r(r), m(m), n(SmallPower(2, m)), k(genK(r, m)),
                                                      d(SmallPower(2, m - r)), xrows(genXrows(R, m)),
                                                      G(genG(R, r, m, xrows)), votingRows(genVotingRows(R, r, m)),
                                                      ribd(genRibd(r, m)) {
      }

    private:
      /**
       * Calculates the dimension of the code of given parameters.
       * @param r The degree of the code
       * @param m The variety of the code
       * @return The dimension of the code
       */
      static long genK(long r, long m);

      /**
       * Calculates the rows which correspond each to the variable `x_i` of the code.
       * @param R The ring over which the code is defined
       * @param m The variety of the code
       * @return The rows which correspond each to the variable `x_i` of the code
       */
      static vector<vector<RingElem>> genXrows(const ring &R, long m);

      /**
       * Calculates a generator matrix of the code.
       * @param R The ring over which the code is defined
       * @param r The degree of the code
       * @param m The variety of the code
       * @param xrows The rows which correspond each to the variable `x_i` of the code
       * @return A generator matrix of the code
       */
      static matrix genG(const ring &R, long r, long m, const vector<vector<RingElem>> &xrows);

      /**
       * Calculates the list of all monomials not in each row.
       * @param R The ring over which the code is defined
       * @param r The degree of the code
       * @param m The variety of the code
       * @return The list of all monomials not in each row
       */
      static vector<vector<vector<RingElem>>> genVotingRows(const ring &R, long r, long m);

      /**
       * Calculates the list of indices of the rows corresponding to monomials of degree `i`.
       * @param r The degree of the code
       * @param m The variety of the code
       * @return The list of indices of the rows corresponding to monomials of degree `i`
       */
      static vector<long> genRibd(long r, long m);
    };

    /**
     * Encodes the given word using the given {@link RM} code.
     * @param rm The {@link RM} code to use
     * @param w The word to encode
     * @return The encoded word
     */
    matrix encodeRM(const RM &rm, const matrix &w);

    /**
     * Decodes the given word using the given {@link RM} code using majority voting.
     * @param rm The {@link RM} code to use
     * @param w The word to decode
     * @return The decoded word
     * @see Raaphorst, S. (2003). Reed-muller codes
     */
    matrix decodeRM(const RM &rm, matrix w);

  }
}

#endif // ECC_RM_H
