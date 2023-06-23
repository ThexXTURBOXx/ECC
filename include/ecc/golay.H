#ifndef ECC_GOLAY_H
#define ECC_GOLAY_H

#include "CoCoA/library.H"
#include "types/linear.H"

using namespace std;

namespace CoCoA {

  /**
   * Returns the part of the generator matrix of the Golay code of order `n` that is not the identity matrix,
   * defined over the given ring.
   * @param R The ring over which the matrix should be defined.
   * @param n The order of the matrix.
   * @return The "non-identity-matrix-part" of the generator matrix of the Golay code of order `n`.
   */
  matrix golayMatPart(const ring &R, long n);

  /**
   * Returns the part of the generator matrix of the Golay code of order `n` that is not the identity matrix,
   * constructing a new appropriate ring.
   * @param n The order of the matrix.
   * @return The "non-identity-matrix-part" of the generator matrix of the Golay code of order `n`.
   */
  matrix golayMatPart(long n);

  /**
   * Returns the generator matrix of the Golay code of order `n`, defined over the given ring.
   * @param R The ring over which the matrix should be defined.
   * @param n The order of the matrix.
   * @return The generator code of the Golay code of order `n`.
   */
  matrix golayMat(const ring &R, long n);

  /**
   * Returns the generator matrix of the Golay code of order `n`, constructing a new appropriate ring.
   * @param n The order of the matrix.
   * @return The generator code of the Golay code of order `n`.
   */
  matrix golayMat(long n);

  /**
   * Holds the parameters of a Golay code.
   * @see Golay, M. J. (1949). Notes on digital coding
   */
  struct Golay {
    /**
     * The part of the generator matrix of the code that is not the identity matrix.
     */
    const matrix A;

    /**
     * The ring over which the code is defined.
     */
    const ring R;

    /**
     * The generator matrix of the code.
     */
    const matrix G;

    /**
     * The part of the generator matrix of the extended version of the code that is not the identity matrix.
     */
    const matrix AExt;

    /**
     * The generator matrix of the extended version of the code.
     */
    const matrix GExt;

    /**
     * The order of the field.
     */
    const long q;

    /**
     * The length of the code.
     */
    const long n;

    /**
     * The dimension of the code.
     */
    const long k;

    /**
     * The minimal Hamming distance of the code.
     */
    const long d;

    /**
     * Constructs a {@link Golay} code.
     * @param n The order of the code
     */
    explicit Golay(const long n)
        : A(golayMatPart(n)), R(RingOf(A)), G(NewDenseMat(ConcatHor(IdentityMat(R, NumRows(A)), A))),
          AExt(IsEven(n) ? A:golayMatPart(R, n + 1)),
          GExt(IsEven(n) ? G:NewDenseMat(ConcatHor(IdentityMat(R, NumRows(AExt)), AExt))),
          q(n > 20 ? 2:3), n(n), k((n + 1) / 2), d(n / (5 - q)) /* magic formula */
    {
    }

    /**
     * Constructs a {@link Golay} code using the given ring.
     * @param R The ring over which the code should be defined
     * @param n The order of the code
     */
    explicit Golay(const ring& R, const long n)
        : A(golayMatPart(R, n)), R(R), G(NewDenseMat(ConcatHor(IdentityMat(R, NumRows(A)), A))),
          AExt(IsEven(n) ? A:golayMatPart(R, n + 1)),
          GExt(IsEven(n) ? G:NewDenseMat(ConcatHor(IdentityMat(R, NumRows(AExt)), AExt))),
          q(n > 20 ? 2:3), n(n), k((n + 1) / 2), d(n / (5 - q)) /* magic formula */
    {
    }
  };

  /**
   * Encodes the given word using the given {@link Golay} code.
   * @param golay The {@link Golay} code to use
   * @param w The word to encode
   * @return The encoded word
   */
  matrix encodeGolay(const Golay &golay, const matrix &w);

  /**
   * Decodes the given word using the given {@link Golay} code.
   * @param golay The {@link Golay} code to use
   * @param w The word to decode
   * @return The decoded word
   */
  matrix decodeGolay(const Golay &golay, const matrix &w);

}

#endif // ECC_GOLAY_H
