#ifndef ECC_HAM_H
#define ECC_HAM_H

#include "CoCoA/library.H"
#include "types/linear.H"

using namespace std;

namespace CoCoA {

  /**
   * Generates the parity check matrix of a Hamming code.
   * @param R The ring over which the code should be defined
   * @param r The number of parity bits
   * @param q The order of the field
   * @return The parity check matrix
   */
  matrix hamH(const ring &R, long r, long q);

  /**
   * Holds the parameters of a {@link Ham} code.
   * @see Shannon, C. E. (1948). A mathematical theory of communication
   * @see Hamming, R. W. (1950). Error detecting and error correcting codes
   */
  struct Ham {
    /**
     * The number of parity bits.
     */
    const long r;

    /**
     * The order of the field.
     */
    const long q;

    /**
     * The ring over which the code is defined.
     */
    const ring R;

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
    const long d = 3;

    /**
     * The parity check matrix.
     */
    const matrix H;

    /**
     * The generator matrix.
     */
    const matrix G;

    /**
     * Constructs a {@link Ham} code.
     * @param r The number of parity bits
     * @param q The order of the field
     */
    Ham(const long r, const long q)
        : r(r), q(q), R(NewZZmod(q)), n((SmallPower(q, r) - 1) / (q - 1)), k(n - r),
          H(hamH(R, r, q)), G(genMat(H)) {
    }
  };

  /**
   * Encodes a message using the given {@link Ham} code.
   * @param ham The {@link Ham} code to use
   * @param w The message to encode
   * @return The encoded message
   */
  matrix encodeHam(const Ham &ham, const matrix &w);

  /**
   * Decodes a message using the given {@link Ham} code.
   * @param ham The {@link Ham} code to use
   * @param w The message to decode
   * @return The decoded message
   */
  matrix decodeHam(const Ham &ham, const matrix &w);

}

#endif // ECC_HAM_H
