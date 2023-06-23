#ifndef ECC_BCH_H
#define ECC_BCH_H

#include "CoCoA/library.H"

using namespace std;

namespace CoCoA {

  /**
   * Holds the parameters of a BCH code.
   * @see Hocquenghem, A. (1959). Codes correcteurs d'erreurs
   * @see Bose, R. C., & Ray-Chaudhuri, D. K. (1960). On a class of error correcting binary group codes
   */
  struct BCH {
    /**
     * The order of the field.
     */
    const long q;

    /**
     * The size of the field.
     */
    const long qn;

    /**
     * The length of the code.
     */
    const long n;

    /**
     * The dimension of the code.
     */
    const long k;

    /**
     * The (designed) distance of the code.
     */
    const long d;

    /**
     * The sense of the code.
     */
    const long c;

    /**
     * The underlying primitive element.
     */
    const RingElem a;

    /**
     * The generator polynomial.
     */
    const RingElem g;

    /**
     * The variable.
     */
    const RingElem x;

    /**
     * Constructs a {@link BCH} code.
     * @param q The order of the field
     * @param qn The size of the field
     * @param n The length of the code
     * @param k The dimension of the code
     * @param d The (designed) distance of the code
     * @param c The sense of the code
     * @param a The underlying primitive element
     * @param g The generator polynomial
     * @param x The variable
     */
    BCH(const long q, const long qn, const long n, const long k,
        const long d, const long c, RingElem a, RingElem g, RingElem x)
        : q(q), qn(qn), n(n), k(k), d(d), c(c),
          a(std::move(a)), g(std::move(g)), x(std::move(x)) {
    }
  };

  /**
   * Constructs a {@link BCH} code based on the given parameters.
   * @param q The order of the field
   * @param d The (designed) distance of the code
   * @param c The sense of the code
   * @param prim A primitive polynomial, dependent on "alpha"
   * @param alpha The primitive element
   * @param x The variable
   * @return The constructed BCH code
   * @see `prim` could be generated using, e.g. Saxena, N. R., & Mccluskey, E. J. (2004). Primitive polynomial
   * generation algorithms implementation and performance analysis
   */
  BCH constructBCH(long q, long d, long c, const string &prim, const string &alpha, const string &x);

  /**
   * Encodes the given polynomial using the given {@link BCH} code.
   * @param bch The {@link BCH} code
   * @param p The polynomial to encode
   * @return The encoded polynomial
   */
  RingElem encodeBCH(const BCH &bch, ConstRefRingElem p);

  /**
   * Decodes the given polynomial using the given {@link BCH} code and the {@link PetersonGorensteinZierler} algorithm.
   * @param bch The {@link BCH} code
   * @param p The polynomial to decode
   * @return The decoded polynomial
   */
  RingElem decodeBCH(const BCH &bch, ConstRefRingElem p);

  /**
   * Decodes the given polynomial using the given {@link BCH} code and the Gröbner basis algorithm.
   * @param bch The {@link BCH} code
   * @param p The polynomial to decode
   * @return The decoded polynomial
   */
  RingElem decodeBCHGroebner(const BCH &bch, ConstRefRingElem p);

}

#endif // ECC_BCH_H
