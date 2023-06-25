#ifndef ECC_LINEAR_H
#define ECC_LINEAR_H

#include "CoCoA/library.H"

using namespace std;

namespace CoCoA {
  namespace ECC {

    /**
     * Calculates the generator matrix of a linear code from its parity check matrix.
     * @param H The parity check matrix
     * @return The generator matrix
     */
    matrix genMat(const matrix &H);

    /**
     * Calculates the parity check matrix of a linear code from its generator matrix.
     * @param G The generator matrix
     * @return The parity check matrix
     */
    matrix checkMat(const matrix &G);

    /**
     * Encodes a message using a generator matrix.
     * @param G The generator matrix
     * @param w The message
     * @return The encoded message
     */
    matrix linEncode(const matrix &G, const matrix &w);

  }
}

#endif // ECC_LINEAR_H
