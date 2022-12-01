#ifndef ECC_LINEAR_H
#define ECC_LINEAR_H

#include "CoCoA/library.H"

using namespace std;

namespace CoCoA {

  matrix genMat(const matrix &H);

  matrix checkMat(const matrix &G);

  matrix linEncode(const matrix &G, const matrix &w);

}

#endif // ECC_LINEAR_H
