#ifndef ECC_HAM_H
#define ECC_HAM_H

#include "CoCoA/library.H"
#include "types/linear.H"

using namespace std;

namespace CoCoA {

  matrix hamH(const ring &R, long r, long q);

  struct Ham {
    const long r;
    const long q;
    const ring R;
    const long n;
    const long k;
    const long d = 3;
    const matrix H;
    const matrix G;

    Ham(const long r, const long q)
        : r(r), q(q), R(NewZZmod(q)), n((SmallPower(q, r) - 1) / (q - 1)), k(n - r),
          H(hamH(R, r, q)), G(genMat(H)) {
    }
  };

  matrix encodeHam(const Ham &ham, const matrix &w);

  matrix decodeHam(const Ham &ham, const matrix &w);

}

#endif // ECC_HAM_H
