#ifndef ECC_GOLAY_H
#define ECC_GOLAY_H

#include "CoCoA/library.H"
#include "types/linear.H"

using namespace std;

namespace CoCoA {

    matrix golayMatPart(long n);

    matrix golayMat(long n);

    struct Golay {
        const matrix A;
        const ring R;
        const matrix G;
        const long q;
        const long n;
        const long k;
        const long d;

        explicit Golay(const long n)
                : A(golayMatPart(n)), R(RingOf(A)), G(NewDenseMat(ConcatHor(IdentityMat(R, NumRows(A)), A))),
                  q(n > 20 ? 2 : 3), n(n), k((n + 1) / 2), d(n / (5 - q)) /* magic formula */  {
        }
    };

    matrix encodeGolay(const Golay &golay, const matrix &w);

    matrix decodeGolay(const Golay &golay, const matrix &w);

}

#endif //ECC_GOLAY_H
