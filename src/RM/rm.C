#include "CoCoA/library.H"
#include "ecc/rm.H"

using namespace std;

namespace CoCoA {

    // This function is much easier when it stays recursive.
    matrix RM(const ring &Fp, const int r, const int m) { // NOLINT(misc-no-recursion)
        if (r == 0) {
            return NewDenseMat(RowMat(vector<RingElem>(SmallPower(2, m), one(Fp))));
        }

        if (r == m) {
            const ConstMatrixView R = RM(Fp, m - 1, m);
            const long l = NumCols(R);
            vector<RingElem> v(l, zero(Fp));
            v[l - 1] = one(Fp);
            return NewDenseMat(ConcatVer(R, RowMat(v)));
        }

        const ConstMatrixView R = RM(Fp, r, m - 1);
        const ConstMatrixView S = RM(Fp, r - 1, m - 1);
        const long c = NumCols(R);
        const long l = NumRows(S);
        return NewDenseMat(BlockMat2x2(R, R, ZeroMat(Fp, l, c), S));
    }

    matrix RM(const int r, const int m) {
        return RM(NewZZmod(2), r, m);
    }

}
