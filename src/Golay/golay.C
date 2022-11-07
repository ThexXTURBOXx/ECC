#include <numeric>
#include "CoCoA/library.H"
#include "ecc/golay.H"
#include "util/utils.H"

using namespace std;

namespace CoCoA {

    matrix G11Mat() {
        const ring R = NewZZmod(3);
        return NewDenseMat(R, {
                {0, 1, 1, 1, 1},
                {1, 0, 1, 2, 2},
                {1, 1, 0, 1, 2},
                {1, 2, 1, 0, 1},
                {1, 2, 2, 1, 0},
                {1, 1, 2, 2, 1}
        });
    }

    matrix G12Mat() {
        const ConstMatrixView G11 = G11Mat();
        const ring &R = RingOf(G11);
        const RingElem O = zero(R);
        const RingElem I = one(R);
        const RingElem Z = 2 * I;
        return NewDenseMat(ConcatHor(G11, ColMat({I, I, Z, Z, I, O})));
    }

    matrix G23Mat() {
        const ring R = NewZZmod(2);
        const RingElem O = zero(R);
        const RingElem I = one(R);
        const ConstMatrixView circ = revCirculantMatrix({I, I, O, I, I, I, O, O, O, I, O});
        return NewDenseMat(ConcatVer(circ, RowMat(vector<RingElem>(11, I))));
    }

    matrix G23toG24(const matrix &G23) {
        const ring &R = RingOf(G23);
        return NewDenseMat(ConcatHor(G23,
                                     ConcatVer(ColMat(vector<RingElem>(11, one(R))), ColMat({zero(R)}))));
    }

    matrix G24Mat() {
        return G23toG24(G23Mat());
    }

    matrix golayMatPart(const long n) {
        switch (n) {
            case 11:
                return G11Mat();
            case 12:
                return G12Mat();
            case 23:
                return G23Mat();
            case 24:
                return G24Mat();
            default:
                CoCoA_THROW_ERROR(ERR::BadArg, __func__);
                return NewDenseMat(ZeroMat(RingZZ(), 1, 1)); // Shut up compiler warnings
        }
    }

    matrix golayMat(const long n) {
        const matrix G = golayMatPart(n);
        return NewDenseMat(ConcatHor(IdentityMat(RingOf(G), NumRows(G)), G));
    }

    matrix encodeGolay(const Golay &gol, const matrix &w) {
        return linEncode(gol.G, w);
    }

    matrix decodeG12(const Golay &gol, const matrix &w) {
        CoCoA_THROW_ERROR(ERR::NYI, __func__);
        return w; // Shut up compiler warnings
    }

    matrix decodeG11(const Golay &gol, const matrix &w) {
        CoCoA_THROW_ERROR(ERR::NYI, __func__);
        return w; // Shut up compiler warnings
    }

    // See D.G. Hoffman
    matrix decodeG24(const Golay &gol, const matrix &w) {
        const bool isExtended = gol.n == 24;
        const matrix G = isExtended ? gol.G : G23toG24(gol.G);

        const ConstMatrixView S = w * transpose(G);
        if (wt(S) <= 3)return w + ConcatHor(S, ZeroMat(gol.R, 1, gol.k));

        const matrix A = isExtended ? gol.A : G23toG24(gol.A);
        for (long j = 0; j < 12; j++) {
            const ConstMatrixView aj = RowMat(A, j);
            if (wt(S + aj) <= 2) return w + ConcatHor(S + aj, e(gol.R, j, one(gol.R), gol.k));
        }

        const ConstMatrixView SA = S * transpose(A);
        if (wt(SA) <= 3) return w + ConcatHor(ZeroMat(gol.R, 1, gol.k), SA);

        for (long j = 0; j < 12; j++) {
            const ConstMatrixView aj = RowMat(A, j);
            if (wt(SA + aj) <= 2) return w + ConcatHor(e(gol.R, j, one(gol.R), gol.k), SA + aj);
        }

        CoCoA_THROW_ERROR("Cannot decode!", "Golay");
        return w; // Shut up compiler warnings
    }

    matrix decodeG23(const Golay &gol, const matrix &w) {
        const matrix w24 = NewDenseMat(ConcatHor(w, RowMat({RingElem(gol.R, IsOdd(wt(w)) ? 0 : 1)})));

        vector<long> cols(23);
        iota(begin(cols), end(cols), 0);
        return NewDenseMat(submat(decodeG24(gol, w24), {0}, cols));
    }

    matrix decodeGolay(const Golay &gol, const matrix &w) {
        switch (gol.n) {
            case 11:
                return decodeG11(gol, w);
            case 12:
                return decodeG12(gol, w);
            case 23:
                return decodeG23(gol, w);
            case 24:
                return decodeG24(gol, w);
            default:
                CoCoA_THROW_ERROR(ERR::BadArg, __func__);
                return w; // Shut up compiler warnings
        }
    }

}
