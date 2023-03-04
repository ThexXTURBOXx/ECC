#include <numeric>
#include "CoCoA/library.H"
#include "ecc/golay.H"
#include "util/utils.H"

using namespace std;

namespace CoCoA {

  matrix G11Mat(const ring &R) {
    return NewDenseMat(R, {{0, 1, 1, 1, 1},
                           {1, 0, 1, 2, 2},
                           {1, 1, 0, 1, 2},
                           {1, 2, 1, 0, 1},
                           {1, 2, 2, 1, 0},
                           {1, 1, 2, 2, 1}});
  }

  matrix G12Mat(const ring &R) {
    const ConstMatrixView G11 = G11Mat(R);
    const RingElem O = zero(R);
    const RingElem I = one(R);
    const RingElem Z = 2 * I;
    return NewDenseMat(ConcatHor(G11, ColMat({I, I, Z, Z, I, O})));
  }

  matrix G23Mat(const ring &R) {
    const RingElem O = zero(R);
    const RingElem I = one(R);
    const ConstMatrixView circ = revCirculantMatrix({I, I, O, I, I, I, O, O, O, I, O});
    return NewDenseMat(ConcatVer(circ, RowMat(vector<RingElem>(11, I))));
  }

  matrix G24Mat(const ring &R) {
    const matrix G23 = G23Mat(R);
    return NewDenseMat(ConcatHor(G23,
                                 ConcatVer(ColMat(vector<RingElem>(11, one(R))), ColMat({zero(R)}))));
  }

  matrix golayMatPart(const ring &R, const long n) {
    switch (n) {
      case 11:
        return G11Mat(R);
      case 12:
        return G12Mat(R);
      case 23:
        return G23Mat(R);
      case 24:
        return G24Mat(R);
      default:
        CoCoA_THROW_ERROR(ERR::BadArg, __func__);
    }
  }

  matrix golayMatPart(const long n) {
    return golayMatPart(NewZZmod(n > 20 ? 2:3), n);
  }

  matrix golayMat(const ring &R, const long n) {
    const matrix G = golayMatPart(R, n);
    return NewDenseMat(ConcatHor(IdentityMat(R, NumRows(G)), G));
  }

  matrix golayMat(const long n) {
    return golayMat(NewZZmod(n > 20 ? 2:3), n);
  }

  matrix encodeGolay(const Golay &gol, const matrix &w) {
    return linEncode(gol.G, w);
  }

  matrix decodeG12(const Golay &gol, const matrix &w) {
    CoCoA_THROW_ERROR(ERR::NYI, __func__);
  }

  matrix decodeG11(const Golay &gol, const matrix &w) {
    CoCoA_THROW_ERROR(ERR::NYI, __func__);
  }

  // See D.G. Hoffman
  matrix decodeG24(const Golay &gol, const matrix &w) {
    const ConstMatrixView S = w * transpose(gol.GExt);
    if (wt(S) <= 3)
      return w + ConcatHor(S, ZeroMat(gol.R, 1, gol.k));

    for (long j = 0; j < 12; j++) {
      const ConstMatrixView aj = RowMat(gol.AExt, j);
      if (wt(S + aj) <= 2)
        return w + ConcatHor(S + aj, e(gol.R, j, one(gol.R), gol.k));
    }

    const ConstMatrixView SA = S * transpose(gol.AExt);
    if (wt(SA) <= 3)
      return w + ConcatHor(ZeroMat(gol.R, 1, gol.k), SA);

    for (long j = 0; j < 12; j++) {
      const ConstMatrixView aj = RowMat(gol.AExt, j);
      if (wt(SA + aj) <= 2)
        return w + ConcatHor(e(gol.R, j, one(gol.R), gol.k), SA + aj);
    }

    CoCoA_THROW_ERROR("Cannot decode!", "Golay");
  }

  matrix decodeG23(const Golay &gol, const matrix &w) {
    const matrix w24 = NewDenseMat(ConcatHor(w, RowMat({RingElem(gol.R, IsEven(wt(w)) ? 1:0)})));

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
    }
  }

}
