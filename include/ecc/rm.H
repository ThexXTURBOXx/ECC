#ifndef ECC_RM_H
#define ECC_RM_H

#include "CoCoA/library.H"

using namespace std;

namespace CoCoA {

  struct RM {
    const ring R;
    const long r;
    const long m;
    const long n;
    const long k;
    const long d;
    const vector<vector<RingElem>> xrows;
    const matrix G;
    const vector<vector<vector<RingElem>>> votingRows;
    const vector<long> ribd;

    RM(const long r, const long m) : RM(NewZZmod(2), r, m) {
    }

    RM(const ring &R, const long r, const long m) : R(R), r(r), m(m), n(SmallPower(2, m)), k(genK(r, m)),
                                                    d(SmallPower(2, m - r)), xrows(genXrows(R, m)),
                                                    G(genG(R, r, m, xrows)), votingRows(genVotingRows(R, r, m)),
                                                    ribd(genRibd(r, m)) {
    }

  private:
    static long genK(long r, long m);

    static vector<vector<RingElem>> genXrows(const ring &R, long m);

    static matrix genG(const ring &R, long r, long m, const vector<vector<RingElem>> &xrows);

    static vector<vector<vector<RingElem>>> genVotingRows(const ring &R, long r, long m);

    static vector<long> genRibd(long r, long m);
  };

  matrix encodeRM(const RM &rm, const matrix &w);

  matrix decodeRM(const RM &rm, matrix w);

}

#endif // ECC_RM_H
