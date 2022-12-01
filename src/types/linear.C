#include <numeric>
#include "CoCoA/library.H"
#include "types/linear.H"

using namespace std;

namespace CoCoA {

  matrix genMat(const matrix &H) {
    const long n = NumCols(H);
    const long k = n - NumRows(H);

    vector<long> rows(n - k);
    iota(begin(rows), end(rows), 0);
    vector<long> cols(k);
    iota(begin(cols), end(cols), 0);

    const MatrixView P = submat(H, rows, cols);
    return NewDenseMat(ConcatHor(IdentityMat(RingOf(H), k), -transpose(P)));
  }

  matrix checkMat(const matrix &G) {
    const long n = NumCols(G);
    const long k = NumRows(G);

    vector<long> rows(k);
    iota(begin(rows), end(rows), 0);
    vector<long> cols(n - k);
    iota(begin(cols), end(cols), k);

    const MatrixView P = submat(G, rows, cols);
    return NewDenseMat(ConcatHor(-transpose(P), IdentityMat(RingOf(G), n - k)));
  }

  matrix linEncode(const matrix &G, const matrix &w) {
    return w * G;
  }

}
