#include <numeric>
#include "CoCoA/library.H"
#include "ecc/ham.H"
#include "util/utils.H"

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

    bool sortHam(const vector<long> &a, const vector<long> &b) {
        const long sa = accumulate(a.cbegin(), a.cend(), 0L,
                                   [](const long a, const long b) { return a + sign(b); });
        const long sb = accumulate(b.cbegin(), b.cend(), 0L,
                                   [](const long a, const long b) { return a + sign(b); });
        if (sa != sb) return sa < sb;

        const size_t len = a.size();
        for (int i = 0; i < len; i++) {
            if (a[i] != b[i]) return a[i] < b[i];
        }
        return false;
    }

    vector<RingElem> mult(const long a, const vector<RingElem> &L) {
        vector<RingElem> ret;
        ret.reserve(L.size());
        transform(L.cbegin(), L.cend(), back_inserter(ret),
                  [a](const auto &l) { return a * l; });
        return ret;
    }

    long isMultiple(const long q, const vector<RingElem> &a, const vector<RingElem> &b) {
        for (long i = 1; i <= q; i++) {
            if (mult(i, a) == b) return i;
        }
        return 0;
    }

    matrix hamH(const ring &R, const long r, const long q) {
        vector<long> elems(q);
        iota(begin(elems), end(elems), 0);

        vector<vector<long>> tup = tuples(elems, r);
        sort(tup.begin(), tup.end(), sortHam);

        vector<vector<RingElem>> cols;
        for (const auto &next: tup) {
            if (any_of(next.cbegin(), next.cend(),
                       [](const auto &c) { return c != 0; })) {

                vector<RingElem> n;
                n.reserve(r);
                transform(next.cbegin(), next.cend(), back_inserter(n),
                          [R](const auto &t) { return RingElem(R, t); });

                if (none_of(cols.cbegin(), cols.cend(),
                            [q, n](const auto &c) { return isMultiple(q, c, n); })) {
                    cols.push_back(n);
                }
            }
        }

        reverse(cols.begin(), cols.end());
        return NewDenseMatTranspose(R, cols);
    }

    matrix encodeHam(const Ham &ham, const matrix &w) {
        return w * ham.G;
    }

    matrix e(const Ham &ham, const long i, const RingElem &b) {
        matrix m = NewDenseMat(ZeroMat(ham.R, 1, ham.n));
        SetEntry(m, 0, i, b);
        return m;
    }

    matrix decodeHam(const Ham &ham, const matrix &w) {
        const matrix S = ham.H * transpose(w);
        if (IsZero(S)) return w;

        const vector<RingElem> Svec = GetCol(S, 0);
        for (long i = 0; i < ham.n; i++) {
            const long b = isMultiple(ham.q, GetCol(ham.H, i), Svec);
            if (b != 0) return w - e(ham, i, RingElem(ham.R, b));
        }
        CoCoA_THROW_ERROR("Cannot decode!", "Ham");
        return w; // Shut up compiler warnings
    }

}
