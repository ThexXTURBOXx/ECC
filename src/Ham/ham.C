#include <numeric>
#include "CoCoA/library.H"
#include "ecc/ham.H"
#include "util/utils.H"

using namespace std;

namespace CoCoA {

    bool sortHam(const vector<long> &a, const vector<long> &b) {
        const long wa = wt(a);
        const long wb = wt(b);
        if (wa != wb) return wa < wb;

        const size_t len = a.size();
        for (int i = 0; i < len; ++i) {
            if (a[i] != b[i]) return a[i] < b[i];
        }
        return false;
    }

    bool isMultiple(const vector<RingElem> &a, const vector<RingElem> &b, const long c) {
        for (size_t i = 0; i < a.size(); ++i) {
            if (a[i] != c * b[i]) return false;
        }
        return true;
    }

    long divide(const vector<RingElem> &a, const vector<RingElem> &b, const long q) {
        for (long i = 1; i <= q; ++i) {
            if (isMultiple(a, b, i)) return i;
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
                            [q, n](const auto &c) { return divide(n, c, q); })) {
                    cols.push_back(n);
                }
            }
        }

        reverse(cols.begin(), cols.end());
        return NewDenseMatTranspose(R, cols);
    }

    matrix encodeHam(const Ham &ham, const matrix &w) {
        return linEncode(ham.G, w);
    }

    matrix decodeHam(const Ham &ham, const matrix &w) {
        const matrix S = ham.H * transpose(w);
        if (IsZero(S)) return w;

        const vector<RingElem> Svec = GetCol(S, 0);
        for (long i = 0; i < ham.n; ++i) {
            const long b = divide(Svec, GetCol(ham.H, i), ham.q);
            if (b != 0) return w - e(ham.R, i, RingElem(ham.R, b), ham.n);
        }
        CoCoA_THROW_ERROR("Cannot decode!", "Ham");
        return w; // Shut up compiler warnings
    }

}
