#include "CoCoA/library.H"
#include "types/cyclic.H"
#include "ecc/bch.H"
#include "util/utils.H"

using namespace std;

namespace CoCoA {

  BCH constructBCH(const long q, const long d, const long c, const string &prim, const string &alpha, const string &x) {
    // F_(q^m)[x]
    const ring P = NewZZmod(q);
    const PolyRing Px = NewPolyRing(P, symbols(alpha));
    RingElem a(Px, alpha);
    const RingElem primPoly = RingElem(Px, prim);
    const ideal I(primPoly);

    RingElem g = one(Px);
    for (long i = c; i <= c + d - 2; ++i)
      g = lcm(g, MinPolyQuot(power(a, i), I, a));

    const PolyRing Rx = NewPolyRing(NewQuotientRing(Px, I), symbols(x));

    const RingElem indet(Rx, x);
    g = PolyRingHom(Px, Rx, ChainCanonicalHom(P, Rx), {indet})(g);
    a = RingElem(Rx, alpha);
    const long qn = SmallPower(q, deg(primPoly));
    const long n = qn - 1; // Assume primitive BCH code, length of codewords in the code
    const long k = n - deg(g); // Amount of information bits in codeword (input length)

    return {q, qn, n, k, d, c, a, g, indet};
  }

  /**
   * Calculates the error locator polynomial using the Peterson-Gorenstein-Zierler algorithm.
   * @param p The generator polynomial
   * @param x The variable
   * @param s The necessary syndromes
   * @param v The number of errors that can be corrected
   * @return The error locator polynomial
   * @see Gorenstein, D., Peterson, W. W., & Zierler, N. (1960). Two-error correcting Bose-Chaudhuri codes are
   * quasi-perfect
   */
  RingElem PetersonGorensteinZierler(ConstRefRingElem p, ConstRefRingElem x, const vector<RingElem> &s, long v) {
    const ring &Px = owner(p);

    if (all_of(s.cbegin(), s.cend(),
               [](const RingElem &r) {
                 return IsZero(r);
               }))
      return zero(Px);

    matrix M = NewDenseMat(Px, v, v);
    for (long i = 0; i < v; ++i)
      for (long j = 0; j < v; ++j)
        SetEntry(M, i, j, s[i + j]);

    bool zd;
    do {
      // Resizing suffices since the entries do not change
      M->myResize(v, v);
      zd = IsZeroDet(M);
      if (zd && v == 1)
        CoCoA_THROW_ERROR("Cannot decode!", "BCH");
      --v;
    } while (zd && v != 0);

    const long w = v + 1;
    matrix V = NewDenseMat(Px, w, 1);
    for (long i = 0; i < w; ++i)
      SetEntry(V, i, 0, s[w + i]);

    // Sadly NYI: const matrix P = LinSolve(M, -V);
    const matrix P = -inverse(M) * V;

    RingElem ret = one(Px);
    for (long i = 0; i < w; ++i)
      ret += P(w - i - 1, 0) * power(x, i + 1);
    return ret;
  }

  /**
   * Calculates the error values using the Forney algorithm.
   * @param bch The {@link BCH} code
   * @param s The list of necessary syndromes
   * @param e The error locator polynomial
   * @param roots The roots of the error locator polynomial
   * @param x The variable
   * @return The error polynomial
   * @see Forney, G. (1965). On decoding BCH codes
   */
  RingElem Forney(const BCH &bch, const vector<RingElem> &s, ConstRefRingElem e, const vector<long> &roots,
                  ConstRefRingElem x) {
    const ring &Px = owner(e);

    RingElem S = zero(Px);
    for (long i = 0; i < s.size(); ++i)
      S += s[i] * power(x, i);

    const RingElem O = NR(S * e, {power(x, bch.d - 1)});
    const RingElem ed = deriv(e, x);

    RingElem ret = zero(Px);
    for (long k: roots) {
      const long revK = (k - bch.qn + 1) % (bch.qn - 1);
      const RingHom eval = PolyAlgebraHom(
          Px, Px, {power(bch.a, revK)});
      ret += -power(x, -revK) * (power(bch.a, -revK) * eval(O)) / (power(bch.a, -bch.c * revK) * eval(ed));
    }
    return ret;
  }

  RingElem encodeBCH(const BCH &bch, ConstRefRingElem p) {
    return sysEncodeCyclic(bch.g, p, bch.x, bch.n, bch.k);
  }

  RingElem decodeBCH(const BCH &bch, ConstRefRingElem p) {
    const ring &Px = owner(p);
    const long t = (bch.d - 1) / 2;

    // Calculate syndromes
    vector<RingElem> s(bch.d - 1, zero(Px));
    for (long j = 0; j <= bch.d - 2; ++j) {
      const RingHom eval = PolyAlgebraHom(
          Px, Px, {power(bch.a, bch.c + j)});
      s[j] = eval(p);
    }

    // Calculate error locator polynomial using the
    // Peterson-Gorenstein-Zierler algorithm
    const RingElem e = PetersonGorensteinZierler(p, bch.x, s, t);
    if (IsZero(e))
      return p;

    // Factor error locator polynomial using Chien Search
    const vector<long> roots = ChienSearch(e, bch.a, bch.qn, bch.x);
    if (roots.empty())
      return p;

    // Calculate error values using Forney's algorithm and correct errors
    if (bch.q == 2) {
      RingElem f = p;
      for (long j: roots)
        f -= power(bch.x, -((j - bch.qn + 1) % (bch.qn - 1)));
      return f;
    } else {
      return p - Forney(bch, s, e, roots, bch.x);
    }
  }

  RingElem decodeBCHGroebner(const BCH &bch, ConstRefRingElem p) {
    return decodeCyclicGroebner(bch.g, p, bch.x, bch.a, bch.q, bch.n, bch.qn);
  }

}
