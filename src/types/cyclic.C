#include "CoCoA/library.H"
#include "types/cyclic.H"
#include "util/utils.H"

using namespace std;

namespace CoCoA {
  namespace ECC {
    RingElem dualPoly(ConstRefRingElem g, const long n) {
      const ring P = owner(g);
      const RingElem h = (IndetPower(P, UnivariateIndetIndex(g), n) - one(P)) / g;
      return reverse(h);
    }

    RingElem sysEncodeCyclic(ConstRefRingElem g, ConstRefRingElem p, ConstRefRingElem x, const long n, const long k) {
      const RingElem px = p * power(x, n - k);
      const RingElem r = NR(px, {g});
      return px - r;
    }

    RingElem decodeCyclicGroebner(ConstRefRingElem g, ConstRefRingElem p, ConstRefRingElem x, ConstRefRingElem a,
                                  const long q, const long n, const long qn) {
      const ring& Px = owner(p);
      const ring& P = CoeffRing(Px);
      const vector<long> J = ChienSearch(g, a, qn, x); // Complete defining set
      const size_t js = J.size();
      const RingElem aP = EvalHom(Px, one(P))(a);

      // Calculate syndromes
      vector<RingElem> s(js, zero(P));
      for (long j = 0; j < js; ++j) {
        const RingHom eval = EvalHom(Px, power(aP, J[j]));
        s[j] = eval(p);
      }

      // All defining syndromes zero -> no error
      if (all_of(s.cbegin(), s.cend(),
                 [](const RingElem& r) {
                   return IsZero(r);
                 }))
        return p;

      long v = 1;
      ideal I({one(Px)});
      do {
        vector<RingElem> S;

        if (q == 2) {
          const ring Rx = NewPolyRing(P, NewSymbols(v), xel);
          const RingElem& Rx1 = one(Rx);

          for (long j = 0; j < js; ++j) {
            RingElem e = RingElem(Rx, -s[j]);
            for (long m = 0; m < v; ++m) {
              e += IndetPower(Rx, m, J[j]);
            }
            S.push_back(e);
          }
          for (long m = 0; m < v; ++m) {
            S.push_back(IndetPower(Rx, m, n) - Rx1);
          }
        } else {
          const ring Rx = NewPolyRing(P, NewSymbols(2 * v), xel);
          const RingElem& Rx1 = one(Rx);

          for (long j = 0; j < js; ++j) {
            RingElem e = RingElem(Rx, -s[j]);
            for (long m = 0; m < v; ++m) {
              e += indet(Rx, v + m) * IndetPower(Rx, m, J[j]);
            }
            S.push_back(e);
          }
          for (long m = 0; m < v; ++m) {
            S.push_back(IndetPower(Rx, v + m, q) - indet(Rx, v + m));
            S.push_back(IndetPower(Rx, m, n) - Rx1);
          }
        }

        I = ideal(S);
        ++v;
      } while (IsOne(I));

      vector<RingElem> G = ReducedGBasis(I);
      const ring Rx = RingOf(I);
      const RingElem aR = RingElem(Rx, aP);
      const RingElem x1 = indet(Rx, 0);
      const RingElem gx1 = getUniPoly(G, 0, IndetPower(Rx, 0, v + 1));

      RingElem f = zero(Px);
      if (deg(gx1) > v)
        CoCoA_THROW_ERROR("Cannot decode!", __func__);

      const vector<long> roots = ChienSearch(gx1, aP, qn, x1);
      if (q == 2) {
        for (auto& j : roots) {
          f += power(x, j);
        }
      } else {
        vector<RingElem> evalPts;
        for (long i = 0; i < v - 1; ++i) {
          evalPts.push_back(power(aR, roots[i]));
        }
        for (long i = v - 1; i < 2 * (v - 1); ++i) {
          evalPts.push_back(indet(Rx, i));
        }
        const RingHom eval = PolyAlgebraHom(Rx, Rx, evalPts);
        for_each(G.begin(), G.end(),
                 [eval](auto& g) {
                   g = eval(g);
                 });
        // TODO: Do we need another reduction here, i.e., G = ReducedGBasis(ideal(G)); ???
        for (long m = 0; m < v - 1; ++m) {
          const long r = ChienSearchSingleRoot(getUniPoly(G, (v - 1) + m, one(Rx)),
                                               aR, qn, indet(Rx, (v - 1) + m));
          f += power(a, r) * power(x, roots[m]);
        }
      }

      return p - f;
    }
  }
}
