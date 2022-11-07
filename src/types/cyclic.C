#include "CoCoA/library.H"
#include "types/cyclic.H"

using namespace std;

namespace CoCoA {

    RingElem sysEncodeCyclic(ConstRefRingElem g, ConstRefRingElem p, ConstRefRingElem x, const long n, const long k) {
        const RingElem px = p * power(x, n - k);
        const RingElem r = NR(px, {g});
        return px - r;
    }

}
