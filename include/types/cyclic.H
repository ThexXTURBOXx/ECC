#ifndef ECC_CYCLIC_H
#define ECC_CYCLIC_H

#include "CoCoA/library.H"

using namespace std;

namespace CoCoA {

    RingElem sysEncodeCyclic(ConstRefRingElem g, ConstRefRingElem p, ConstRefRingElem x, long n, long k);

}

#endif //ECC_CYCLIC_H
