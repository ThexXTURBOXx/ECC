#ifndef ECC_FUZZY_H
#define ECC_FUZZY_H

#include <functional>
#include <numeric>
#include <utility>
#include <SHA256.h>

#include "CoCoA/library.H"
#include "util/utils.H"

using namespace std;

namespace CoCoA {

  struct HelperData {
    matrix s;
    matrix x;

    HelperData() : HelperData(NewDenseMat(RingZZ(), 0, 0), NewDenseMat(RingZZ(), 0, 0)) {
    }

    HelperData(const matrix &s, const matrix &x) : s(s), x(x) {
    }
  };

  struct FuzzyExtractor {
    using ECCFn = std::function<matrix(const matrix &)>;

    FuzzyExtractor(ECCFn encode, ECCFn decode, long messageBits, long parityBits)
        : encode(std::move(encode)), decode(std::move(decode)), messageBits(messageBits), parityBits(parityBits),
          helperDataSet(false), hd(HelperData()) {
    }

    FuzzyExtractor(ECCFn encode, ECCFn decode, long messageBits, long parityBits, HelperData hd)
        : encode(std::move(encode)), decode(std::move(decode)), messageBits(messageBits), parityBits(parityBits),
          helperDataSet(true), hd(std::move(hd)) {
    }

    matrix generateHelperData(const matrix &w);

    matrix extract(const matrix &wd);

    HelperData getHelperData() {
      return {hd}; // Copy to avoid modifying
    }

  private:
    ECCFn encode;
    ECCFn decode;
    long messageBits;
    long parityBits;
    bool helperDataSet;
    HelperData hd;

    static matrix strongExtract(const matrix &w);
  };

}

#endif //ECC_FUZZY_H
