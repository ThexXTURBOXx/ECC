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

  /**
   * Holds the helper data for the {@link FuzzyExtractor}.
   */
  struct HelperData {
    /**
     * The first part of the helper data.
     */
    matrix s;

    /**
     * The second part of the helper data.
     */
    matrix x;

    /**
     * Creates a new {@link HelperData} object with empty data.
     */
    HelperData() : HelperData(NewDenseMat(RingZZ(), 0, 0), NewDenseMat(RingZZ(), 0, 0)) {
    }

    /**
     * Creates a new {@link HelperData} object with the given data.
     * @param s The first part of the helper data
     * @param x The second part of the helper data
     */
    HelperData(const matrix &s, const matrix &x) : s(s), x(x) {
    }
  };

  /**
   * Provides a fuzzy extractor for a given ECC.
   */
  struct FuzzyExtractor {
    /**
     * Function signature for encoding and decoding procedures.
     */
    using ECCFn = std::function<matrix(const matrix &)>;

    /**
     * Creates a new {@link FuzzyExtractor} object with the given parameters and no helper data.
     * @param encode The encoding function
     * @param decode The decoding function
     * @param messageBits The number of message bits
     * @param parityBits The number of parity bits
     */
    FuzzyExtractor(ECCFn encode, ECCFn decode, long messageBits, long parityBits)
        : encode(std::move(encode)), decode(std::move(decode)), messageBits(messageBits), parityBits(parityBits),
          helperDataSet(false), hd(HelperData()) {
    }

    /**
     * Creates a new {@link FuzzyExtractor} object with the given parameters and previously generated helper data.
     * @param encode The encoding function
     * @param decode The decoding function
     * @param messageBits The number of message bits
     * @param parityBits The number of parity bits
     * @param hd {@link HelperData} object with the helper data to use
     */
    FuzzyExtractor(ECCFn encode, ECCFn decode, long messageBits, long parityBits, HelperData hd)
        : encode(std::move(encode)), decode(std::move(decode)), messageBits(messageBits), parityBits(parityBits),
          helperDataSet(true), hd(std::move(hd)) {
    }

    /**
     * Generates helper data for the given word, can only be called once.
     * @param w The word to generate helper data for
     * @return The extracted word/protected secret
     */
    matrix generateHelperData(const matrix &w);

    /**
     * Tries to extract the word/protected secret from the given word.
     * @param wd The word to extract the secret from
     * @return The extracted word/protected secret iff it is valid
     */
    matrix extract(const matrix &wd);

    /**
     * Returns a copy of the {@link HelperData} object.
     * @return A copy of the {@link HelperData} object
     */
    HelperData getHelperData() {
      return {hd}; // Copy to avoid modifying
    }

  private:
    /**
     * The encoding function.
     */
    ECCFn encode;

    /**
     * The decoding function.
     */
    ECCFn decode;

    /**
     * The amount of message bits.
     */
    long messageBits;

    /**
     * The amount of parity bits.
     */
    long parityBits;

    /**
     * Whether the {@link HelperData} object has been set.
     */
    bool helperDataSet;

    /**
     * The {@link HelperData} object.
     */
    HelperData hd;

    /**
     * The underlying strong extractor, uses SHA-256.
     */
    static matrix strongExtract(const matrix &w);
  };

}

#endif //ECC_FUZZY_H
