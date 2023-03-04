#include <cstring>
#include <functional>
#include "CoCoA/library.H"
#include "util/utils.H"
#include "ecc/bch.H"
#include "ecc/golay.H"
#include "ecc/ham.H"
#include "ecc/rm.H"
#include "fuzzy/fuzzy.H"
#include "types/cyclic.H"

using namespace std;
using namespace std::placeholders;

//----------------------------------------------------------------------
const string description = "This file provides simple examples of Fuzzy Extractors.\n";
//----------------------------------------------------------------------

namespace CoCoA {

  void testBCHFuzzy() {
    cout << "========================== BCH ==========================" << endl
         << endl;

    PolyRing Px = NewPolyRing(NewZZmod(2), symbols("alpha"));
    ideal I(RingElem(Px, "alpha^4+alpha+1"));
    PolyRing Rx = NewPolyRing(NewQuotientRing(Px, I), symbols("x"));
    RingElem x = RingElem(Rx, "x");

    // BCH code from QR-codes
    // n = 15: length of codewords in the code
    // k = 5: amount of information bits in codeword (input length)
    long q = 2;
    long qn = SmallPower(q, 4);
    long n = 15;
    long k = 5;
    long d = 7;
    long c = 1;
    RingElem a(Rx, "alpha");
    RingElem g = toPolynomial("10100110111", x);
    BCH bch(q, qn, n, k, d, c, a, g, x);
    RingElem w = toPolynomial("110111001010100", x);

    // std::bind is also a possibility here, but Clang-Tidy complains!
    FuzzyExtractor ext([bch](auto &&w) {
                         return toMatrix(encodeBCH(bch, toPolynomial(std::forward<decltype(w)>(w), bch.x)),
                                         bch.n, bch.x);
                       },
                       [bch](auto &&w) {
                         return toMatrix(decodeBCH(bch, toPolynomial(std::forward<decltype(w)>(w), bch.x)),
                                         bch.n, bch.x);
                       },
                       5, 10);

    matrix original = ext.generateHelperData(toMatrix(w, n, x));
    cout << toString(original) << endl;
    cout << toString(ext.extract(toMatrix(w, n, x))) << endl;
    cout << toString(ext.extract(toMatrix(w + power(x, 2), n, x))) << endl;
    cout << toString(ext.extract(toMatrix(w + power(x, 13) + power(x, 5), n, x))) << endl;
    cout << toString(ext.extract(toMatrix(w + power(x, 8) + power(x, 2) + 1, n, x))) << endl;

    cout << "---" << endl;

    // Fuzzy extractor using Groebner basis decoding

    // std::bind is also a possibility here, but Clang-Tidy complains!
    FuzzyExtractor ext2([bch](auto &&w) {
                          return toMatrix(encodeBCH(bch, toPolynomial(std::forward<decltype(w)>(w), bch.x)),
                                          bch.n, bch.x);
                        },
                        [bch](auto &&w) {
                          return toMatrix(decodeCyclicGroebner(bch.g,
                                                               toPolynomial(std::forward<decltype(w)>(w), bch.x),
                                                               bch.x, bch.a, bch.q, bch.n, bch.qn), bch.n, bch.x);
                        },
                        5, 10, ext.getHelperData());

    cout << toString(original) << endl;
    cout << toString(ext2.extract(toMatrix(w, n, x))) << endl;
    cout << toString(ext2.extract(toMatrix(w + power(x, 2), n, x))) << endl;
    cout << toString(ext2.extract(toMatrix(w + power(x, 13) + power(x, 5), n, x))) << endl;
    cout << toString(ext2.extract(toMatrix(w + power(x, 8) + power(x, 2) + 1, n, x))) << endl;
  }

  void testGolayFuzzy() {
    cout << "========================== Golay ==========================" << endl
         << endl;

    const Golay g24(24);
    const matrix w = toMatrix("101100101101000011101111", g24.R);

    // std::bind is also a possibility here, but Clang-Tidy complains!
    FuzzyExtractor ext([g24](auto &&w) {
                         return encodeGolay(g24, std::forward<decltype(w)>(w));
                       },
                       [g24](auto &&w) {
                         return decodeGolay(g24, std::forward<decltype(w)>(w));
                       },
                       12, 12);

    matrix original = ext.generateHelperData(w);
    cout << toString(original) << endl;
    cout << toString(ext.extract(w)) << endl;
    cout << toString(ext.extract(w + toMatrix("100000000000000000000000", g24.R))) << endl;
    cout << toString(ext.extract(w + toMatrix("000100000000000001000000", g24.R))) << endl;
    cout << toString(ext.extract(w + toMatrix("000000001110000000000000", g24.R))) << endl;
    cout << toString(ext.extract(w + toMatrix("000000000000011100000000", g24.R))) << endl;
  }

  void testHamFuzzy() {
    cout << "========================== Ham ==========================" << endl
         << endl;

    const Ham ham(3, 2);
    const matrix w = toMatrix("0101100", ham.R);

    // std::bind is also a possibility here, but Clang-Tidy complains!
    FuzzyExtractor ext([ham](auto &&w) {
                         return encodeHam(ham, std::forward<decltype(w)>(w));
                       },
                       [ham](auto &&w) {
                         return decodeHam(ham, std::forward<decltype(w)>(w));
                       },
                       4, 3);

    matrix original = ext.generateHelperData(w);
    cout << toString(original) << endl;
    cout << toString(ext.extract(w)) << endl;
    cout << toString(ext.extract(w + toMatrix("0010000", ham.R))) << endl;
    cout << toString(ext.extract(w + toMatrix("0000100", ham.R))) << endl;
    cout << toString(ext.extract(w + toMatrix("0100010", ham.R))) << endl;
  }

  void testRMFuzzy() {
    cout << "========================== RM ==========================" << endl
         << endl;

    RM rm(2, 5);
    const matrix w = toMatrix("11011100101110111001110011111011", rm.R);

    // std::bind is also a possibility here, but Clang-Tidy complains!
    FuzzyExtractor ext([rm](auto &&w) {
                         return encodeRM(rm, std::forward<decltype(w)>(w));
                       },
                       [rm](auto &&w) {
                         return decodeRM(rm, std::forward<decltype(w)>(w));
                       },
                       16, 16);

    matrix original = ext.generateHelperData(w);
    cout << toString(original) << endl;
    cout << toString(ext.extract(w)) << endl;
    cout << toString(ext.extract(w + toMatrix("00000000000000000000000000100000", rm.R))) << endl;
    cout << toString(ext.extract(w + toMatrix("00010000000000000000000010000000", rm.R))) << endl;
    cout << toString(ext.extract(w + toMatrix("00000000001000010000000000000001", rm.R))) << endl;
    cout << toString(ext.extract(w + toMatrix("10000001000000000000000000010000", rm.R))) << endl;
  }

  void exampleFuzzy(const int argc, const char *argv[]) {
    bool all = false, bch = false, golay = false, ham = false, rm = false;

    if (argc <= 1) {
      all = true;
    } else {
      for (int i = 1; i < argc; ++i) {
        const char *arg = argv[i];
        if (strcmp(arg, "bch") == 0) {
          bch = true;
        } else if (strcmp(arg, "ham") == 0) {
          ham = true;
        } else if (strcmp(arg, "golay") == 0) {
          golay = true;
        } else if (strcmp(arg, "rm") == 0) {
          rm = true;
        }
      }
    }

    if (all || bch) {
      cout << endl
           << endl
           << endl;
      testBCHFuzzy();
    }
    if (all || golay) {
      cout << endl
           << endl
           << endl;
      testGolayFuzzy();
    }
    if (all || ham) {
      cout << endl
           << endl
           << endl;
      testHamFuzzy();
    }
    if (all || rm) {
      cout << endl
           << endl
           << endl;
      testRMFuzzy();
    }
  }

}
