#include <cstring>
#include "CoCoA/library.H"
#include "util/utils.H"
#include "ecc/bch.H"
#include "ecc/golay.H"
#include "ecc/ham.H"
#include "ecc/rm.H"

using namespace std;

//----------------------------------------------------------------------
const string description = "This file provides simple examples of ECCs.\n";
//----------------------------------------------------------------------

namespace CoCoA {

  void testBCHECC() {
    cout << "========================== BCH ==========================" << endl
         << endl;

    // Setup (BCH code from QR-codes)
    long q = 2;
    long d = 7;
    long c = 1;
    BCH bch = constructBCH(q, d, c, "alpha^4+alpha+1", "alpha", "x");
    cout << bch.g << endl;

    // Encode

    // Example word to encode
    RingElem p = toPolynomial("11011", bch.x);

    // Non-systematic encoding
    cout << toString(bch.g * p, bch.n, bch.x) << endl;

    // Systematic encoding
    cout << toString(p * power(bch.x, bch.n - bch.k), bch.n, bch.x) << endl; // Shifted original poly
    RingElem sent = encodeBCH(bch, p);
    cout << toString(sent, bch.n, bch.x) << endl;

    // Decoding

    RingElem recv1 = sent;
    RingElem recv2 = sent + power(bch.x, 2);
    RingElem recv3 = sent + power(bch.x, 13) + power(bch.x, 5);
    RingElem recv4 = sent + power(bch.x, 8) + power(bch.x, 2) + 1;

    cout << "---" << endl;

    RingElem dec1 = decodeBCH(bch, recv1);
    RingElem dec2 = decodeBCH(bch, recv2);
    RingElem dec3 = decodeBCH(bch, recv3);
    RingElem dec4 = decodeBCH(bch, recv4);

    cout << toString(dec1, bch.n, bch.x) << endl;
    cout << toString(dec2, bch.n, bch.x) << endl;
    cout << toString(dec3, bch.n, bch.x) << endl;
    cout << toString(dec4, bch.n, bch.x) << endl;

    dec1 = decodeBCHGroebner(bch, recv1);
    dec2 = decodeBCHGroebner(bch, recv2);
    dec3 = decodeBCHGroebner(bch, recv3);
    dec4 = decodeBCHGroebner(bch, recv4);

    cout << toString(dec1, bch.n, bch.x) << endl;
    cout << toString(dec2, bch.n, bch.x) << endl;
    cout << toString(dec3, bch.n, bch.x) << endl;
    cout << toString(dec4, bch.n, bch.x) << endl;

    cout << "=============================" << endl;

    // 97-Article%20Text-328-1-10-20180907.pdf
    // Setup (Example ternary BCH Code)

    q = 3;
    d = 5;
    c = 1;
    BCH bch2 = constructBCH(q, d, c, "alpha^2+alpha+2", "alpha", "x");
    cout << bch2.g << endl;

    // Encode

    // Example word to encode
    p = toPolynomial("201", bch2.x);
    sent = encodeBCH(bch2, p);
    cout << toString(sent, bch2.n, bch2.x) << endl;

    // Decoding

    recv1 = sent;
    recv2 = sent + 2 * power(bch2.x, 2);
    recv3 = sent + power(bch2.x, 5) + 2;

    cout << "---" << endl;

    dec1 = decodeBCH(bch2, recv1);
    dec2 = decodeBCH(bch2, recv2);
    dec3 = decodeBCH(bch2, recv3);

    cout << toString(dec1, bch2.n, bch2.x) << endl;
    cout << toString(dec2, bch2.n, bch2.x) << endl;
    cout << toString(dec3, bch2.n, bch2.x) << endl;

    dec1 = decodeBCHGroebner(bch2, recv1);
    dec2 = decodeBCHGroebner(bch2, recv2);
    dec3 = decodeBCHGroebner(bch2, recv3);

    cout << toString(dec1, bch2.n, bch2.x) << endl;
    cout << toString(dec2, bch2.n, bch2.x) << endl;
    cout << toString(dec3, bch2.n, bch2.x) << endl;

    cout << "=============================" << endl;

    // Setup (Example ternary code)

    q = 3;
    d = 7;
    c = 1;
    BCH bch3 = constructBCH(q, d, c, "alpha^3+2*alpha+1", "alpha", "x");
    cout << bch3.g << endl;

    // Encode

    // Example word to encode
    p = toPolynomial("2011211221102", bch3.x);
    sent = encodeBCH(bch3, p);
    cout << toString(sent, bch3.n, bch3.x) << endl;

    // Decoding

    recv1 = sent;
    recv2 = sent + 2 * power(bch3.x, 2);
    recv3 = sent + power(bch3.x, 5) + 2;
    recv4 = sent + 2 * power(bch3.x, 7) + power(bch3.x, 4) + 1;

    cout << "---" << endl;

    dec1 = decodeBCH(bch3, recv1);
    dec2 = decodeBCH(bch3, recv2);
    dec3 = decodeBCH(bch3, recv3);
    dec4 = decodeBCH(bch3, recv4);

    cout << toString(dec1, bch3.n, bch3.x) << endl;
    cout << toString(dec2, bch3.n, bch3.x) << endl;
    cout << toString(dec3, bch3.n, bch3.x) << endl;
    cout << toString(dec4, bch3.n, bch3.x) << endl;

    dec1 = decodeBCHGroebner(bch3, recv1);
    dec2 = decodeBCHGroebner(bch3, recv2);
    dec3 = decodeBCHGroebner(bch3, recv3);
    dec4 = decodeBCHGroebner(bch3, recv4);

    cout << toString(dec1, bch3.n, bch3.x) << endl;
    cout << toString(dec2, bch3.n, bch3.x) << endl;
    cout << toString(dec3, bch3.n, bch3.x) << endl;
    cout << toString(dec4, bch3.n, bch3.x) << endl;
  }

  void testGolayECC() {
    cout << "========================== Golay ==========================" << endl
         << endl;

    const Golay g11(11);
    cout << "[" << g11.n << "," << g11.k << "," << g11.d << "]" << endl;
    cout << g11.G << endl;
    const matrix w1 = toMatrix("101201", g11.R);
    const matrix u1 = encodeGolay(g11, w1);
    cout << toString(u1) << endl;
    /*cout << toString(decodeGolay(g11, u1)) << endl;
    cout << toString(decodeGolay(g11, u1 + toMatrix("10000000000", g11.R))) << endl;
    cout << toString(decodeGolay(g11, u1 + toMatrix("00010020000", g11.R))) << endl;
    cout << toString(decodeGolay(g11, u1 + toMatrix("00020010020", g11.R))) << endl;
    cout << toString(decodeGolay(g11, u1 + toMatrix("00000011200", g11.R))) << endl;*/

    const Golay g12(12);
    cout << "[" << g12.n << "," << g12.k << "," << g12.d << "]" << endl;
    cout << g12.G << endl;
    const matrix w2 = toMatrix("101201", g12.R);
    const matrix u2 = encodeGolay(g12, w2);
    cout << toString(u2) << endl;
    /*cout << toString(decodeGolay(g12, u2)) << endl;
    cout << toString(decodeGolay(g12, u2 + toMatrix("100000000000", g12.R))) << endl;
    cout << toString(decodeGolay(g12, u2 + toMatrix("000100200000", g12.R))) << endl;
    cout << toString(decodeGolay(g12, u2 + toMatrix("000200100200", g12.R))) << endl;
    cout << toString(decodeGolay(g12, u2 + toMatrix("000000112000", g12.R))) << endl;*/

    const Golay g23(23);
    cout << "[" << g23.n << "," << g23.k << "," << g23.d << "]" << endl;
    cout << g23.G << endl;
    const matrix w3 = toMatrix("101100101101", g23.R);
    const matrix u3 = encodeGolay(g23, w3);
    cout << toString(u3) << endl;
    cout << toString(decodeGolay(g23, u3)) << endl;
    cout << toString(decodeGolay(g23, u3 + toMatrix("10000000000000000000000", g23.R))) << endl;
    cout << toString(decodeGolay(g23, u3 + toMatrix("00000010000000100000000", g23.R))) << endl;
    cout << toString(decodeGolay(g23, u3 + toMatrix("00000000000011100000000", g23.R))) << endl;
    cout << toString(decodeGolay(g23, u3 + toMatrix("00000000000111000000000", g23.R))) << endl;

    const Golay g24(24);
    cout << "[" << g24.n << "," << g24.k << "," << g24.d << "]" << endl;
    cout << g24.G << endl;
    const matrix w4 = toMatrix("101100101101", g24.R);
    const matrix u4 = encodeGolay(g24, w4);
    cout << toString(u4) << endl;
    cout << toString(decodeGolay(g24, u4)) << endl;
    cout << toString(decodeGolay(g24, u4 + toMatrix("100000000000000000000000", g24.R))) << endl;
    cout << toString(decodeGolay(g24, u4 + toMatrix("000000100000001000000000", g24.R))) << endl;
    cout << toString(decodeGolay(g24, u4 + toMatrix("000000000000111000000000", g24.R))) << endl;
    cout << toString(decodeGolay(g24, u4 + toMatrix("000000000001110000000000", g24.R))) << endl;
  }

  void testHamECC() {
    cout << "========================== Ham ==========================" << endl
         << endl;

    const Ham ham(3, 2);
    cout << ham.H << endl;
    cout << ham.G << endl;
    const matrix w = toMatrix("0001", ham.R);
    const matrix u = encodeHam(ham, w);
    cout << toString(u) << endl;
    cout << toString(decodeHam(ham, u)) << endl;
    cout << toString(decodeHam(ham, u + toMatrix("0010000", ham.R))) << endl;
    cout << toString(decodeHam(ham, u + toMatrix("0000100", ham.R))) << endl;
    cout << toString(decodeHam(ham, u + toMatrix("0100010", ham.R))) << endl;

    const Ham ham2(3, 3);
    cout << ham2.H << endl;
    cout << ham2.G << endl;
    const matrix w2 = toMatrix("0210210000", ham2.R);
    const matrix u2 = encodeHam(ham2, w2);
    cout << toString(u2) << endl;
    cout << toString(decodeHam(ham2, u2)) << endl;
    cout << toString(decodeHam(ham2, u2 + toMatrix("0020000000000", ham2.R))) << endl;
    cout << toString(decodeHam(ham2, u2 + toMatrix("0000002000000", ham2.R))) << endl;
    cout << toString(decodeHam(ham2, u2 + toMatrix("0000000100000", ham2.R))) << endl;
    cout << toString(decodeHam(ham2, u2 + toMatrix("0100000000020", ham2.R))) << endl;

    const Ham ham3(2, 5);
    cout << ham3.H << endl;
    cout << ham3.G << endl;
    const matrix w3 = toMatrix("3410", ham3.R);
    const matrix u3 = encodeHam(ham3, w3);
    cout << toString(u3) << endl;
    cout << toString(decodeHam(ham3, u3)) << endl;
    cout << toString(decodeHam(ham3, u3 + toMatrix("300000", ham3.R))) << endl;
    cout << toString(decodeHam(ham3, u3 + toMatrix("000001", ham3.R))) << endl;
    cout << toString(decodeHam(ham3, u3 + toMatrix("004000", ham3.R))) << endl;
    cout << toString(decodeHam(ham3, u3 + toMatrix("002030", ham3.R))) << endl;
  }

  void testRMECC() {
    cout << "========================== RM ==========================" << endl
         << endl;

    RM rm(2, 4);
    const matrix w = toMatrix("10011011001", rm.R);
    const matrix u = encodeRM(rm, w);
    cout << toString(u) << endl;
    cout << toString(decodeRM(rm, u)) << endl;
    cout << toString(decodeRM(rm, u + toMatrix("0000000000100000", rm.R))) << endl;
    cout << toString(decodeRM(rm, u + toMatrix("0001000000000000", rm.R))) << endl;
    cout << toString(decodeRM(rm, u + toMatrix("0000000000000100", rm.R))) << endl;
    cout << toString(decodeRM(rm, u + toMatrix("0000000000000001", rm.R))) << endl;

    RM rm2(2, 5);
    const matrix w2 = toMatrix("1001101100101100", rm2.R);
    const matrix u2 = encodeRM(rm2, w2);
    cout << toString(u2) << endl;
    cout << toString(decodeRM(rm2, u2)) << endl;
    cout << toString(decodeRM(rm2, u2 + toMatrix("00000000000000000000000000100000", rm2.R))) << endl;
    cout << toString(decodeRM(rm2, u2 + toMatrix("00010000000000000000000010000000", rm2.R))) << endl;
    cout << toString(decodeRM(rm2, u2 + toMatrix("00000000001000010000000000000001", rm2.R))) << endl;
    cout << toString(decodeRM(rm2, u2 + toMatrix("10000001000000000000000000010000", rm2.R))) << endl;
  }

  void exampleECC(const int argc, const char *argv[]) {
    bool all = false, bch = false, golay = false, ham = false, rm = false;

    if (argc <= 0) {
      all = true;
    } else {
      for (int i = 0; i < argc; ++i) {
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
      testBCHECC();
    }
    if (all || golay) {
      cout << endl
           << endl
           << endl;
      testGolayECC();
    }
    if (all || ham) {
      cout << endl
           << endl
           << endl;
      testHamECC();
    }
    if (all || rm) {
      cout << endl
           << endl
           << endl;
      testRMECC();
    }
  }

}
