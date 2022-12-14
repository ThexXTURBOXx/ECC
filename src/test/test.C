#include <cstring>
#include "CoCoA/library.H"
#include "util/utils.H"
#include "ecc/bch.H"
#include "types/cyclic.H"
#include "ecc/golay.H"
#include "ecc/ham.H"
#include "ecc/rm.H"

using namespace std;

//----------------------------------------------------------------------
const string description = "This file provides simple examples of ECCs.\n";
//----------------------------------------------------------------------

namespace CoCoA {

  void testBCH() {
    cout << "========================== BCH ==========================" << endl
         << endl;

    cout << bchGenPoly(1, 7, 2, "alpha^4+alpha+1") << endl;
    cout << bchGenPoly(1, 5, 3, "alpha^2+alpha+2") << endl;
    cout << bchGenPoly(1, 7, 3, "alpha^3+2*alpha+1") << endl;

    // Setup

    // F_(2^q)[x]
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
    BCH bch(q, qn, n, k, d, c, a, g);

    // Encode

    // Example word to encode
    RingElem p = toPolynomial("11011", x);

    // Non-systematic encoding
    cout << toString(g * p, n, x) << endl;

    // Systematic encoding
    cout << toString(p * power(x, n - k), n, x) << endl;
    RingElem sent = encodeBCH(bch, p, x);
    cout << toString(sent, n, x) << endl;

    // Decoding

    RingElem recv1 = sent;
    RingElem recv2 = sent + power(x, 2);
    RingElem recv3 = sent + power(x, 13) + power(x, 5);
    RingElem recv4 = sent + power(x, 8) + power(x, 2) + 1;

    cout << "---" << endl;

    RingElem dec1 = decodeBCH(bch, recv1, x);
    RingElem dec2 = decodeBCH(bch, recv2, x);
    RingElem dec3 = decodeBCH(bch, recv3, x);
    RingElem dec4 = decodeBCH(bch, recv4, x);

    cout << toString(dec1, n, x) << endl;
    cout << toString(dec2, n, x) << endl;
    cout << toString(dec3, n, x) << endl;
    cout << toString(dec4, n, x) << endl;

    dec1 = decodeCyclicGroebner(g, recv1, x, a, q, n, qn);
    dec2 = decodeCyclicGroebner(g, recv2, x, a, q, n, qn);
    dec3 = decodeCyclicGroebner(g, recv3, x, a, q, n, qn);
    dec4 = decodeCyclicGroebner(g, recv4, x, a, q, n, qn);

    cout << toString(dec1, n, x) << endl;
    cout << toString(dec2, n, x) << endl;
    cout << toString(dec3, n, x) << endl;
    cout << toString(dec4, n, x) << endl;

    cout << "=============================" << endl;

    // 97-Article%20Text-328-1-10-20180907.pdf
    // Setup

    Px = NewPolyRing(NewZZmod(3), symbols("alpha"));
    I = ideal(RingElem(Px, "alpha^2+alpha+2"));
    Rx = NewPolyRing(NewQuotientRing(Px, I), symbols("x"));
    x = RingElem(Rx, "x");

    // Example ternary BCH Code
    q = 3;
    qn = SmallPower(3, 2);
    n = 8;
    k = 3;
    a = RingElem(Rx, "alpha");
    d = 5;
    c = 1;
    g = toPolynomial("121102", x);
    BCH bch2(q, qn, n, k, d, c, a, g);

    // Encode

    // Example word to encode
    p = toPolynomial("201", x);
    sent = encodeBCH(bch2, p, x);
    cout << toString(sent, n, x) << endl;

    // Decoding

    recv1 = sent;
    recv2 = sent + 2 * power(x, 2);
    recv3 = sent + power(x, 5) + 2;

    cout << "---" << endl;

    dec1 = decodeBCH(bch2, recv1, x);
    dec2 = decodeBCH(bch2, recv2, x);
    dec3 = decodeBCH(bch2, recv3, x);

    cout << toString(dec1, n, x) << endl;
    cout << toString(dec2, n, x) << endl;
    cout << toString(dec3, n, x) << endl;

    dec1 = decodeCyclicGroebner(g, recv1, x, a, q, n, qn);
    dec2 = decodeCyclicGroebner(g, recv2, x, a, q, n, qn);
    dec3 = decodeCyclicGroebner(g, recv3, x, a, q, n, qn);

    cout << toString(dec1, n, x) << endl;
    cout << toString(dec2, n, x) << endl;
    cout << toString(dec3, n, x) << endl;

    cout << "=============================" << endl;

    // Setup

    Px = NewPolyRing(NewZZmod(3), symbols("alpha"));
    I = ideal(RingElem(Px, "alpha^3+2*alpha+1"));
    Rx = NewPolyRing(NewQuotientRing(Px, I), symbols("x"));
    x = RingElem(Rx, "x");

    // Example ternary BCH Code
    q = 3;
    qn = SmallPower(3, 3);
    n = 26;
    k = 13;
    a = RingElem(Rx, "alpha");
    d = 7;
    c = 1;
    g = toPolynomial("1100002001221", x);
    BCH bch3(q, qn, n, k, d, c, a, g);

    // Encode

    // Example word to encode
    p = toPolynomial("2011211221102", x);
    sent = encodeBCH(bch3, p, x);
    cout << toString(sent, n, x) << endl;

    // Decoding

    recv1 = sent;
    recv2 = sent + 2 * power(x, 2);
    recv3 = sent + power(x, 5) + 2;
    recv4 = sent + 2 * power(x, 7) + power(x, 4) + 1;

    cout << "---" << endl;

    dec1 = decodeBCH(bch3, recv1, x);
    dec2 = decodeBCH(bch3, recv2, x);
    dec3 = decodeBCH(bch3, recv3, x);
    dec4 = decodeBCH(bch3, recv4, x);

    cout << toString(dec1, n, x) << endl;
    cout << toString(dec2, n, x) << endl;
    cout << toString(dec3, n, x) << endl;
    cout << toString(dec4, n, x) << endl;

    dec1 = decodeCyclicGroebner(g, recv1, x, a, q, n, qn);
    dec2 = decodeCyclicGroebner(g, recv2, x, a, q, n, qn);
    dec3 = decodeCyclicGroebner(g, recv3, x, a, q, n, qn);
    dec4 = decodeCyclicGroebner(g, recv4, x, a, q, n, qn);

    cout << toString(dec1, n, x) << endl;
    cout << toString(dec2, n, x) << endl;
    cout << toString(dec3, n, x) << endl;
    cout << toString(dec4, n, x) << endl;
  }

  void testGolay() {
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

  void testHam() {
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

  void testRM() {
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

  void example(const int argc, const char *argv[]) {
    bool all = false, bch = false, golay = false, ham = false, rm = false;

    if (argc <= 1) {
      all = true;
    } else {
      for (int i = 1; i < argc; ++i) {
        const char *arg = argv[i];
        if (strcmp(arg, "bch")==0) {
          bch = true;
        } else if (strcmp(arg, "ham")==0) {
          ham = true;
        } else if (strcmp(arg, "golay")==0) {
          golay = true;
        } else if (strcmp(arg, "rm")==0) {
          rm = true;
        }
      }
    }

    if (all || bch) {
      cout << endl
           << endl
           << endl;
      testBCH();
    }
    if (all || golay) {
      cout << endl
           << endl
           << endl;
      testGolay();
    }
    if (all || ham) {
      cout << endl
           << endl
           << endl;
      testHam();
    }
    if (all || rm) {
      cout << endl
           << endl
           << endl;
      testRM();
    }
  }

  // DO NOT EDIT LINES BELOW HERE

  void program(const int argc, const char *argv[]) {
    GlobalManager CoCoAFoundations(UseNonNegResidues);
    SignalWatcher MonitorInterrupt(SIGINT); // you must also call CheckForInterrupt every so often

    cout << description << endl;
    cout << boolalpha; // so that bools print out as true/false

    example(argc, argv);
  }

}

int main(const int argc, const char *argv[]) {
  try {
    CoCoA::program(argc, argv);
    return 0;
  }
  catch (const CoCoA::InterruptReceived &intr) {
    cerr << endl
         << "------------------------------" << endl
         << ">>>  CoCoALib interrupted  <<<" << endl
         << "------------------------------" << endl
         << "-->>  " << intr << "  <<--" << endl;
    return 2;
  }
  catch (const CoCoA::ErrorInfo &err) {
    cerr << "***ERROR***  UNCAUGHT CoCoA error";
    ANNOUNCE(cerr, err);
  }
  catch (const std::exception &exc) {
    cerr << "***ERROR***  UNCAUGHT std::exception: " << exc.what() << endl;
  }
  catch (...) {
    cerr << "***ERROR***  UNCAUGHT UNKNOWN EXCEPTION" << endl;
  }

  CoCoA::BuildInfo::PrintAll(cerr);
  return 1;
}
