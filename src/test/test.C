#include <cstring>
#include "CoCoA/library.H"
#include "util/utils.H"
#include "ecc/bch.H"
#include "ecc/ham.H"
#include "ecc/rm.H"

using namespace std;

//----------------------------------------------------------------------
const string description = "This file provides simple examples of ECCs.\n";
//----------------------------------------------------------------------

namespace CoCoA {

    void testBCH() {
        cout << "========================== BCH ==========================" << endl << endl;

        cout << bchGenPoly(1, 7, 2, "alpha^4+alpha+1") << endl;
        cout << bchGenPoly(1, 5, 3, "alpha^2+alpha+2") << endl;

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
        RingElem g = toPolynomial("10100110111", n - k + 1, x);
        BCH bch(q, qn, n, k, d, c, a, g);

        // Encode

        // Example word to encode
        RingElem p = toPolynomial("11011", k, x);

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
        g = toPolynomial("121102", n - k + 1, x);
        BCH bch2(q, qn, n, k, d, c, a, g);

        // Encode

        // Example word to encode
        p = toPolynomial("201", k, x);
        sent = encodeBCH(bch2, p, x);
        cout << toString(sent, n, x) << endl;

        // Decoding

        recv1 = sent;
        recv2 = sent + 2 * power(x, 2);
        recv3 = sent + power(x, 5) + 1;

        cout << "---" << endl;

        dec1 = decodeBCH(bch2, recv1, x);
        dec2 = decodeBCH(bch2, recv2, x);
        dec3 = decodeBCH(bch2, recv3, x);

        cout << toString(dec1, n, x) << endl;
        cout << toString(dec2, n, x) << endl;
        cout << toString(dec3, n, x) << endl;
    }

    void testHam() {
        cout << "========================== Ham ==========================" << endl << endl;

        const Ham ham(3, 2);
        cout << ham.H << endl;
        cout << ham.G << endl;
        const RingElem O = zero(ham.R);
        const RingElem I = one(ham.R);
        const matrix w = NewDenseMat(RowMat({O, O, O, I}));
        const matrix u = encodeHam(ham, w);
        cout << u << endl;
        cout << decodeHam(ham, u) << endl;
        cout << decodeHam(ham, u + RowMat({O, O, I, O, O, O, O})) << endl;
        cout << decodeHam(ham, u + RowMat({O, O, O, O, I, O, O})) << endl;
        cout << decodeHam(ham, u + RowMat({O, I, O, O, O, I, O})) << endl;

        const Ham ham2(3, 3);
        cout << ham2.H << endl;
        cout << ham2.G << endl;
        const RingElem O2 = zero(ham2.R);
        const RingElem I2 = one(ham2.R);
        const RingElem Z2 = 2 * I2;
        const matrix w2 = NewDenseMat(RowMat({O2, Z2, I2, O2, Z2, I2, O2, O2, O2, O2}));
        const matrix u2 = encodeHam(ham2, w2);
        cout << u2 << endl;
        cout << decodeHam(ham2, u2) << endl;
        cout << decodeHam(ham2, u2 + RowMat({O2, O2, Z2, O2, O2, O2, O2, O2, O2, O2, O2, O2, O2})) << endl;
        cout << decodeHam(ham2, u2 + RowMat({O2, O2, O2, O2, O2, O2, Z2, O2, O2, O2, O2, O2, O2})) << endl;
        cout << decodeHam(ham2, u2 + RowMat({O2, O2, O2, O2, O2, O2, O2, I2, O2, O2, O2, O2, O2})) << endl;
        cout << decodeHam(ham2, u2 + RowMat({O2, O2, O2, O2, O2, O2, Z2, I2, O2, O2, O2, O2, O2})) << endl;

        const Ham ham3(2, 5);
        cout << ham3.H << endl;
        cout << ham3.G << endl;
        const RingElem O3 = zero(ham3.R);
        const RingElem I3 = one(ham3.R);
        const RingElem Z3 = 2 * I3;
        const RingElem E3 = 3 * I3;
        const RingElem A3 = 4 * I3;
        const matrix w3 = NewDenseMat(RowMat({E3, A3, I3, O3}));
        const matrix u3 = encodeHam(ham3, w3);
        cout << u3 << endl;
        cout << decodeHam(ham3, u3) << endl;
        cout << decodeHam(ham3, u3 + RowMat({E3, O3, O3, O3, O3, O3})) << endl;
        cout << decodeHam(ham3, u3 + RowMat({O3, O3, O3, O3, O3, I3})) << endl;
        cout << decodeHam(ham3, u3 + RowMat({O3, O3, A3, O3, O3, O3})) << endl;
        cout << decodeHam(ham3, u3 + RowMat({O3, O3, Z3, O3, E3, O3})) << endl;
    }

    void testRM() {
        cout << "========================== RM ==========================" << endl << endl;

        cout << RM(4, 7) << endl;
    }

    void example(const int argc, const char *argv[]) {
        bool all = false, bch = false, ham = false, rm = false;

        if (argc <= 1) {
            all = true;
        } else {
            for (int i = 1; i < argc; i++) {
                const char *arg = argv[i];
                if (strcmp(arg, "bch") == 0) {
                    bch = true;
                } else if (strcmp(arg, "ham") == 0) {
                    ham = true;
                } else if (strcmp(arg, "rm") == 0) {
                    rm = true;
                }
            }
        }

        if (all || bch) {
            cout << endl << endl << endl;
            testBCH();
        }
        if (all || ham) {
            cout << endl << endl << endl;
            testHam();
        }
        if (all || rm) {
            cout << endl << endl << endl;
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
    } catch (const CoCoA::InterruptReceived &intr) {
        cerr << endl
             << "------------------------------" << endl
             << ">>>  CoCoALib interrupted  <<<" << endl
             << "------------------------------" << endl
             << "-->>  " << intr << "  <<--" << endl;
        return 2;
    } catch (const CoCoA::ErrorInfo &err) {
        cerr << "***ERROR***  UNCAUGHT CoCoA error";
        ANNOUNCE(cerr, err);
    } catch (const std::exception &exc) {
        cerr << "***ERROR***  UNCAUGHT std::exception: " << exc.what() << endl;
    } catch (...) {
        cerr << "***ERROR***  UNCAUGHT UNKNOWN EXCEPTION" << endl;
    }

    CoCoA::BuildInfo::PrintAll(cerr);
    return 1;
}
