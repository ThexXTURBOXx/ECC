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
        cout << "========================== RM ==========================" << endl << endl;

        cout << RM(2, 5) << endl;
    }

    void example(const int argc, const char *argv[]) {
        bool all = false, bch = false, ham = false, rm = false;

        if (argc <= 1) {
            all = true;
        } else {
            for (int i = 1; i < argc; ++i) {
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
