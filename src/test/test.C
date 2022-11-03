#include "CoCoA/library.H"
#include "util/utils.H"
#include "ecc/bch.H"

using namespace std;

//----------------------------------------------------------------------
const string description = "This file provides simple examples of BCH-Codes.\n";
//----------------------------------------------------------------------

namespace CoCoA {

    void example() {
        cout << constructGenPoly(1, 7, 2, "alpha^4+alpha+1") << endl;
        cout << constructGenPoly(1, 5, 3, "alpha^2+alpha+2") << endl;

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
        RingElem sent = encode(bch, p, x);
        cout << toString(sent, n, x) << endl;

        // Decoding

        RingElem recv1 = sent;
        RingElem recv2 = sent + power(x, 2);
        RingElem recv3 = sent + power(x, 13) + power(x, 5);
        RingElem recv4 = sent + power(x, 8) + power(x, 2) + 1;

        cout << "---" << endl;

        RingElem dec1 = decode(bch, recv1, x);
        RingElem dec2 = decode(bch, recv2, x);
        RingElem dec3 = decode(bch, recv3, x);
        RingElem dec4 = decode(bch, recv4, x);

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
        sent = encode(bch2, p, x);
        cout << toString(sent, n, x) << endl;

        // Decoding

        recv1 = sent;
        recv2 = sent + 2 * power(x, 2);
        recv3 = sent + power(x, 5) + 1;

        cout << "---" << endl;

        dec1 = decode(bch2, recv1, x);
        dec2 = decode(bch2, recv2, x);
        dec3 = decode(bch2, recv3, x);

        cout << toString(dec1, n, x) << endl;
        cout << toString(dec2, n, x) << endl;
        cout << toString(dec3, n, x) << endl;
    }

    // DO NOT EDIT LINES BELOW HERE

    void program() {
        GlobalManager CoCoAFoundations(UseNonNegResidues);
        SignalWatcher MonitorInterrupt(SIGINT); // you must also call CheckForInterrupt every so often

        cout << description << endl;
        cout << boolalpha; // so that bools print out as true/false

        example();
    }

}

int main() {
    try {
        CoCoA::program();
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
