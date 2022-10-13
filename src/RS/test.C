#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string description =
        "This file provides simple examples for the Reed-Solomon Code.\n";
//----------------------------------------------------------------------

namespace CoCoA {

    void example() {
        cout << "Hello world" << endl;

        auto R = NewPolyRing(RingZZ(), symbols("x"));
        auto f = RingElem(R, "x^2-5");
        cout << IsLRSDegenerate(f) << endl;
    }

// DO NOT EDIT LINES BELOW HERE

    void program() {
        GlobalManager CoCoAFoundations;
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
