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
const string description = "This file provides simple examples for this software.\n";
//----------------------------------------------------------------------

namespace CoCoA {

  // Forward decl -- start
  void exampleECC(const int argc, const char *argv[]);

  void exampleFuzzy(const int argc, const char *argv[]);
  // Forward decl -- end

  void program(const int argc, const char *argv[]) {
    GlobalManager CoCoAFoundations(UseNonNegResidues);
    SignalWatcher MonitorInterrupt(SIGINT); // you must also call CheckForInterrupt every so often

    cout << description << endl;
    cout << boolalpha; // so that bools print out as true/false

    if (argc < 2) {
      CoCoA_THROW_ERROR("Too few arguments!", "test");
    }

    if (strcmp(argv[1], "ecc") == 0) {
      exampleECC(argc - 2, argv + 2);
    } else if (strcmp(argv[1], "fuzzy") == 0) {
      exampleFuzzy(argc - 2, argv + 2);
    } else {
      cout << "" << endl;
    }
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
