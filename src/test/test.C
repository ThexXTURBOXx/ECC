#include <cstring>
#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string description = "This file provides simple examples for this software.\n";
//----------------------------------------------------------------------

namespace CoCoA {
  namespace ECC {
    // Forward decl -- start
    void exampleECC(int argc, const char* argv[]);

    void exampleFuzzy(int argc, const char* argv[]);
    // Forward decl -- end

    /**
     * Delegates the program to the appropriate function.
     * @param argc The number of arguments
     * @param argv The arguments
     */
    void program(const int argc, const char* argv[]) {
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
}

/**
 * Entry point of the test program.
 * @param argc The number of arguments
 * @param argv The arguments
 * @return The exit code
 */
int main(const int argc, const char* argv[]) {
  try {
    CoCoA::ECC::program(argc, argv);
    return 0;
  } catch (const CoCoA::InterruptReceived& intr) {
    cerr << endl
      << "------------------------------" << endl
      << ">>>  CoCoALib interrupted  <<<" << endl
      << "------------------------------" << endl
      << "-->>  " << intr << "  <<--" << endl;
    return 2;
  }
  catch (const CoCoA::ErrorInfo& err) {
    cerr << "***ERROR***  UNCAUGHT CoCoA error";
    ANNOUNCE(cerr, err);
  }
  catch (const std::exception& exc) {
    cerr << "***ERROR***  UNCAUGHT std::exception: " << exc.what() << endl;
  }
  catch (...) {
    cerr << "***ERROR***  UNCAUGHT UNKNOWN EXCEPTION" << endl;
  }

  CoCoA::BuildInfo::PrintAll(cerr);
  return 1;
}
