# ECC

Error correction code implementations in CoCoALib

## Setup (tested using CygWin and Ubuntu)

1. Download the [CoCoALib](https://cocoa.dima.unige.it/cocoa/cocoalib/) source code (tested with `0.99810`, snapshots
   available [here](https://cocoa.dima.unige.it/cocoa/cocoalib/snapshot.shtml))
2. Extract the archive into `cocoalib/`
3. Install [GMP](https://gmplib.org/) using your favorite package manager or compile it yourself
4. Install the CoCoALib library and header files: `./make_cocoalib.sh`
5. Setup CMake using the `cmake` command and appropriate arguments
6. Build the resulting executables using `cmake --build`
