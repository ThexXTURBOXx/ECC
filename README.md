# ECC

Error correction code implementations in CoCoALib

## Setup (tested using CygWin and Ubuntu)

1. Download the [CoCoALib](https://cocoa.dima.unige.it/cocoa/cocoalib/) source code (tested with `0.99815`, snapshots
   available [here](https://cocoa.dima.unige.it/cocoa/cocoalib/snapshot.shtml))
2. Extract the archive into `cocoalib/`
3. Maybe, you need to apply the patch mentioned in [Fix LaTeX errors](#fix-latex-errors)
4. Install [GMP](https://gmplib.org/) using your favorite package manager or compile it yourself
5. Install a valid LaTeX distribution like [TeXLive](https://www.tug.org/texlive/) or [MikTeX](https://miktex.org/)
6. Install the latest version of [txt2tags 2.7-dev](https://txt2tags.org/) along with a compatible [Python](https://www.python.org/) version
7. Install the CoCoALib library and header files (depending on your setup, you might need to run this as an admin): `./make_cocoalib.sh`
8. Setup CMake using the `cmake` command and appropriate arguments
9. Build the resulting executables using `cmake --build`

## Fix LaTeX errors

**Update: CoCoALib >0.99817 should not have this issue anymore!**

In my tests, at least on Windows, there may be LaTeX compilation errors in the CoCoALib documentation that look similar to this:
```
Doing LaTeX-->PDF
[[CoCoALib:LaTeX-->PDF]]  Pass 1 of 3
!!!LaTeX error!!! -- see file tex/CoCoALib.log
make[1]: *** [Makefile:40: CoCoALib.pdf] Error 1
```
In order to fix those, open the file `cocoalib/doc/aux-txt2tags/cocoalib-doc.sty` and uncomment the following line:
```tex
%\newcommand{\htmladdnormallink}[2]{\href{#2}{#1}}
```
