#ifndef ECC_UTILS_H
#define ECC_UTILS_H

#include "CoCoA/library.H"

using namespace std;

namespace CoCoA {

  template<class T>
  T getOr(const vector<T> &vec, size_t i, T defaultVal);

  static int parseNum(char c);

  char toChar(long i);

  RingElem toPolynomial(const string &str, ConstRefRingElem x);

  RingElem toPolynomial(const matrix &mat, ConstRefRingElem x);

  matrix toMatrix(const string &str, const ring &R);

  matrix toMatrix(ConstRefRingElem p, long n, ConstRefRingElem x);

  string toString(ConstRefRingElem p, long n, ConstRefRingElem x);

  string toString(const matrix &m);

  template<class T>
  vector<vector<T>> tuples(const vector<T> &set, long tupleSize);

  template<class T>
  void subsets(const vector<T> &set, long setSize, vector<vector<T>> &ret);

  matrix e(const ring &R, long i, const RingElem &b, long n);

  template<class T>
  vector<T> cycShift(const vector<T> &vec, long s);

  matrix revCirculantMatrix(const vector<RingElem> &firstRow);

  long wt(const vector<long> &v);

  long wt(const ConstMatrixView &m);

  // f: Polynomial to find roots of
  // a: Primitive element in field
  // qn: Amount of elements in field
  // x: Variable x
  vector<long> ChienSearch(ConstRefRingElem f, ConstRefRingElem a, long qn, ConstRefRingElem x);

  long ChienSearchSingleRoot(ConstRefRingElem f, ConstRefRingElem a, long qn, ConstRefRingElem x);

  RingElem getUniPoly(const vector<RingElem> &G, long IndetIndex, ConstRefRingElem fallback);

  long binom(long n, long k);

  long GetLength(const matrix &m);

  // As of CoCoALib 0.99815, this has been incorporated into CoCoALib itself
  /*bool IsPrimitivePoly(ConstRefRingElem f);*/

  RingElem BruteForcePrimPoly(const ring &Px, long n, long IndetIndex);

}

#endif // ECC_UTILS_H
