#ifndef ECC_UTILS_H
#define ECC_UTILS_H

#include "CoCoA/library.H"

using namespace std;

/**
 * Namespace of CoCoALib.
 */
namespace CoCoA {

  /**
   * Returns the `i`-th element of a vector, or the default value if `i` is out of bounds.
   * @tparam T The generic type of the vector
   * @param vec The vector
   * @param i The index
   * @param defaultVal The default value
   * @return The `i`-th element of the vector, or the default value if `i` is out of bounds
   */
  template<class T>
  T getOr(const vector<T> &vec, size_t i, T defaultVal);

  /**
   * Parses a character to a number.
   * @param c The character to parse
   * @return The parsed number
   */
  static int parseNum(char c);

  /**
   * Converts a number to a character.
   * @param i The number to convert
   * @return The converted character
   */
  char toChar(long i);

  /**
   * Converts a string to a polynomial.
   * @param str The string to convert
   * @param x The variable with respect to which the polynomial is defined
   * @return The converted polynomial
   */
  RingElem toPolynomial(const string &str, ConstRefRingElem x);

  /**
   * Converts a vector to a polynomial.
   * @param mat The vector to convert
   * @param x The variable with respect to which the polynomial is defined
   * @return The converted polynomial
   */
  RingElem toPolynomial(const matrix &mat, ConstRefRingElem x);

  /**
   * Converts a string to a vector.
   * @param str The string to convert
   * @param R The ring over which the vector should be defined
   * @return The converted vector
   */
  matrix toMatrix(const string &str, const ring &R);

  /**
   * Converts a polynomial to a vector.
   * @param p The polynomial to convert
   * @param n The size of the resulting vector
   * @param x The variable with respect to which the polynomial is defined
   * @return The converted matrix
   */
  matrix toMatrix(ConstRefRingElem p, long n, ConstRefRingElem x);

  /**
   * Converts a polynomial to a string.
   * @param p The polynomial to convert
   * @param n The size of the resulting string
   * @param x The variable with respect to which the polynomial is defined
   * @return The converted string
   */
  string toString(ConstRefRingElem p, long n, ConstRefRingElem x);

  /**
   * Converts a vector to a string.
   * @param m The vector to convert
   * @return The converted string
   */
  string toString(const matrix &m);

  /**
   * Calculates the set of all tuples of a given size of a given set.
   * @tparam T The generic type of the vector
   * @param set The set of which the tuples should be calculated
   * @param tupleSize The size of the tuples
   * @return The set of all tuples of the given size of the given set
   */
  template<class T>
  vector<vector<T>> tuples(const vector<T> &set, long tupleSize);

  /**
   * Calculates the set of all subsets of a given size of a given set.
   * @tparam T The generic type of the vector
   * @param set The set of which the subsets should be calculated
   * @param setSize The size of the subsets
   * @param ret Output parameter for the set of all subsets of the given size of the given set
   * @see This is an improved version of `examples/ex-MVT-Simplicial.C` in CoCoALib
   */
  template<class T>
  void subsets(const vector<T> &set, long setSize, vector<vector<T>> &ret);

  /**
   * Calculates the standard basis vector of the given size with the given element at the given index.
   * @param R The underlying ring
   * @param i The index of the only non-zero element
   * @param b The value of the only non-zero element
   * @param n The size of the resulting vector
   * @return The standard basis vector of the given size with the given element at the given index
   */
  matrix e(const ring &R, long i, const RingElem &b, long n);

  /**
   * Calculates the cyclic right shift of the given vector by the given amount.
   * @tparam T The generic type of the vector
   * @param vec The vector to shift
   * @param s The amount to shift (negative values shift left)
   * @return The cyclic right shift of the given vector by the given amount
   */
  template<class T>
  vector<T> cycShift(const vector<T> &vec, long s);

  /**
   * Calculates a reverse circulant matrix with the given first row.
   * @param firstRow The first row of the matrix
   * @return The reverse circulant matrix with the given first row
   */
  matrix revCirculantMatrix(const vector<RingElem> &firstRow);

  /**
   * Calculates the weight of the given vector, i.e., the number of non-zero elements.
   * @param v The vector
   * @return The weight of the given vector
   */
  long wt(const vector<long> &v);

  /**
   * Calculates the weight of the given matrix, i.e., the number of non-zero elements.
   * @param m The matrix
   * @return The weight of the given matrix
   */
  long wt(const ConstMatrixView &m);

  /**
   * Does a Chien search for roots of the given polynomial over the given field.
   * @param f The polynomial to search for roots
   * @param a A primitive element of the field
   * @param qn The size of the field
   * @param x The variable of the polynomial
   * @return The powers of the primitive element that are roots of the polynomial
   * @see Chien, R. (1964). Cyclic decoding procedures for Bose-Chaudhuri-Hocquenghem codes
   * @see {@link ChienSearchSingleRoot}
   */
  vector<long> ChienSearch(ConstRefRingElem f, ConstRefRingElem a, long qn, ConstRefRingElem x);

  /**
   * Does a Chien search for roots of the given polynomial over the given field, but stops after finding the first root.
   * @param f The polynomial to search for roots
   * @param a A primitive element of the field
   * @param qn The size of the field
   * @param x The variable of the polynomial
   * @return A single power of the primitive element that is a root of the polynomial
   * @see Chien, R. (1964). Cyclic decoding procedures for Bose-Chaudhuri-Hocquenghem codes
   * @see {@link ChienSearch}
   */
  long ChienSearchSingleRoot(ConstRefRingElem f, ConstRefRingElem a, long qn, ConstRefRingElem x);

  /**
   * Returns a polynomial in the given set of polynomials that is univariate in the given variable.
   * @param G The set of polynomials
   * @param IndetIndex The index of the variable
   * @param fallback The polynomial to return if no polynomial in the set is univariate in the given variable
   * @return A polynomial in the given set of polynomials that is univariate in the given variable
   */
  RingElem getUniPoly(const vector<RingElem> &G, long IndetIndex, ConstRefRingElem fallback);

  /**
   * Returns the binomial coefficient of the given two numbers as a `long`.
   * @param n The first number (top)
   * @param k The second number (bottom)
   * @return The binomial coefficient of the given two numbers as a `long`
   */
  long binom(long n, long k);

  /**
   * Returns the "length" of the given matrix, i.e., the maximum of the number of rows and columns.
   * @param m The matrix
   * @return The length of the given matrix
   */
  long GetLength(const matrix &m);

  // As of CoCoALib 0.99815, this has been incorporated into CoCoALib itself
  /*
   * Checks if f is a primitive polynomial in ZZ/(p).
   * @param f The polynomial to check
   * @return Whether f is primitive
   */
  /*bool IsPrimitivePoly(ConstRefRingElem f);*/

  /**
   * Returns a primitive polynomial of the given degree over the given field through brute force/randomised search.
   * @param Px The polynomial ring
   * @param n The desired degree of the polynomial
   * @param IndetIndex The index of the variable
   * @return A primitive polynomial of the given degree over the given field
   */
  RingElem BruteForcePrimPoly(const ring &Px, long n, long IndetIndex);

}

#endif // ECC_UTILS_H
