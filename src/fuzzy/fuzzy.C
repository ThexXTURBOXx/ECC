#include "CoCoA/library.H"
#include "fuzzy/fuzzy.H"

using namespace std;

namespace CoCoA {

  matrix FuzzyExtractor::generateHelperData(const matrix &w) {
    if (helperDataSet) CoCoA_THROW_ERROR("Helper data already set", __func__);
    const long respCols = GetLength(w);

    matrix x = NewDenseMat(RingOf(w), 1, respCols);
    for (int i = 0; i < respCols; ++i)
      SetEntry(x, 0, i, RandomBool());

    // SS
    matrix k = NewDenseMat(RingOf(w), 1, respCols - parityBits);
    for (int i = 0; i < respCols - parityBits; ++i)
      SetEntry(k, 0, i, RandomBool());
    const matrix r = encode(k);
    const matrix s = r + w;

    hd = HelperData(s, x);

    // Ext
    // Strong Extract to receive R
    return strongExtract(w + x);
  }

  matrix FuzzyExtractor::extract(const matrix &wd) {
    vector<long> cols(messageBits);
    iota(begin(cols), end(cols), 0);

    // Rec
    const matrix rd = wd + hd.s;
    const matrix k = NewDenseMat(submat(decode(rd), {0}, cols));
    const matrix r = encode(k);
    const matrix w = hd.s + r;

    // Ext
    // Strong Extract to receive R
    return strongExtract(hd.x + w);
  }

  matrix FuzzyExtractor::strongExtract(const matrix &w) {
    long len = GetLength(w);
    uint8_t toHash[len];
    for (long i = 0; i < len; ++i)
      toHash[i] = IsOne(w(0, i));

    SHA256 sha256;
    sha256.update(ToString(w)); // As long as ToString does not change, let's use that
    const uint8_t *hash = sha256.digest();

    matrix ext = NewDenseMat(RingOf(w), 1, 32 * 8);
    for (int i = 0; i < 32; ++i) {
      const uint8_t num = hash[i];
      for (int j = 7; j >= 0; --j)
        SetEntry(ext, 0, i * 8 + j, (num >> j) & 1);
    }

    delete hash;
    return ext;
  }

}
