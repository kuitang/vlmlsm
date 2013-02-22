#include <algorithm>
#include "mexcpp.h"
#include "matrix.h"
#include "MinSum.h"
#include "BetheApprox.h"

using namespace mexcpp;

// Test the marginalize, entropy and binaryEntropy functions
void mexFunction(int nOut, mxArray *pOut[], int nIn, const mxArray *pIn[]) {
  const char *usage = "Usage: [margs, margEnts, binEnts] = AuxTests_mex(vals, alphas) -- where vals is N-vector of values in [0, 1] to test. Alphas is N^2 long.";
  if (nIn != 2 || nOut != 3) {
    mexErrMsgIdAndTxt("AuxTests_mex:args", usage);
  }

  Mat<double> vals(pIn[0]);
  mxAssert(vals.N == 1 || vals.M == 1, "vals must be vector");
  size_t N = vals.N * vals.M;

  Mat<double> alphas(pIn[1]);
  mxAssert(alphas.length == N*N, "alphas must have length N*N (where N = length(vector))");


  // test binary entropy
  Mat<double> binEnts(N, 1);
  std::transform(vals.re, vals.re + N, binEnts.re, binaryEntropy);
  pOut[2] = binEnts;

  Mat<double> margs(4, N*N);
  Mat<double> margEnts(N*N, 1);

  size_t ind = 0;
  for (size_t i = 0; i < N; i++) {
    for (size_t j = 0; j < N; j++) {
      marginalize(alphas[ind], vals[i], vals[j], &margs[ind * 4]);
      margEnts[ind] = entropy<4>(&margs[ind * 4]);
      ind++;
    }
  }

  pOut[0] = margs;
  pOut[1] = margEnts;
}

