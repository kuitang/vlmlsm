#include "mexcpp.h"
#include "BetheApprox.h"
using namespace mexcpp;

enum {
  oA,
  oB,
  nO,
};

enum {
  iTheta,
  iW,
  iEpsilon,
  iMethod, /* 0 -> BBP ; 1 -> MK */
  nI,
};

void mexFunction(int nOut, mxArray *pOut[], int nIn, const mxArray *pIn[]) {
  const char *usage = "Usage: [A, B] = BBP_MK_mex.cpp(theta, W, epsilon, method) where method = 0 -> BBP; 1 -> MK";
  if (nIn != nI || nOut != nO) {
    mexErrMsgIdAndTxt("BBP_MK_mex:args", usage);
  }

  Mat<double>  theta           (pIn[0]);
  cscMatrix    W = extractCSC  (pIn[1]);
  double       epsilon = scalar<double>(pIn[2]);
  double       method  = scalar<double>(pIn[3]);

  int nNodes = theta.length;

  Mat<double> Amat(nNodes, 1);
  double *A = Amat.re;

  Mat<double> Bmat(nNodes, 1);
  double *B = Bmat.re;

  // Defaults as in BetheApprox_mex.cpp; TODO: Change
  double bbThresh = 1e-4;
  int maxIter = 10000;

  double alpha[W.nzMax];

  if (method == 0.0) {
    propogateBetheBound(nNodes, theta, W, bbThresh, maxIter, A, B, alpha);
  } else {
    mooijBound(nNodes, theta, W, bbThresh, maxIter, A, B, alpha);
  }

  pOut[oA] = Amat;
  pOut[oB] = Bmat;
}

