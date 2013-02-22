#include "mexcpp.h"
#include "matrix.h"
#include "BetheApprox.h"
using namespace mexcpp;

enum {
  oLogZ,           /* 0 */
  oOneMarginals,   /* 1 */
  oTwoMarginals,   /* 2 */
  oMisc            /* 3 */
};

// Compressed sparse column layout. (MATLAB's format)
//
// jc[j] is the linear index of the first nonzero element in column j and
// jc[j+1]-1 is the linear index of the last nonzero element in column j.
//
// ir[i] is the ROW of the node at linear index i.
cscMatrix extractCSC(const mxArray *M) {
  if (!mxIsDouble(M) || !mxIsSparse(M)) {
    mexErrMsgIdAndTxt("extractCSC:M", "matrix M must be double and sparse");
  }

  return { .N = mxGetN(M),
           .M = mxGetM(M),
           .nzMax = mxGetNzmax(M),
           .pr = mxGetPr(M),
           .ir = mxGetIr(M),
           .jc = mxGetJc(M) };
}

void mexFunction(int nOut, mxArray *pOut[], int nIn, const mxArray *pIn[]) {
  const char *usage = "Usage: [logZ, oneMarginals, twoMarginals, misc] = BetheApprox_mex(theta, W, epsilon, opts)";
  if (nIn != 4 || nOut != 4) {
    mexErrMsgIdAndTxt("BetheApprox_mex:args", usage);
  }

  Mat<double> theta                      (pIn[0]);
  cscMatrix   W = extractCSC             (pIn[1]);
  double      epsilon = scalar<double>(pIn[2]);
  StructMat   opts                       (pIn[3]);

  Entry       opt(opts[0]);

  int nNodes = theta.length;
  mwSize nEdges = W.nzMax / 2; // symmetric mat

  // DO COMPUTATION
  double A[nNodes];
  double B[nNodes];
  double alpha[W.nzMax];

  // default thresh
  double bbThresh = 0.002;
  int maxIter = 50000; // something really huge

  propogateBetheBound(nNodes, theta, W, bbThresh, maxIter, A, B, alpha);
  double intervalSz = getIntervalSz(nNodes, A, B, W, alpha, epsilon);

  // DEBUGGING: COMPARE THE BOUND PROPAGATION (TODO)
  auto misc = StructMat(1, 1, {"A", "B", "intervalSz", "Vm", "elMat", "e"});

  Mat<double> Amat(nNodes, 1);
  std::copy(A, A + nNodes, Amat.re);
  misc[0].set("A", Amat);

  auto Bmat = Mat<double>(nNodes, 1);
  std::copy(B, B + nNodes, Bmat.re);
  misc[0].set("B", Bmat);
  misc[0].setS("intervalSz", intervalSz);

  pOut[oMisc] = misc;

  MinSum m(nNodes, nEdges, mexErrMsgTxt);
  makeBetheMinSum(nNodes, theta, W, A, B, alpha, intervalSz, m);

  // Output the potentials
  CellMat<Mat<double>> Vm(m.potentials.size(), 1);
  for (int i = 0; i < m.potentials.size(); i++) {
    Potential &p = *m.potentials[i];
    Mat<double> potMat(p.nLo, p.nHi);
    std::copy(&p(0,0), &p(0,0) + p.nLo * p.nHi, potMat.re);
    Vm.set(i, static_cast<mxArray *>(potMat));
  }
  misc[0].set("Vm", Vm);

  m.validate();

  mexPrintf("MINIMIZING -- READ IT NOW!\n");

  std::vector<int> x;
  double energy[3];
  double maxFlow;
  m.minimize(x, energy, maxFlow);

  misc[0].set("e", Mat<double>(1, 3, energy));

  // Output the edge list
  Mat<double> elMat(m.debugEdges.size(), 4);
  for (int i = 0; i < m.debugEdges.size(); i++) {
    const auto &de = m.debugEdges[i];
    elMat(i,0) = de.src;
    elMat(i,1) = de.dst;
    elMat(i,2) = de.fw;
    elMat(i,3) = de.rw;
  }
  misc[0].set("elMat", elMat);

  pOut[oLogZ] = scalar<double>(-energy[0]);

  Mat<double> oneMarginals(nNodes, 1);
  pOut[oOneMarginals] = oneMarginals;
  for (mwIndex n = 0; n < nNodes; n++) {
    if (x[n] == m.nodes[n]->nStates - 1) { // if we're at the end
      oneMarginals[n] = 1 - B[n];
    } else {
      oneMarginals[n] = A[n] + x[n] * intervalSz;
    }
  }

  mexPrintf("%s:%d -- oneMarginals success\n", __FILE__, __LINE__);

  // Recover marginals
  mwSize dims[] = {2, 2, nEdges};
  // TODO: Write wrapper into mexcpp
  pOut[oTwoMarginals] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
  double *pTwoMarg = mxGetPr(pOut[oTwoMarginals]);

  mwIndex i, off;
  mwIndex mIdx = 0;
  for (mwIndex j = 0; j < nNodes; j++) {
    for (mwIndex wIdx = W.jc[j]; wIdx < W.jc[j+1]; wIdx++) {
      i = W.ir[wIdx];

      // Upper triangular only
      if (j > i) {
        off = 4 * mIdx; // we vary along the 3rd dimension.
        marginalize(alpha[wIdx], oneMarginals[i], oneMarginals[j], pTwoMarg + off);
        mIdx++;
      }
    }
  }
  mexPrintf("%s:%d -- twoMarginals success\n", __FILE__, __LINE__);
}
