#include <mex.h>
#include "graph.h"

/*
 * [e, cut] = BK(ivec, jvec, ijvec, nnodes)
 *
 * Runs BK on the graph whose edges are given by (ivec, jvec, ijvec, jivec) with nnodes nodes.
 * By convention, node nnodes + 1 is s and nnodes + 2 is t. ivec and jvec are 1-indexed.
 *
 * ijvec and jivec specify bidirectional edge weights. ijvec is the edge ivec -> jvec and
 * jivec is jvec -> ivec. NOTE: For s -> x or x -> t edges, jivec is ignored.
 *
 * e is a scalar output for the maximum flow. cut is a 1 x N logical vector indicating
 * the nodes assigned to the SOURCE.
 */

enum {
  iivec = 0,
  ijvec,
  iijvec,
  ijivec,
  innodes,
  nI
};

enum {
  oE = 0,
  oCut,
  nO
};

#define MAX(a, b) ((a) > (b)) ? (a) : (b)
#define mxGetLength(v) MAX(mxGetM(v), mxGetN(v))

void mexFunction(int nOut, mxArray *pOut[], int nIn, const mxArray *pIn[]) {
  // Setup and definitions
  if (nIn != nI) {
    mexErrMsgIdAndTxt("BK_mex:n_inputs:", "Must have %d inputs", nI);
  }
  if (nOut > nO) {
    mexErrMsgIdAndTxt("BK_mex:n_output:", "Must have at most %d outputs", nO);
  }
  /*
  if (mxGetNumberOfDimensions(pIn[iivec]) != 1 ||
      mxGetNumberOfDimensions(pIn[ijvec]) != 1 ||
      mxGetNumberOfDimensions(pIn[iijvec]) != 1) {
    mexErrMsgIdAndTxt("BK_mex:dim", "input vectors must be one dimensional");
  }
  */

  /*
  if ((mxGetLength(pIn[iivec]) != mxGetLength(pIn[ijvec])) ||
      (mxGetLength(pIn[ijvec]) != mxGetLength(pIn[iijvec])) ||
      (mxGetLength(pIn[iijvec]) != mxGetLength(pIn[ijivec]))) {
    mexErrMsgIdAndTxt("BK_mex:eq", "input vectors must have same size");
  }
  */

  int nEdges = mxGetLength(pIn[iivec]);
  int nNodes = mxGetScalar(pIn[innodes]);

  // we're now on 0-indexing.
  int sNode = nNodes;
  int tNode = nNodes + 1;

  double *iVec = mxGetPr(pIn[iivec]);
  double *jVec = mxGetPr(pIn[ijvec]);
  double *ijVec = mxGetPr(pIn[iijvec]);
  double *jiVec = mxGetPr(pIn[ijivec]);

  // Translate sparse entries to a BK graph
  typedef Graph<double, double, double> dGraph;
  dGraph *gp = new dGraph(nNodes, nEdges, mexErrMsgTxt);
  mexPrintf("nNodes = %d\n", nNodes);
  gp->add_node(nNodes);

  int i, j;
  double ij, ji;
  for (int n = 0; n < nEdges; n++) {
    // 0-based indexing
    i = iVec[n];
    j = jVec[n];
    ij = ijVec[n];
    ji = jiVec[n];
    if (i == sNode) {
      // In MATLAB, we added an entry (sNode, j). graph.h expects a call
      // (j, w, 0) (since we have only a flow from the source)
      gp->add_tweights(j, ij, 0);
      mexPrintf("[BK_mex] Added s -> %d edge; weight %g\n", j, ij);
    }
    else if (j == tNode) {
      gp->add_tweights(i, 0, ij);
      mexPrintf("[BK_mex] Added %d -> t edge; weight %g\n", i, ij);
    }
    else {
      if (i < 0 || i >= nNodes) {
        mexPrintf("i = %d\n", i);
        mexErrMsgIdAndTxt("BK_mex:ibounds", "i out of bounds");
      }
      if (j < 0 || j >= nNodes) {
        mexPrintf("j = %d\n", j);
        mexErrMsgIdAndTxt("BK_mex:jbounds", "j out of bounds");
      }
      gp->add_edge(i, j, ij, ji);
      mexPrintf("[BK_mex] Added %d -> %d edge; weights %g, %g\n", i, j, ij, ji);
    }
  }

  /****************************************************************
   * optimize!
   */
  mexPrintf("[BK_mex] before maxflow\n");
  double e = gp->maxflow();
  mexPrintf("[BK_mex] after maxflow\n");


  /****************************************************************
   * read results
   */
//    pOut[oE] = mxCreateDoubleScalar(e); // output the flow


  pOut[oE] = mxCreateDoubleScalar(e);
  mexPrintf("[BK_mex] before cut\n");
  if (nOut > 1) {
    pOut[oCut] = mxCreateLogicalMatrix(1, nNodes);
    mxLogical *pl = mxGetLogicals(pOut[oCut]);

    for (int n = 0; n < nNodes; n++) {
      pl[n] = gp->what_segment(n);
    }
  }
  mexPrintf("[BK_mex] after cut\n");
}
