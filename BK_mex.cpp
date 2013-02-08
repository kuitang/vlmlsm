#include <mex.h>
#include "graph.h"

/*
 * [e, cut] = BK(iVec, jVec, ijVec, jiVec, nNodes, sNode, tNode)
 *
 * Runs BK on the graph whose edges are given by (iVec, jVec, ijVec, jiVec)
 * with nNodes nodes.
 *
 * ijVec and jivec specify bidirectional edge weights. ijVec is the edge
 * ivec -> jvec and jivec is jvec -> ivec.
 *
 * NOTE: For s -> x or x -> t edges, jivec is ignored.
 *
 * e is a scalar output for the maximum flow. cut is a 1 x (N+1) logical vector indicating
 * the nodes assigned to the SINK. Note the extra dimension to be agnostic with 0/1 indexing.
 */

enum {
  iiVec = 0,
  ijVec,
  iijVec,
  ijiVec,
  inNodes,
  isNode,
  itNode,
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
  if (mxGetNumberOfDimensions(pIn[iiVec]) != 1 ||
      mxGetNumberOfDimensions(pIn[ijVec]) != 1 ||
      mxGetNumberOfDimensions(pIn[iijVec]) != 1) {
    mexErrMsgIdAndTxt("BK_mex:dim", "input vectors must be one dimensional");
  }
  */

  /*
  if ((mxGetLength(pIn[iiVec]) != mxGetLength(pIn[ijVec])) ||
      (mxGetLength(pIn[ijVec]) != mxGetLength(pIn[iijVec])) ||
      (mxGetLength(pIn[iijVec]) != mxGetLength(pIn[ijiVec]))) {
    mexErrMsgIdAndTxt("BK_mex:eq", "input vectors must have same size");
  }
  */

  int nEdges = mxGetLength(pIn[iiVec]);
  int nNodes = mxGetScalar(pIn[inNodes]);
  int sNode  = mxGetScalar(pIn[isNode]);
  int tNode  = mxGetScalar(pIn[itNode]);

  double *iVec = mxGetPr(pIn[iiVec]);
  double *jVec = mxGetPr(pIn[ijVec]);
  double *ijVec = mxGetPr(pIn[iijVec]);
  double *jiVec = mxGetPr(pIn[ijiVec]);

  // Translate sparse entries to a BK graph
  typedef Graph<double, double, double> dGraph;
  dGraph *gp = new dGraph(nNodes, nEdges, mexErrMsgTxt);
  //mexPrintf("nNodes = %d\n", nNodes);
  // IMPORTANT! Enables one-indexing transparently from MATLAB. (We can
  // effectively ignore node 0 if we want.)
  gp->add_node(nNodes + 1);

  int i, j;
  double ij, ji;
  for (int n = 0; n < nEdges; n++) {
    // 0-based indexing
    i = iVec[n];
    j = jVec[n];
    ij = ijVec[n];
    ji = jiVec[n];
    // Allow for empty edges.
    if ((i == 0  && j == 0) || (ij == 0.0 && ji == 0.0)) {
      continue;
    }
    if (i == sNode) {
      // In MATLAB, we added an entry (sNode, j). graph.h expects a call
      // (j, w, 0) (since we have only a flow from the source)
      gp->add_tweights(j, ij, 0);
      //mexPrintf("[BK_mex] Added s -> %d edge; weight %g\n", j, ij);
    }
    else if (j == tNode) {
      gp->add_tweights(i, 0, ij);
      //mexPrintf("[BK_mex] Added %d -> t edge; weight %g\n", i, ij);
    }
    else {
      if (i < 0 || i > nNodes) {
        mexPrintf("i, j = %d, %d\n", i, j);
        mexErrMsgIdAndTxt("BK_mex:ibounds", "i out of bounds");
      }
      if (j < 0 || j > nNodes) {
        mexPrintf("i, j = %d, %d\n", i, j);
        mexErrMsgIdAndTxt("BK_mex:jbounds", "j out of bounds");
      }
      gp->add_edge(i, j, ij, ji);
      //mexPrintf("[BK_mex] Added %d -> %d edge; weights %g, %g\n", i, j, ij, ji);
    }
  }

  /****************************************************************
   * optimize!
   */
  //mexPrintf("[BK_mex] before maxflow\n");
  double e = gp->maxflow();
  //mexPrintf("[BK_mex] after maxflow\n");


  /****************************************************************
   * read results
   */
//    pOut[oE] = mxCreateDoubleScalar(e); // output the flow


  pOut[oE] = mxCreateDoubleScalar(e);
  //mexPrintf("[BK_mex] before cut\n");
  if (nOut > 1) {
    // IMPORTANT: Support both 0 and 1 indexing. (Depending on convention,
    // ignore either the head or tail.)
    pOut[oCut] = mxCreateLogicalMatrix(1, nNodes + 1);
    mxLogical *pl = mxGetLogicals(pOut[oCut]);

    for (int n = 0; n <= nNodes; n++) {
      pl[n] = gp->what_segment(n) == dGraph::SINK;
    }
  }
  //mexPrintf("[BK_mex] after cut\n");
}
