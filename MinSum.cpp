#include <cmath>
#include <cstdio>
#include <numeric>
#include "MinSum.h"

// TODO: errMsg function

const double TEN_EPS = 2.2204e-15;
const double INF = 1e100;

// Checks that potentials are submodular and that the neighbors matrix
// is symmetric and potentials are transposed when needed.

bool MinSum::isSubModular(Potential *p) {
  Potential &q = *p;
  char errMsg[1000];
  for (size_t lo = 0; lo < q.nLo - 1; lo++) {
    for (size_t hi = 0; hi < q.nHi - 1; hi++) {
      if (q(lo,hi) + q(lo+1, hi+1) > q(lo, hi+1) + q(lo+1, hi)) {
        snprintf(errMsg, 1000, "isSubModular: GEQ Potential %p entry (%lu,%lu) failed (%g + %g > %g + %g)", p, lo, hi, q(lo,hi), q(lo+1,hi+1), q(lo,hi+1), q(lo+1,hi));
        errFunc(errMsg);
        return false;
      }
    }
  }

  return true;
}

MinSum::FailCode MinSum::validate() {
  for (Potential *p : potentials) {
    if (!isSubModular(p)) { return FailCode::NOT_TRANSPOSED_POTENTIAL;  }
  }

  for (size_t src = 0; src < nNodes; src++) {
    for (const Edge &e : neighbors[src]) {
      size_t dst = e.dst;
      if ( (src > dst && !e.srcHigher) ||
           (src < dst && e.srcHigher) ) { return FailCode::NOT_TRANSPOSED_POTENTIAL; }
      if (src == dst) { return FailCode::SELF_LOOP; }
    }
  }

  return FailCode::SUCCESS;
}

void MinSum::minimize(std::vector<int> &x, double *energy, double &maxFlow) {
  // Compute offsets to convert (r, k) indices to linear
  std::vector<int> offsets;
  offsets.push_back(0);
  for (int r = 1; r < nNodes; r++) {
    offsets.push_back(offsets[r - 1] + nodes[r - 1]->nStates - 1);
  }
  assert(offsets.size() == nNodes);

  int nBKNodes = offsets.back() + nodes[nNodes - 1]->nStates - 1;

  for (int r = 0; r < nNodes; r++) {
    //mexPrintf("%s:%d -- Node %d had %d represented states. Offset is %d\n",
    //          __FILE__, __LINE__, r, nodes[r]->nStates - 1, offsets[r]);
  }

  // upper bound on number of neighbors:
  // one edge per node for s/t (but those don't count)
  // one edge per state for the zero/infinity (captured by offset)
  // variable neighbors for neighbors in underlying graph
  int nBKneighbors = offsets.back() + nodes.back()->nStates;

  for (int r = 0; r < nNodes; r++) {
    for (const Edge &e : neighbors[r]) {
      int rr = e.dst;
      int rStates  = nodes[r]->nStates - 2;
      int rrStates = nodes[rr]->nStates - 2;
      nBKneighbors += rStates * rrStates;
    }
  }

  // Make our graph
  dGraph *gp = new dGraph(nBKNodes, nBKneighbors, errFunc);
  //mexPrintf("%s:%d -- nBKNodes = %d\n", __FILE__, __LINE__, nBKNodes);

  // Add nodes
  gp->add_node(nBKNodes);

  // Zero/Infinity weights
  int loNode, hiNode;
  for (int r = 0; r < nodes.size(); r++) {
    for (int k = 0; k < nodes[r]->nStates - 2; k++) {
      loNode = offsets[r] + k;
      hiNode = offsets[r] + k + 1;

      gp->add_edge(loNode, hiNode, 0.0, INF);
      debugEdges.emplace_back(loNode + 1, hiNode + 1, 0.0, INF);

      //mexPrintf("%s:%d -- EDGE %d -> %d FW = %g RW = %g\n",
      //          __FILE__, __LINE__, loNode, hiNode, 0.0, INF);
    }
  }

  // Unary source/sink weights
  double qrk;
  for (int r = 0; r < nNodes; r++) {
    Node &nr = *nodes[r];

    for (int k = 0; k < nodes[r]->nStates - 1; k++) {
      qrk = 0;

      for (const Edge &e : neighbors[r]) {
        int end = nodes[e.dst]->nStates - 1;
        qrk += e.w * (e.p(k,0) + e.p(k,end) - e.p(k+1,0) - e.p(k+1,end));
        //mexPrintf("%s:%d qrk for r = %d, k = %d, rr = %d is now %g\n",
        //          __FILE__, __LINE__, r, k, e.dst, qrk);
        //mexPrintf("%s:%d Components: %g %g %g %g\n", __FILE__, __LINE__, e.p(k,1), e.p(k,end), e.p(k+1,1), e.p(k+1,end));
      }

      qrk = qrk / 2.0;
      qrk = qrk + nr(k) - nr(k+1);
      //mexPrintf("%s:%d qrk for r = %d, k = %d, FINAL is now %g\n", __FILE__, __LINE__, r, k, qrk);
      //mexPrintf("%s:%d Unary term was %g\n", __FILE__, __LINE__, nr(k) - nr(k+1));
      int node = offsets[r] + k;

      if (qrk > 0) {
        gp->add_tweights(node, qrk, 0);
        debugEdges.emplace_back(-1, node + 1, qrk, 0);
        //mexPrintf("%s:%d -- S EDGE s -> %d W = %g\n",
        //          __FILE__, __LINE__, node, qrk);
      } else {
        gp->add_tweights(node, 0, -qrk);
        debugEdges.emplace_back(node + 1, -2, -qrk, 0);
        //mexPrintf("%s:%d -- T EDGE %d -> t W = %g\n",
        //          __FILE__, __LINE__, node, -qrk);
      }
    }
  }

  // Pairwise weights
  int rr, rNode, rrNode;
  double inside, arr;
  for (int r = 0; r < nNodes; r++) {
    for (const Edge &e : neighbors[r]) {
      rr = e.dst;

      for (int k = 0; k < nodes[r]->nStates - 1; k++) {
        rNode = offsets[r] + k;

        for (int kk = 0; kk < nodes[rr]->nStates - 1; kk++) {
          rrNode = offsets[rr] + kk;
          inside = e.p(k,kk) + e.p(k+1,kk+1) - e.p(k+1,kk) - e.p(k,kk+1);

          // We may have small perturbations
          if (fabs(inside) < TEN_EPS) {
            inside = 0.0;
          }

          arr = -(e.w * inside) / 2.0;
          // Will be true if problem was submodular.
          assert(arr >= 0);

          gp->add_edge(rNode, rrNode, arr, 0.0);
          debugEdges.emplace_back(rNode + 1, rrNode + 1, arr, 0.0);
          //mexPrintf("%s:%d -- EDGE %d -> %d FW = %g RW = %g\n",
          //          __FILE__, __LINE__, rNode, rrNode, arr, 0.0);
        }
      }
    }
  }

  //mexPrintf("%s:%d -- Before calling gp->maxflow()\n", __FILE__, __LINE__);
  maxFlow = gp->maxflow();

  std::vector<bool> cut(nBKNodes);
  for (int n = 0; n < nBKNodes; n++) {
    cut[n] = gp->what_segment(n) == dGraph::SINK;
  }


  //mexPrintf("%s:%d -- Before assigning cut\n", __FILE__, __LINE__);
  x.resize(nNodes);
  for (int r = 0; r < nNodes; r++) {
    int nst = nodes[r]->nStates;
    int k = 0;
    while (k < nst - 1 && !cut[offsets[r] + k]) {
      k++;
    }
    x[r] = k;
  }

  // compute energy
  energy[1] = 0.0;
  energy[2] = 0.0;
  for (size_t src = 0; src < nNodes; src++) {
    for (Edge &e : neighbors[src]) {
      size_t dst = e.dst;
      if (src > dst) { // count just once
        energy[2] += e.w * e.p(x[src], x[dst]);
      }
    }

    energy[1] += (*nodes[src])(x[src]);
  }
  energy[0] = energy[1] + energy[2];
}

