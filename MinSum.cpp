#include <cmath>
#include <cstdio>
#include <numeric>
#include "MinSum.h"

// TODO: errMsg function

const double TEN_EPS = 2.2204e-15;
const double INF = 1e100;

// Checks that potentials are submodular and that the neighbors matrix
// is symmetric and potentials are transposed when needed.

bool MinSum::isSubModular(Potential &q) {
  char errMsg[1000];
  for (int lo = 0; lo < q.nLo - 1; lo++) {
    for (int hi = 0; hi < q.nHi - 1; hi++) {
      double diff = q(lo,hi) + q(lo+1, hi+1) - q(lo, hi+1) - q(lo+1, hi);
      if (diff - 5 * TEN_EPS > 0) {
        snprintf(errMsg, 1000, "isSubModular: GEQ Potential %p entry (%d,%d) failed (%g + %g > %g + %g) (diff = %g)\n", &q, lo, hi, q(lo,hi), q(lo+1,hi+1), q(lo,hi+1), q(lo+1,hi), diff);

        errFunc(errMsg);
        return false;
      }
    }
  }

  return true;
}

Node &MinSum::addNode(size_t n, int nStates_) {
  // Recall that all nodes were preallocated
  mxAssert(n < nNodes, "No space for this node!");

  Node &node = nodes[n];
  node.nStates = nStates_;
  node.vals = alloc(nStates_);
  return node;
}

// Potentials are always indexed lo on the rows, hi on the columns.
void MinSum::addEdge(size_t src, size_t dst, double w, Potential *ep) {
  size_t lo = std::min(src, dst);
  size_t hi = std::max(src, dst);

  mxAssert(lo >= 0 && lo < nNodes && hi >= 0 && hi < nNodes, "addEdge index out of bounds");

  neighbors[lo].emplace_back(hi, /*srcHigher*/ false, w, ep);
  neighbors[hi].emplace_back(lo, /*srcHigher*/ true,  w, ep);

  nEdges++;
}

// should be fast because this happens a lot
Potential &MinSum::addPotential(int nLo, int nHi) {
  mxAssert(nLo > 0 && nHi > 0, "Potential was allocated without entries!");

  double *vals = alloc(nLo * nHi);
  potentials.emplace_back(nLo, nHi, vals);
  return potentials.back();
}

double *MinSum::alloc(size_t nDoubles) {
  size_t nextMemUsed = memUsed + nDoubles;
  if (nextMemUsed >= mem.size()) {
    double *oldBase = mem.data();
    mem.resize(2 * nextMemUsed);
    double *newBase = mem.data();

    printFunc("%s:%d -- Realloc; old base was %lp; new base is %lp; new size (in doubles) is %d\n",
               __FILE__, __LINE__, oldBase, newBase, mem.size());

    if (newBase != oldBase) { // translate all old pointers
      for (Node &node : nodes) {
        if (node.vals != nullptr) {
          ptrdiff_t off = node.vals - oldBase;
          node.vals = newBase + off;
        }
      }
      for (Potential &pot : potentials) {
        if (pot.vals != nullptr) {
          ptrdiff_t off = pot.vals - oldBase;
          pot.vals = newBase + off;
        }
      }
    }
  }

  double *ret = &mem[memUsed];
  memUsed = nextMemUsed;
  return ret;
}


MinSum::FailCode MinSum::validate() {
  for (Potential &p : potentials) {
    if (!isSubModular(p)) { return FailCode::NOT_TRANSPOSED_POTENTIAL;  }
  }

  for (int src = 0; src < nNodes; src++) {
    for (const Edge &e : neighbors[src]) {
      int dst = e.dst;
      if ( (src > dst && !e.srcHigher) ||
           (src < dst && e.srcHigher) ) { return FailCode::NOT_TRANSPOSED_POTENTIAL; }
      if (src == dst) { return FailCode::SELF_LOOP; }
    }
  }

  for (const Node &n : nodes) {
    if (n.nStates < 2) {
      mxAssert(false, "Node had fewer than two states!");
      return FailCode::INSUFFICIENT_STATES;
    }
  }

  mxAssert(nNodes == nodes.size(), "Not enough nodes!");

  return FailCode::SUCCESS;
}

void MinSum::minimize(std::vector<int> &x, double *energy, double &maxFlow) {
  // Compute offsets to convert (r, k) indices to linear
  std::vector<int> offsets;
  offsets.push_back(0);
  for (int r = 1; r < nNodes; r++) {
    offsets.push_back(offsets[r - 1] + nodes[r - 1].nStates - 1);
  }
  mxAssert(offsets.size() == nNodes, "Offsets and nodes epic fail");

  int nBKNodes = offsets.back() + nodes.back().nStates - 1;

//  for (int r = 0; r < nNodes; r++) {
//    mexPrintf("%s:%d -- Node %d@%p had %d represented states. Offset is %d\n",
//              __FILE__, __LINE__, r, nodes[r], nodes[r]->nStates - 1, offsets[r]);
//  }

  // upper bound on number of neighbors:
  // one edge per node for s/t (but those don't count)
  // one edge per state for the zero/infinity (captured by offset)
  // variable neighbors for neighbors in underlying graph
  int nBKneighbors = offsets.back() + nodes.back().nStates;

  for (int r = 0; r < nNodes; r++) {
    for (const Edge &e : neighbors[r]) {
      int rr = e.dst;
      int rStates  = nodes[r].nStates - 2;
      int rrStates = nodes[rr].nStates - 2;
      nBKneighbors += rStates * rrStates;
    }
  }

  // Make our graph
  dGraph *gp = new dGraph(nBKNodes, nBKneighbors, errFunc);
  mexPrintf("%s:%d -- nBKNeighbors = %d\n", __FILE__, __LINE__, nBKneighbors);
  mexPrintf("%s:%d -- nBKNodes = %d\n", __FILE__, __LINE__, nBKNodes);

  // Add nodes
  gp->add_node(nBKNodes);

  // Zero/Infinity weights
  int loNode, hiNode;
  for (int r = 0; r < nodes.size(); r++) {
    Node &currNode = nodes[r];
    //mexPrintf("%s:%d -- r = %d, currNode->nStates = %d\n", __FILE__, __LINE__, r, currNode->nStates);

    for (int k = 0; k < currNode.nStates - 2; k++) {
      //mexPrintf("%s:%d -- r = %d, k = %d, currNode->nStates = %d\n",
      //          __FILE__, __LINE__, r, k, currNode->nStates);
      loNode = offsets[r] + k;
      hiNode = offsets[r] + k + 1;

      mxAssert(loNode >= 0 && hiNode >= 0, "No negative nodes!");
      mxAssert(loNode < nBKNodes && hiNode < nBKNodes, "Node idx out of bounds!");
      gp->add_edge(loNode, hiNode, 0.0, INF);
      debugEdges.emplace_back(loNode + 1, hiNode + 1, 0.0, INF);

      //mexPrintf("%s:%d -- EDGE %d -> %d FW = %g RW = %g\n",
      //          __FILE__, __LINE__, loNode, hiNode, 0.0, INF);
    }
  }

  // Unary source/sink weights
  double qrk;
  for (int r = 0; r < nNodes; r++) {
    for (int k = 0; k < nodes[r].nStates - 1; k++) {
      qrk = 0;

      for (const Edge &e : neighbors[r]) {
        int end = nodes[e.dst].nStates - 1;
        qrk += e.w * (e.p(k,0) + e.p(k,end) - e.p(k+1,0) - e.p(k+1,end));
        //mexPrintf("%s:%d qrk for r = %d, k = %d, rr = %d is now %g\n",
        //          __FILE__, __LINE__, r, k, e.dst, qrk);
        //mexPrintf("%s:%d Components: %g %g %g %g\n", __FILE__, __LINE__, e.p(k,1), e.p(k,end), e.p(k+1,1), e.p(k+1,end));
      }

      qrk = qrk / 2.0;
      qrk = qrk + nodes[r](k) - nodes[r](k+1);
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

      for (int k = 0; k < nodes[r].nStates - 1; k++) {
        rNode = offsets[r] + k;

        for (int kk = 0; kk < nodes[rr].nStates - 1; kk++) {
          rrNode = offsets[rr] + kk;
          //mexPrintf("%s:%d -- r = %d, rr = %d, k = %d, kk = %d, nodes[r]->nStates = %d, nodes[rr]->nStates = %d", __FILE__, __LINE__, r, rr, k, kk, nodes[r]->nStates, nodes[rr]->nStates);
          inside = e.p(k,kk) + e.p(k+1,kk+1) - e.p(k+1,kk) - e.p(k,kk+1);

          // We may have small perturbations
          if (fabs(inside) < 10 * TEN_EPS) {
            inside = 0.0;
          }

          arr = -(e.w * inside) / 2.0;
          // Will be true if problem was submodular.
          if (arr < 0) {
            mexPrintf("%s:%d -- arr = %g\n", __FILE__, __LINE__, arr);
          }
          mxAssert(arr >= 0, "not submodular!");

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
    int nst = nodes[r].nStates;
    int k = 0;
    while (k < nst - 1 && !cut[offsets[r] + k]) {
      k++;
    }
    x[r] = k;
  }

  // compute energy
  energy[1] = 0.0;
  energy[2] = 0.0;
  for (int src = 0; src < nNodes; src++) {
    for (Edge &e : neighbors[src]) {
      int dst = e.dst;
      if (src < dst) { // count just once
        energy[2] += e.w * e.p(x[src], x[dst]);
      }
    }

    energy[1] += nodes[src](x[src]);
  }
  energy[0] = energy[1] + energy[2];

  delete gp;
}

