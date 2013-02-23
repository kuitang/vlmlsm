#pragma once
#include <algorithm>
#include <cstdio>
#include <vector>
#include <utility>
#include "graph.h"
#include "mex.h"
#include "matrix.h"

struct MinSum;

// Assumes that the potential data is stored immediately after the node.
struct alignas(double) Node {
  // Used a signed type to detection < subtraction fail.
  int64_t nStates;

  Node(int64_t nStates_) : nStates(nStates_) { }

  // Node must live exactly where it was allocated, or bad things will happen.
  // (We will touch memory that we haven't allocated.)
  Node(const Node &other) = delete;
  Node &operator=(const Node &other) = delete;

  double &operator()(int i) {
    mxAssert(i < nStates, "Node operator() out of bounds");
    double *p = reinterpret_cast<double *>(this);
    p += sizeof(Node) / sizeof(double) + i; // units of doubles
    //mexPrintf("%s:%d -- Node %lx at %d is %g at %lx\n", __FILE__, __LINE__, this, i, *p, p);
    return *p;
  }
};

static_assert(sizeof(Node) % sizeof(double) == 0, "Node cannot be stored in double array");

struct alignas(double) Potential {
  int64_t nLo;
  int64_t nHi;

  Potential(int nLo_, int nHi_) : nLo(nLo_), nHi(nHi_) { }

  Potential(const Potential &other) = delete;
  Potential &operator=(const Potential &other) = delete;

  // Column major
  double &operator()(int lo, int hi) {
    mexPrintf("%s:%d -- Potential %lx (nLo = %d, nHi = %d), at (%d, %d)\n", __FILE__, __LINE__, this, nLo, nHi, lo, hi);
    mxAssert(lo < nLo && hi < nHi, "Potential::operator() out of bounds");
    double *p = reinterpret_cast<double *>(this);
    p += sizeof(Potential) / sizeof(double) + nLo*hi + lo;
    return *p;
  }
};

static_assert(sizeof(Potential) % sizeof(double) == 0, "Potential cannot be stored in double array");

struct Edge {
  // If E \in edges[i], then E.dst = i.
  int dst;
  const bool srcHigher;
  double w;
  Potential *potential;

  Edge(int dst_, bool srcHigher_, double w_, Potential *potential_) :
    dst(dst_),
    srcHigher(srcHigher_),
    w(w_),
    potential(potential_) { }

  // Evaluate the potential, transposing if necessary.
  double &p(size_t srcState, size_t dstState) const {
    // I'm hoping the branch predictor catches this
    if (srcHigher) {
      return (*potential)(dstState, srcState);
    } else {
      return (*potential)(srcState, dstState);
    }
  }
};

struct DebugEdge {
  int src;
  int dst;
  double fw;
  double rw;

  DebugEdge(int s, int d, double f, double r) : src(s), dst(d), fw(f), rw(r) { }
};

struct MinSum {
  int nNodes;
  int nEdges;

  // Undirected edge
  std::vector<std::vector<Edge>> neighbors;
  std::vector<Potential *> potentials;
  std::vector<Node *> nodes;
  std::vector<DebugEdge> debugEdges;

  void (*errFunc)(const char *);

  // from BK
  typedef Graph<double, double, double> dGraph;

  // Fixed length allocation
  std::vector<double> mem;
  size_t memUsed;

  // Allocate from mem aligned to type T, with elements at end for
  // nDoubles doubles. To use with Node and Potential.
  template <class T>
  T *alloc(size_t nDoubles) {
    static_assert(sizeof(T) % sizeof(double) == 0, "T cannot be stored in double array");

    size_t structDoubles = sizeof(T) / sizeof(double);
    size_t nextMemUsed   = memUsed + structDoubles + nDoubles;

    if (nextMemUsed >= mem.size()) {
      // Perhaps round by powers of 2?
      mem.resize(2 * nextMemUsed);
    }

    size_t off = memUsed;
    memUsed = nextMemUsed;

    return reinterpret_cast<T *>(mem.data() + off);
  }

  inline Node &getNode(size_t i) { return *nodes[i]; }
  inline Potential &getPotential(size_t i) { return *potentials[i]; }

  // Construct MinSum with an initial number of nodes, a initial memory size
  // to store nodes and potentials of initMem DOUBLES (not bytes) and an
  // error message. Note that nNodes must be fixed ahead of time.
  //
  // The fixed nNodes restriction is not difficult to remove, but incurs a
  // small runtime cost. Perhaps remove in future.
  MinSum(size_t nNodes_,
         size_t nEdges_=1000,
         void (*errFunc_)(const char*) = (void (*)(const char *)) puts,
         size_t initMem=10000) :
    nNodes(nNodes_),
    nodes(nNodes_),
    neighbors(nNodes_),
    nEdges(nEdges_),
    errFunc(errFunc_),
    mem(initMem),
    memUsed(0)
  {
    potentials.reserve(nEdges);
  }

  Node &addNode(size_t n, int nStates_) {
    nodes[n] = new (alloc<Node>(nStates_)) Node(nStates_);
    return *(nodes[n]);
  }

  // Potentials are always indexed lo on the rows, hi on the columns.
  void addEdge(size_t src, size_t dst, double w, Potential *ep) {
    size_t lo = std::min(src, dst);
    size_t hi = std::max(src, dst);

    mxAssert(lo >= 0 && hi < nNodes, "addEdge index out of bounds");

    neighbors[lo].emplace_back(hi, /*srcHigher*/ false, w, ep);
    neighbors[hi].emplace_back(lo, /*srcHigher*/ true,  w, ep);

    nEdges++;
  }

  // should be fast because this happens a lot
  Potential &addPotential(int nLo, int nHi) {
    mxAssert(nLo > 0 && nHi > 0, "Potential was allocated without entries!");
    // Allocate from our memory stash
    Potential *p = new (alloc<Potential>(nLo * nHi)) Potential(nLo, nHi);
    potentials.push_back(p);
    return *p;
  }

  enum class FailCode {
    SUCCESS = 0,
    NOT_SUBMODULAR,
    NOT_SYMMETRIC_W,
    NOT_TRANSPOSED_POTENTIAL,
    SELF_LOOP,
    INSUFFICIENT_STATES,
  };

  bool isSubModular(Potential *p);
  FailCode validate();
  // energy is 3-vector; [0] is total, [1] is unary, [2] is pairwise.
  void minimize(std::vector<int> &x, double *energy, double &maxFlow);
};

