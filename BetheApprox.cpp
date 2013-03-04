#include "BetheApprox.h"
#include <algorithm>
#include <limits>
#include "mex.h"

// Return true if convergence threshold met
// Output are A, B, and alpha, which must be preallocated!
bool propogateBetheBound(size_t nNodes,
                         const double *theta,
                         const cscMatrix &W,
                         double thresh,
                         int maxIter,
                         double *A, // outputs
                         double *B,
                         double *alpha) {
  double posW[nNodes], negW[nNodes];
  for (size_t j = 0; j < nNodes; j++) {
    posW[j] = 0;
    negW[j] = 0;

    for (int idx = W.jc[j]; idx < W.jc[j+1]; idx++) {
      double w = W.pr[idx];
      if (w > 0) {
        posW[j] += w;
      } else {
        negW[j] -= w;
      }

      alpha[idx] = exp(fabs(w)) - 1;
    }

    A[j] = sigmoid(theta[j] - negW[j]);
    B[j] = 1 - sigmoid(theta[j] + posW[j]);
  }

  bool converged = false;
  double oldA[nNodes], oldB[nNodes];
  double L, U;

  for (int iter = 0; !converged; iter++) {
    std::copy(A, A + nNodes, oldA);
    std::copy(B, B + nNodes, oldB);

    for (size_t j = 0; j < nNodes; j++) {
      L = 1.0;
      U = 1.0;
      for (int idx = W.jc[j]; idx < W.jc[j+1]; idx++) {
        int i = W.ir[idx];
        double a = alpha[idx];

        if (W.pr[idx] > 0) {
          L = L * (1 + a*A[i] / (1 + a*(1 - B[j])*(1 - A[i])));
          U = U * (1 + a*B[i] / (1 + a*(1 - A[j])*(1 - B[i])));
        } else {
          L = L * (1 + a*B[i] / (1 + a*(1 - B[j])*(1 - B[i])));
          U = U * (1 + a*A[i] / (1 + a*(1 - A[j])*(1 - A[i])));
        }
      }

      A[j] = 1 / (1 + exp(-theta[j] + negW[j]) / L);
      B[j] = 1 / (1 + exp(theta[j]  + posW[j]) / U);
    }

    converged = oneNormConverged(nNodes, A, oldA, thresh) &&
                oneNormConverged(nNodes, B, oldB, thresh);
  }

  return converged;
}

inline int degree(const cscMatrix &W, int j) {
  return W.jc[j+1] - W.jc[j];
}

// For node n, find points in the interval [A[n], 1 - B[n]] such that the
// distance between two consecutive points is at most intervalSz.
//
// Returns a cscMatrix. Remember to free[] the pointers yourself.
cscMatrix calcIntervals(size_t nNodes, const double *A, const double *B, double intervalSz) {
  // Precalculate my number of intervals
  size_t totPoints = 0;
  size_t maxPoints = 0;

  // Number of endpoints is ceil(intervalLength / intervalSz) + 1 (for endpoint)
  //
  // i.e.
  // -----------
  // |  |  |  ||
  //
  // With floating point numbers, we will assume that the endpoint is never
  // reached by stepping with intervalSz.
  //
  // We might have intervalLength < intervalSz. In that case, we still need
  // the left endpoint, hence the max.
  for (size_t n = 0; n < nNodes; n++) {
    double intervalLength = 1 - B[n] - A[n];
    size_t points = std::max((int) std::ceil(intervalLength / intervalSz), 1) + 1;
    if (points > maxPoints) { maxPoints = points; }
    totPoints += points;
  }

  // maxPoints x nNodes sparse matrix
  cscMatrix m;
  m.M = maxPoints;
  m.N = nNodes;
  m.nzMax = totPoints;
  m.pr = new double[totPoints];
  m.ir = new size_t[totPoints];
  m.jc = new size_t[m.N+1];

  size_t ind = 0;
  for (size_t n = 0; n < nNodes; n++) {
    m.jc[n] = ind;
    size_t k = 0;
    mxAssert(A[n] < 1 - B[n], "bounds failed.");

    for (double q = A[n]; q < 1 - B[n]; q += intervalSz) {
      mxAssert(ind < totPoints, "You fucked up calculating totPoints!");
      m.ir[ind] = k;
      m.pr[ind] = q;
      ind++;
      k++;
    }

    // Add the endpoint
    m.ir[ind] = k;
    m.pr[ind] = 1 - B[n];
    ind++;
  }
  // Last column; one-past-end
  // We don't need to add 1, because the ind++ above means that at the last
  // iterations, we were already one past the end.
  m.jc[nNodes] = ind;

  return m;
}

double getIntervalSz(size_t nNodes,
                     const double *A,
                     const double *B,
                     const cscMatrix &W,
                     const double *alpha,
                     double epsilon) {
  // (Eq 16)
  double eta[nNodes];
  for (size_t n = 0; n < nNodes; n++) {
    eta[n] = std::min(A[n], B[n]);
  }

  double aMax = -std::numeric_limits<double>::min();
  double aa, aij;
  for (size_t j = 0; j < nNodes; j++) {
    for (int idx = W.jc[j]; idx < W.jc[j+1]; idx++) {
      int i = W.ir[idx];

      // Upper triangular
      if (j > i) {
        aij  = alpha[idx];
        aa = aij*(aij + 1) / (4*(2*aij + 1)*eta[i]*eta[j]*(1 - eta[i])*(1 - eta[j]));
        if (aa > aMax) {
          aMax = aa;
        }
      }
    }
  }

  double bb, bInner, alphaSum;
  double bMax = -std::numeric_limits<double>::min();
  for (size_t j = 0; j < nNodes; j++) {
    alphaSum = 0.0;
    // loop over the neighbors
    for (size_t idx = W.jc[j]; idx < W.jc[j+1]; idx++) {
      double a = alpha[idx];
      alphaSum += (a + 1)*(a + 1) / (2*a + 1);
    }

    bInner = 1 - degree(W, j) + alphaSum;
    bb = 1 / (eta[j] * (1 - eta[j])) * bInner;
    if (bb > bMax) {
      bMax = bb;
    }
  }

  double Omega = std::max(aMax, bMax);
  // cast inside to force double division
  double density = (W.nzMax + nNodes) / double(W.N * W.M);

  return sqrt(2*epsilon / (nNodes * Omega * sqrt(density)));
}

void makeBetheMinSum(size_t nNodes,
                     const double *theta,
                     const cscMatrix &W,
                     const double *A,
                     const double *B,
                     const double *alpha,
                     double intervalSz,
                     MinSum &m) {
  mxAssert(m.nNodes == nNodes, "nNodes was a lie.");

  cscMatrix intervals = calcIntervals(nNodes, A, B, intervalSz);

  // Unaries: Term two of (Eq 4)
  for (size_t j = 0; j < nNodes; j++) {
    int nStates = intervals.jc[j+1] - intervals.jc[j];
    mxAssert(nStates >= 2, "states did not exceed two!");
    Node &nj = m.addNode(j, nStates);

    int degMinusOne = degree(W, j) - 1;

    for (size_t k = 0; k < nStates; k++) {
      double q = intervals.pr[intervals.jc[j] + k];
      nj(k) = -theta[j] * q + degMinusOne * binaryEntropy(q);
      //mexPrintf("%s:%d -- Node %d[%d] is %g\n", __FILE__, __LINE__, j, k, nj(k));
    }
  }

  // We know ahead of time there will be nnz distinct potentials.
  m.potentials.reserve(W.nzMax / 2);

  // Pairwise: (Eq 5)
  for (size_t hi = 0; hi < nNodes; hi++) {
    int nHiStates = m.nodes[hi].nStates;
    mxAssert(nHiStates >= 0, "nHiStates cannot be negative.");

    for (size_t nodeIdx = W.jc[hi]; nodeIdx < W.jc[hi+1]; nodeIdx++) {
      int lo = W.ir[nodeIdx];

      if (hi > lo) {
        int nLoStates = m.nodes[lo].nStates;
        mxAssert(nLoStates >= 0, "nLoStates cannot be negative.");

        // Only look at the upper triangular (minus diagonal)
        double w = W.pr[nodeIdx];
        double aij = alpha[nodeIdx];

        mxAssert(aij != 0, "alpha cannot be zero.");

        Potential &potential = m.addPotential(nLoStates, nHiStates);
        m.addEdge(lo, hi, 1.0, &potential);

        for (size_t iqLo = 0; iqLo < nLoStates; iqLo++) {
          double qLo = intervals.pr[intervals.jc[lo] + iqLo];

          for (size_t iqHi = 0; iqHi < nHiStates; iqHi++) {
            double qHi = intervals.pr[intervals.jc[hi] + iqHi];
            double marginals[4];
            double xi = marginalize(aij, qLo, qHi, marginals);
            potential(iqLo, iqHi) = -w*xi - entropy<4>(marginals);

            //mexPrintf("%s:%d -- Potential at %lx entry (%d, %d) is %g\n",
            //          __FILE__, __LINE__, &potential, iqLo, iqHi, potential(iqLo, iqHi));
          }
        }

      }
    }
  }

  deleteCscMatrix(intervals);
}

