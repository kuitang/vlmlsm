#pragma once
#include <mex.h>
#include <boost/iterator/iterator_facade.hpp>

struct SparseEntry {
  mwIndex r;
  mwIndex c;
  double w;
};

class SparseMatrixIter :
  public boost::iterator_facade<
    SparseMatrixIter,
    SparseEntry,
    boost::bidirectional_iterator_tag> {
public:
  SparseMatrixIter() = delete;
  SparseMatrixIter(const mxArray *M)
private:
    mxArray *M;
    // Linear index (into ir and pr)
    mwIndex idx;
    mwIndex col;
    friend class boost::iterator_core_access;

    // Boost iterator_facade requirements
    const reference dereference() const {
    }

    bool equal(const SparseMatrixIter &x) const {
      return M == x.M && idx == x.idx && col == x.col;
    }

    void increment() {
    }

    void decrement() {
    }


};

