This code exactly optimizes functions of the form

  x = \argmin \sum_i D_i(x_i) + \sum_ij w_ij V_m (x_i, x_j)

Using Schlensinger and Flach's multi-label to binary reduction of
the min sum problem and Boykov and Kolmogorov's graph cut code. It
corrects and extends the original implementation of Shai Bagon. In
particular, nodes need not have the same label cardinality.

See the header on MultiLabelSubModular.m for additional copyright
and citation information.

