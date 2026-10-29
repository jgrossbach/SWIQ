// C++ implementation of the per-marker part of get_importances().
// Reproduces the R code bit for bit (same operations, same order of
// floating-point operations), so p-values are identical to the R version.
//
// Compile once per R session with:  Rcpp::sourceCpp("swiq_importances.cpp")

// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <cmath>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

namespace {

// ranger forest stored as flat arrays, one block of nodes per tree
struct FlatForest {
  std::vector<size_t> offset;  // index of node 0 of each tree
  std::vector<size_t> n_nodes;
  std::vector<int> left, right, var;
  std::vector<double> val;
};

FlatForest flatten_forest(const List& child_nodeIDs, const List& split_varIDs,
                          const List& split_values) {
  FlatForest f;
  const R_xlen_t T = child_nodeIDs.size();
  f.offset.resize(T);
  f.n_nodes.resize(T);
  size_t total = 0;
  for (R_xlen_t t = 0; t < T; ++t) {
    NumericVector sv = split_values[t];
    f.offset[t] = total;
    f.n_nodes[t] = sv.size();
    total += sv.size();
  }
  f.left.resize(total);
  f.right.resize(total);
  f.var.resize(total);
  f.val.resize(total);
  for (R_xlen_t t = 0; t < T; ++t) {
    List cn = child_nodeIDs[t];
    std::vector<int> l = as<std::vector<int> >(cn[0]);
    std::vector<int> r = as<std::vector<int> >(cn[1]);
    std::vector<int> v = as<std::vector<int> >(split_varIDs[t]);
    std::vector<double> s = as<std::vector<double> >(split_values[t]);
    const size_t off = f.offset[t], nn = f.n_nodes[t];
    if (l.size() != nn || r.size() != nn || v.size() != nn)
      stop("Inconsistent node vectors in ranger forest");
    for (size_t k = 0; k < nn; ++k) {
      f.left[off + k] = l[k];
      f.right[off + k] = r[k];
      f.var[off + k] = v[k];
      f.val[off + k] = s[k];
    }
  }
  return f;
}

// Split decision exactly as in ranger's Tree::predict()
inline bool goes_left(double value, double split_value, bool ordered) {
  if (ordered) {
    return value <= split_value;
  }
  size_t factorID = std::floor(value) - 1;
  size_t splitID = std::floor(split_value);
  return !(splitID & (1ULL << factorID));
}

// Walk with the original data from the root until a terminal node or the
// first node that splits on marker m. Everything above that node is the
// same for every allele of m.
inline size_t walk_to_marker(const FlatForest& f, size_t off, const double* X,
                             size_t n, size_t s, int m,
                             const std::vector<int>& is_ord, bool& terminal) {
  size_t node = 0;
  while (true) {
    const int l = f.left[off + node], r = f.right[off + node];
    if (l == 0 && r == 0) {
      terminal = true;
      return node;
    }
    const int v = f.var[off + node];
    if (v == m) {
      terminal = false;
      return node;
    }
    const double value = X[s + (size_t)v * n];
    node = goes_left(value, f.val[off + node], is_ord[v]) ? l : r;
  }
}

// Continue from `node` with column m replaced by `mval`
inline size_t walk_from(const FlatForest& f, size_t off, size_t node,
                        const double* X, size_t n, size_t s, int m, double mval,
                        const std::vector<int>& is_ord) {
  while (true) {
    const int l = f.left[off + node], r = f.right[off + node];
    if (l == 0 && r == 0) return node;
    const int v = f.var[off + node];
    const double value = (v == m) ? mval : X[s + (size_t)v * n];
    node = goes_left(value, f.val[off + node], is_ord[v]) ? l : r;
  }
}

// Equality as used by R's unique(): NaN equals NaN, NA equals NA,
// but NaN does not equal NA
inline bool r_equal(double a, double b) {
  if (!ISNAN(a) && !ISNAN(b)) return a == b;
  if (R_IsNA(a) && R_IsNA(b)) return true;
  if (R_IsNaN(a) && R_IsNaN(b)) return true;
  return false;
}

}  // namespace

// Which variables are split on in which tree; same as
// (1:p) %in% (treeInfo(rf, t)$splitvarID + 1) for every tree t
// [[Rcpp::export]]
LogicalMatrix swiq_split_var_use(const List& child_nodeIDs,
                                 const List& split_varIDs, int p) {
  const R_xlen_t T = child_nodeIDs.size();
  LogicalMatrix out(p, T);
  for (R_xlen_t t = 0; t < T; ++t) {
    List cn = child_nodeIDs[t];
    std::vector<int> l = as<std::vector<int> >(cn[0]);
    std::vector<int> v = as<std::vector<int> >(split_varIDs[t]);
    for (size_t k = 0; k < l.size(); ++k) {
      // treeInfo() treats a node as terminal when its left child is 0
      if (l[k] != 0 && v[k] >= 0 && v[k] < p) out(v[k], t) = TRUE;
    }
  }
  return out;
}

// Per-marker importances; returns a p x n matrix (markers x samples),
// identical to do.call("rbind", imps) in the R version.
// [[Rcpp::export]]
NumericMatrix swiq_marker_importances(const List& child_nodeIDs,
                                      const List& split_varIDs,
                                      const List& split_values,
                                      const LogicalVector& is_ordered,
                                      const NumericMatrix& X,
                                      const NumericVector& y,
                                      const IntegerMatrix& preds_real,
                                      const LogicalMatrix& oob,
                                      const LogicalMatrix& m_use,
                                      const List& alleles,
                                      const NumericVector& scale,
                                      int num_threads) {
  const size_t n = X.nrow();
  const int p = X.ncol();
  const FlatForest f = flatten_forest(child_nodeIDs, split_varIDs, split_values);
  const size_t T = f.offset.size();

  if ((size_t)y.size() != n) stop("y and X differ in number of samples");
  if ((size_t)preds_real.nrow() != n || (size_t)preds_real.ncol() != T)
    stop("preds_real must be n x number of trees");
  if ((size_t)oob.nrow() != n || (size_t)oob.ncol() != T)
    stop("oob must be n x number of trees");
  if (m_use.nrow() != p || (size_t)m_use.ncol() != T)
    stop("m_use must be number of predictors x number of trees");
  if (alleles.size() != p || scale.size() != p)
    stop("alleles and scale need one entry per predictor");
  if (is_ordered.size() != p) stop("is_ordered needs one entry per predictor");

  // plain C++ copies of everything the threads read
  std::vector<int> is_ord(p);
  for (int j = 0; j < p; ++j) is_ord[j] = is_ordered[j] == TRUE;
  const double* Xp = X.begin();
  const double* yp = y.begin();
  const int* pr = preds_real.begin();
  const int* oobp = oob.begin();

  std::vector<std::vector<double> > allele_vals(p);
  for (int j = 0; j < p; ++j) allele_vals[j] = as<std::vector<double> >(alleles[j]);

  std::vector<std::vector<int> > tree_list(p);
  for (int j = 0; j < p; ++j)
    for (size_t t = 0; t < T; ++t)
      if (m_use(j, t) == TRUE) tree_list[j].push_back((int)t);

  // prediction of each terminal node (node_pred) and number of OOB samples
  // whose original data end in each terminal node (node_size_oob)
  const size_t total_nodes = f.left.size();
  std::vector<double> node_pred(total_nodes, NA_REAL);
  std::vector<int> node_size_oob(total_nodes, 0);
  size_t max_nodes = 0;
  for (size_t t = 0; t < T; ++t) {
    const size_t off = f.offset[t];
    if (f.n_nodes[t] > max_nodes) max_nodes = f.n_nodes[t];
    for (size_t k = 0; k < f.n_nodes[t]; ++k)
      if (f.left[off + k] == 0) node_pred[off + k] = f.val[off + k];
    for (size_t s = 0; s < n; ++s)
      if (oobp[s + t * n]) node_size_oob[off + pr[s + t * n]]++;
  }

  NumericMatrix out(p, n);
  double* outp = out.begin();
  const double* scalep = scale.begin();

#ifdef _OPENMP
  if (num_threads < 1) num_threads = 1;
#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
#endif
  for (int m = 0; m < p; ++m) {
    const std::vector<int>& trees = tree_list[m];
    const size_t Tm = trees.size();
    if (Tm == 0) continue;  // marker never used: importance 0
    const std::vector<double>& a = allele_vals[m];
    const size_t K = a.size();
    if (K == 0) continue;
    const double* xm = Xp + (size_t)m * n;

    // allele of each sample (last match wins, as in the R loop) and counts
    std::vector<int> own(n, -1);
    std::vector<int> cnt(K, 0);
    for (size_t s = 0; s < n; ++s)
      for (size_t i = 0; i < K; ++i)
        if (xm[s] == a[i]) {
          own[s] = (int)i;
          cnt[i]++;
        }

    std::vector<int> nodes(K * n), pops(K * n);
    std::vector<double> err(K * n);
    std::vector<double> sums(max_nodes);
    std::vector<int> counts(max_nodes);
    std::vector<long double> rans(n, 0.0L);  // as in R's rowMeans()

    for (size_t jt = 0; jt < Tm; ++jt) {
      const size_t t = trees[jt];
      const size_t off = f.offset[t];
      const size_t nn = f.n_nodes[t];
      const int* oob_t = oobp + t * n;
      const int* pr_t = pr + t * n;

      // terminal node of every sample for every allele of marker m
      for (size_t s = 0; s < n; ++s) {
        bool terminal;
        const size_t u = walk_to_marker(f, off, Xp, n, s, m, is_ord, terminal);
        for (size_t i = 0; i < K; ++i) {
          if (terminal) {
            nodes[i * n + s] = (int)u;
          } else if (own[s] == (int)i) {
            nodes[i * n + s] = pr_t[s];  // sample already carries allele i
          } else {
            nodes[i * n + s] = (int)walk_from(f, off, u, Xp, n, s, m, a[i], is_ord);
          }
        }
      }

      for (size_t i = 0; i < K; ++i) {
        const int* nd = &nodes[i * n];
        double* e = &err[i * n];
        int* pp = &pops[i * n];
        for (size_t s = 0; s < n; ++s) pp[s] = node_size_oob[off + nd[s]];

        if (cnt[i] == 0) {
          for (size_t s = 0; s < n; ++s) e[s] = 0.0;
          continue;
        }
        for (size_t s = 0; s < n; ++s) {
          const double d = node_pred[off + nd[s]] - yp[s];
          e[s] = d * d;
        }
        // subtract the mean error of the OOB samples in the same node
        for (size_t k = 0; k < nn; ++k) {
          sums[k] = 0.0;
          counts[k] = 0;
        }
        for (size_t s = 0; s < n; ++s)
          if (oob_t[s]) {
            sums[nd[s]] += e[s];
            counts[nd[s]]++;
          }
        for (size_t s = 0; s < n; ++s)
          e[s] = e[s] - sums[nd[s]] / (double)counts[nd[s]];
      }

      for (size_t s = 0; s < n; ++s) {
        int allele_sum = 0;
        for (size_t i = 0; i < K; ++i) allele_sum += pops[i * n + s];
        double blended = 0.0;
        for (size_t i = 0; i < K; ++i) {
          // volatile stops the compiler from fusing multiply and add (FMA),
          // which R never does and which would change the rounding
          volatile double prod =
              err[i * n + s] * ((double)pops[i * n + s] / (double)allele_sum);
          blended = blended + prod;
        }
        const double real = own[s] >= 0 ? err[(size_t)own[s] * n + s] : NA_REAL;
        double diff = blended - real;
        if (!oob_t[s]) {
          diff = 0.0;
        } else {
          bool same = true;
          for (size_t i = 1; i < K; ++i)
            if (!r_equal(err[i * n + s], err[s])) {
              same = false;
              break;
            }
          if (same) diff = 0.0;
        }
        rans[s] += diff;
      }
    }

    for (size_t s = 0; s < n; ++s) {
      const double row_mean = (double)(rans[s] / (long double)Tm);
      outp[m + s * (size_t)p] = row_mean * scalep[m];
    }
  }
  return out;
}
