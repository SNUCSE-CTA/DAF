#ifndef MAXIMUM_MATCHING_H_
#define MAXIMUM_MATCHING_H_

#include <limits>

#include "global/global.h"
#include "include/backtrack_helper.h"
#include "include/candidate_space.h"
#include "include/data_graph.h"
#include "include/query_graph.h"

namespace daf {
class MaximumMatching {
 public:
  inline MaximumMatching(const DataGraph &data, const QueryGraph &query,
                         const CandidateSpace &cs,
                         BacktrackHelper *backtrack_helpers);
  inline ~MaximumMatching();

  MaximumMatching &operator=(const MaximumMatching &) = delete;
  MaximumMatching(const MaximumMatching &) = delete;

  inline Size ComputeMaximumMatching(Label label);

  inline void AddToU(Size u);
  inline void AddToV(Size v);

  inline bool IsAddedToV(Size v);

  inline void Clear(Size *nec_distinct_cands, Size *num_nec_distinct_cands);

 private:
  const DataGraph &data_;
  const QueryGraph &query_;
  const CandidateSpace &cs_;
  BacktrackHelper *backtrack_helpers_;

  Size NIL = 0;
  Size INF = std::numeric_limits<Size>::max();

  Size *pair_U;
  Size *pair_V;

  Size *cand_to_v;

  Size *dist;
  Size *queue;

  Size *nec_index;

  Size size_U;
  Size size_V;
  Size max_label_counter;

  inline bool BFS();
  inline bool DFS(Size u);
};

inline MaximumMatching::MaximumMatching(const DataGraph &data,
                                        const QueryGraph &query,
                                        const CandidateSpace &cs,
                                        BacktrackHelper *backtrack_helpers)
    : data_(data),
      query_(query),
      cs_(cs),
      backtrack_helpers_(backtrack_helpers) {
  pair_U = new Size[query_.GetNumVertices() + 1];
  pair_V = new Size[data_.GetMaxLabelFrequency() + 1];

  dist = new Size[query_.GetNumVertices() + 1];
  queue = new Size[query_.GetNumVertices() + 1];

  nec_index = new Size[query_.GetNumVertices() + 1];

  cand_to_v = new Size[data_.GetNumVertices()];

  max_label_counter = data_.GetMaxLabelFrequency();
  size_U = 0;
  size_V = 0;

  std::fill(cand_to_v, cand_to_v + data_.GetNumVertices(), -1);
}

inline MaximumMatching::~MaximumMatching() {
  delete[] pair_U;
  delete[] pair_V;

  delete[] cand_to_v;

  delete[] dist;
  delete[] queue;

  delete[] nec_index;
}

inline Size MaximumMatching::ComputeMaximumMatching(Label label) {
  Size maximum_matching = 0;

  std::fill(pair_U, pair_U + query_.GetNumVertices() + 1, NIL);
  std::fill(pair_V, pair_V + data_.GetLabelFrequency(label) + 1, NIL);

  // compute maximum matching with Hopcroft-Karp algorithm
  while (BFS()) {
    for (Size u = 1; u <= size_U; ++u) {
      if (pair_U[u] == NIL) {
        if (DFS(u)) maximum_matching += 1;
      }
    }
  }

  return maximum_matching;
}

inline void MaximumMatching::AddToU(Size u) {
  size_U += 1;
  nec_index[size_U] = u;
}

inline void MaximumMatching::AddToV(Size v) { /* code */ }

inline bool MaximumMatching::IsAddedToV(Size v) { /* code */ }

inline void MaximumMatching::Clear(Size *nec_distinct_cands,
                                   Size *num_nec_distinct_cands) { /* code */ }

inline bool MaximumMatching::BFS() { /* code */ }

inline bool MaximumMatching::DFS(Size u) { /* code */ }
}  // namespace daf

#endif  // MAXIMUM_MATCHING_H_
