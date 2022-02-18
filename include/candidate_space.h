#ifndef CANDIDATE_SPACE_H_
#define CANDIDATE_SPACE_H_

#include <utility>

#include "global/global.h"
#include "include/dag.h"
#include "include/data_graph.h"
#include "include/query_graph.h"

namespace daf {
class CandidateSpace {
 public:
  CandidateSpace(const DataGraph &data, const QueryGraph &query,
                 const DAG &dag);
  ~CandidateSpace();

  CandidateSpace &operator=(const CandidateSpace &) = delete;
  CandidateSpace(const CandidateSpace &) = delete;

  bool BuildCS();

  inline Size GetCandidateSetSize(Vertex u) const;
  inline Vertex GetCandidate(Vertex u, Size v_idx) const;

  inline Size GetCandidateStartOffset(Vertex u, Size u_adj_idx,
                                      Size v_idx) const;
  inline Size GetCandidateEndOffset(Vertex u, Size u_adj_idx, Size v_idx) const;
  inline Size GetCandidateIndex(Size idx) const;

 private:
  const DataGraph &data_;
  const QueryGraph &query_;
  const DAG &dag_;

  Size *candidate_set_size_;
  Vertex **candidate_set_;
  Size ***candidate_offsets_;
  Vertex *linear_cs_adj_list_;

  QueryDegree *num_visit_cs_;
  Vertex *visited_candidates_;
  Size *cand_to_cs_idx_;
  Size num_visitied_candidates_;

  Size num_cs_edges_;

  bool FilterByTopDownWithInit();
  bool FilterByBottomUp();
  bool FilterByTopDown();

  void ConstructCS();
  bool InitRootCandidates();

  void ComputeNbrInformation(Vertex u, Size *max_nbr_degree,
                             uint64_t *label_set);
};

inline Size CandidateSpace::GetCandidateSetSize(Vertex u) const {
  return candidate_set_size_[u];
}

inline Vertex CandidateSpace::GetCandidate(Vertex u, Size v_idx) const {
  return candidate_set_[u][v_idx];
}

inline Size CandidateSpace::GetCandidateStartOffset(Vertex u, Size u_adj_idx,
                                                    Size v_idx) const {
  return candidate_offsets_[u][u_adj_idx][v_idx];
}

inline Size CandidateSpace::GetCandidateEndOffset(Vertex u, Size u_adj_idx,
                                                  Size v_idx) const {
  return candidate_offsets_[u][u_adj_idx][v_idx + 1];
}

inline Size CandidateSpace::GetCandidateIndex(Size idx) const {
  return linear_cs_adj_list_[idx];
}

}  // namespace daf

#endif  // CANDIDATE_SPACE_H_
