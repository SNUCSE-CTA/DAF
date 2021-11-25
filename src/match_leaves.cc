#include "include/match_leaves.h"

namespace daf {
namespace {
constexpr uint64_t Factorial(Size n) {
  return n == 1 ? 1 : n * Factorial(n - 1);
}
}  // namespace

MatchLeaves::MatchLeaves(const DataGraph &data, const QueryGraph &query,
                         const CandidateSpace &cs, Vertex *mapped_query_vtx,
                         BacktrackHelper *helpers,
                         SearchTreeNode **mapped_nodes)
    : data_(data),
      query_(query),
      cs_(cs),
      backtrack_mapped_query_vtx(mapped_query_vtx),
      backtrack_helpers_(helpers),
      backtrack_mapped_nodes_(mapped_nodes) {
  nec_distinct_cands_ = new Vertex[data_.GetNumVertices()];
  cand_to_nec_idx_ = new std::vector<Size>[data_.GetNumVertices()];
  nec_ranking_ = new Size[query_.GetNumNECLabel()];

  sum_nec_cands_size_ = new Size[query_.GetNumNECLabel()];
  sum_nec_size_ = new Size[query_.GetNumNECLabel()];

  num_remaining_cands_ = new Size[query_.GetNumVertices()];
  num_remaining_nec_vertices_ = new Size[query_.GetNumVertices()];

  num_nec_distinct_cands_ = 0;

  maximum_matching_ =
      new MaximumMatching(data_, query_, cs_, backtrack_helpers_);
}

MatchLeaves::~MatchLeaves() { /* code */ }

uint64_t MatchLeaves::Match(uint64_t limit) { /* code */ }

uint64_t MatchLeaves::Combine() { /* code */ }

void MatchLeaves::ReserveVertex(Vertex represent,
                                BacktrackHelper *repr_helper) { /* code */ }

void MatchLeaves::ClearMemoryForBacktrack() { /* code */ }
}  // namespace daf
