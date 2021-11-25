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
      backtrack_mapped_nodes_(mapped_nodes) { /* code */ }

MatchLeaves::~MatchLeaves() { /* code */ }

uint64_t MatchLeaves::Match(uint64_t limit) { /* code */ }

uint64_t MatchLeaves::Combine() { /* code */ }

void MatchLeaves::ReserveVertex(Vertex represent,
                                BacktrackHelper *repr_helper) { /* code */ }

void MatchLeaves::ClearMemoryForBacktrack() { /* code */ }
}  // namespace daf
