#include "include/backtrack.h"

#include <iostream>

namespace daf {
Backtrack::Backtrack(const DataGraph &data, const QueryGraph &query,
                     const CandidateSpace &cs)
    : data_(data), query_(query), cs_(cs) {
  num_embeddings_ = 0;
  num_backtrack_calls_ = 0;
  backtrack_depth_ = 0;

  mapped_query_vtx_ = new Vertex[data_.GetNumVertices()];
  node_stack_ = new SearchTreeNode[query_.GetNumNonLeafVertices() + 1];
  mapped_nodes_ = new SearchTreeNode *[query_.GetNumVertices()];
  helpers_ = new BacktrackHelper[query_.GetNumVertices()];

  extendable_queue_ = new Ordering(query_.GetNumVertices());

  if (query_.GetNumNonLeafVertices() < query_.GetNumVertices()) {
    match_leaves_ = new MatchLeaves(data_, query_, cs_, mapped_query_vtx_,
                                    helpers_, mapped_nodes_);
  } else {
    match_leaves_ = nullptr;
  }

  std::fill(mapped_query_vtx_, mapped_query_vtx_ + data_.GetNumVertices(), -1);

  for (Vertex u = 0; u < query_.GetNumVertices(); ++u) {
    helpers_[u].Initialize(query_.GetNumVertices(), query_.GetDegree(u),
                           cs_.GetCandidateSetSize(u), u);
  }
}

Backtrack::~Backtrack() { /* code */ }

uint64_t Backtrack::FindMatches(uint64_t limit) { /* code */ }

Vertex Backtrack::GetRootVertex() { /* code */ }

void Backtrack::InitializeNodeStack() { /* code */ }

void Backtrack::ComputeExtendable(Vertex u, Vertex u_nbr, Size u_nbr_idx,
                                  Size cs_v_idx) { /* code */ }

void Backtrack::ComputeDynamicAncestor(Vertex ancsetor, Vertex child) { /* code */ }

bool Backtrack::ComputeExtendableForAllNeighbors(SearchTreeNode *cur_node,
                                                 Size cs_v_idx) { /* code */ }

void Backtrack::ReleaseNeighbors(SearchTreeNode *cur_node) { /* code */ }
}  // namespace daf
