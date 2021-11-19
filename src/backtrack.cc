#include "include/backtrack.h"

#include <iostream>

namespace daf {
Backtrack::Backtrack(const DataGraph &data, const QueryGraph &query,
                     const CandidateSpace &cs)
    : data_(data), query_(query), cs_(cs) { /* code */ }

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
