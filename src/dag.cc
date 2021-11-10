#include "include/dag.h"

#include <algorithm>

namespace daf {
DAG::DAG(const DataGraph &data, const QueryGraph &query)
    : data_(data), query_(query) {
  // initialize
  bfs_sequence_ = new Vertex[query_.GetNumVertices()];

  num_children_ = new Size[query.GetNumVertices()];
  num_parents_ = new Size[query.GetNumVertices()];

  children_ = new Vertex *[query_.GetNumVertices()];
  parents_ = new Vertex *[query_.GetNumVertices()];

  for (Size i = 0; i < query_.GetNumVertices(); ++i) {
    children_[i] = new Vertex[query_.GetNumVertices()];
    parents_[i] = new Vertex[query_.GetNumVertices()];
  }

  init_cand_size_ = new Size[query_.GetNumVertices()];

  std::fill(num_children_, num_children_ + query.GetNumVertices(), 0);
  std::fill(num_parents_, num_parents_ + query.GetNumVertices(), 0);
}

DAG::~DAG() { /* code */ }

void DAG::BuildDAG() { /* code */ }

Vertex DAG::SelectRootVertex() { /* code */ }
}  // namespace daf
