#include "include/dag.h"

#include <algorithm>

namespace daf {
DAG::DAG(const DataGraph &data, const QueryGraph &query)
    : data_(data), query_(query) { /* code */ }

DAG::~DAG() { /* code */ }

void DAG::BuildDAG() { /* code */ }

Vertex DAG::SelectRootVertex() { /* code */ }
}  // namespace daf
