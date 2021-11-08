#include "include/query_graph.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

namespace daf {
QueryGraph::QueryGraph() {}
QueryGraph::~QueryGraph() {
  delete[] NEC_map_;
  delete[] NEC_elems_;

  if (NEC_start_offs_) delete[] NEC_start_offs_;
  delete[] NEC_size_;
}

bool QueryGraph::LoadAndProcessGraph(const std::string &filename,
                                     const DataGraph &data) { /* code */ }

namespace {
struct NECInfo {
  bool visit = false;
  Vertex representation;
  Size NEC_elems_idx;
};
}  // namespace

void QueryGraph::ExtractResidualStructure() { /* code */ }
}  // namespace daf
