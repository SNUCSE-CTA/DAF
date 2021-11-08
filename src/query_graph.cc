#include "include/query_graph.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

namespace daf {
QueryGraph::QueryGraph() { /* code */ }
QueryGraph::~QueryGraph() { /* code */ }

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
