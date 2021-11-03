#include "include/data_graph.h"

#include <algorithm>
#include <unordered_map>
#include <vector>

namespace daf {
DataGraph::DataGraph() {}

DataGraph::~DataGraph() {
  delete[] transferred_label_;
  delete[] adj_offs_by_label_;
  delete[] offs_by_label_;
  delete[] vertices_sorted_;
  delete[] linear_nbr_bitset_;
  delete[] max_nbr_degree_;
}

void DataGraph::LoadAndProcessGraph(const std::string& filename) { /* code */ }
}  // namespace daf
