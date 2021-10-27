#include "include/graph.h"

#include <fstream>
#include <iostream>

namespace daf {
Graph::Graph() {}

Graph::~Graph() {
  delete[] start_off_;
  delete[] linear_adj_list_;
  delete[] label_;
  delete[] label_frequency_;
  delete[] core_num_;
}

void Graph::LoadRoughGraph(const std::string &filename,
                           std::vector<std::vector<Vertex>> *graph) { /* code */ }

void Graph::computeCoreNum() { /* code */ }
}  // namespace daf
