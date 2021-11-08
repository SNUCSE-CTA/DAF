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
                                     const DataGraph &data) {
  std::vector<std::vector<Vertex>> adj_list;

  LoadRoughGraph(filename, &adj_list);

  max_degree_ = 0;
  num_label_ = 0;
  max_label_ = 0;

  label_frequency_ = new Size[data.GetNumLabels()];
  start_off_ = new Size[GetNumVertices() + 1];
  linear_adj_list_ = new Size[GetNumEdges() * 2];
  core_num_ = new Size[GetNumVertices()];
  std::fill(label_frequency_, label_frequency_ + data.GetNumLabels(), 0);

  Size cur_idx = 0;

  // transfer label & construct adj list and label frequency
  for (Vertex v = 0; v < GetNumVertices(); ++v) {
    Label l = data.GetTransferredLabel(label_[v]);
    if (l < 0) return false;
    label_[v] = l;
    if (l > max_label_) max_label_ = l;
    if (label_frequency_[l] == 0) num_label_ += 1;
    label_frequency_[l] += 1;
    if (static_cast<Size>(adj_list[v].size()) > max_degree_)
      max_degree_ = adj_list[v].size();

    start_off_[v] = cur_idx;
    core_num_[v] = adj_list[v].size();

    std::copy(adj_list[v].begin(), adj_list[v].end(),
              linear_adj_list_ + cur_idx);

    cur_idx += adj_list[v].size();
  }
  start_off_[GetNumVertices()] = num_edge_ * 2;

  // preprocess for query graph
  computeCoreNum();

  is_tree_ = true;
  for (Vertex v = 0; v < GetNumVertices(); ++v) {
    if (GetCoreNum(v) > 1) {
      is_tree_ = false;
      break;
    }
  }

  ExtractResidualStructure();

  return true;
}

namespace {
struct NECInfo {
  bool visit = false;
  Vertex representation;
  Size NEC_elems_idx;
};
}  // namespace

void QueryGraph::ExtractResidualStructure() { /* code */ }
}  // namespace daf
