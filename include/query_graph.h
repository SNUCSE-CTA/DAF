#ifndef QUERY_GRAPH_H_
#define QUERY_GRAPH_H_

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "include/data_graph.h"
#include "include/graph.h"

namespace daf {
struct NECElement;

class QueryGraph : public Graph {
 public:
  explicit QueryGraph(const std::string &filename);
  ~QueryGraph();

  QueryGraph &operator=(const QueryGraph &) = delete;
  QueryGraph(const QueryGraph &) = delete;

  bool LoadAndProcessGraph(const DataGraph &data);

  inline bool IsNECRepresentation(Vertex v) const;
  inline bool IsInNEC(Vertex v) const;

  inline Size GetNumNECLabel() const;
  inline Size GetNECStartOffset(Size i) const;
  inline Size GetNECEndOffset(Size i) const;
  inline const NECElement &GetNECElement(Size i) const;
  inline bool IsTree() const;

  inline Size GetNumNonLeafVertices() const;
  inline Size GetNECSize(Vertex v) const;
  inline Size GetNECRepresentative(Vertex v) const;

  inline Label GetMaxLabel() const;

 private:
  void ExtractResidualStructure();

  Vertex *NEC_map_;
  NECElement *NEC_elems_;
  Size *NEC_start_offs_;
  Size *NEC_size_;

  Size num_NEC_label_;
  Size num_non_leaf_vertices_;
  Label max_label_;

  bool is_tree_;
};

struct NECElement {
  Label label;
  Vertex adjacent;
  Vertex represent;
  Size size;
  Size represent_idx;
};

inline bool QueryGraph::IsNECRepresentation(Vertex v) const {
  return NEC_map_[v] == v;
}

inline bool QueryGraph::IsInNEC(Vertex v) const {
  return NEC_map_[v] != INVALID_VTX;
}

inline Size QueryGraph::GetNumNECLabel() const { return num_NEC_label_; }

inline Size QueryGraph::GetNECStartOffset(Size i) const {
  return NEC_start_offs_[i];
}

inline Size QueryGraph::GetNECEndOffset(Size i) const {
  return NEC_start_offs_[i + 1];
}

inline const NECElement &QueryGraph::GetNECElement(Size i) const {
  return NEC_elems_[i];
}

inline bool QueryGraph::IsTree() const { return is_tree_; }

inline Size QueryGraph::GetNumNonLeafVertices() const {
  return num_non_leaf_vertices_;
}

inline Size QueryGraph::GetNECSize(Vertex v) const { return NEC_size_[v]; }

inline Size QueryGraph::GetNECRepresentative(Vertex v) const {
  if (IsInNEC(v)) {
    return NEC_map_[v];
  } else {
    return v;
  }
}

inline Label QueryGraph::GetMaxLabel() const { return max_label_; }
}  // namespace daf

#endif  // QUERY_GRAPH_H_
