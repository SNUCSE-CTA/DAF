#ifndef GRAPH_H_
#define GRAPH_H_

#include <memory>
#include <string>
#include <vector>

#include "global/global.h"

namespace daf {

class Graph {
 public:
  Graph();
  ~Graph();

  Graph& operator=(const Graph&) = delete;
  Graph(const Graph&) = delete;

  inline Size GetNumLabels() const;
  inline Size GetNumVertices() const;
  inline Size GetNumEdges() const;
  inline Size GetMaxDegree() const;

  inline Label GetLabel(Vertex v) const;
  inline Size GetStartOffset(Vertex v) const;
  inline Size GetEndOffset(Vertex v) const;
  inline Size GetDegree(Vertex v) const;
  inline Size GetCoreNum(Vertex v) const;

  inline Label GetLabelFrequency(Label l) const;

  inline Vertex GetNeighbor(Size i) const;

 protected:
  void LoadRoughGraph(const std::string& filename,
                      std::vector<std::vector<Vertex>>* graph);
  void computeCoreNum();

  Size num_vertex_;
  Size num_edge_;
  Size num_label_;

  Size max_degree_;

  Label* label_;
  Size* start_off_;
  Vertex* linear_adj_list_;
  Size* label_frequency_;

  Size* core_num_;
};

inline Size Graph::GetNumLabels() const { return num_label_; }

inline Size Graph::GetNumVertices() const { return num_vertex_; }

inline Size Graph::GetNumEdges() const { return num_edge_; }

inline Size Graph::GetMaxDegree() const { return max_degree_; }

inline Label Graph::GetLabel(Vertex v) const { /* code */ }

inline Size Graph::GetStartOffset(Vertex v) const { /* code */ }

inline Size Graph::GetEndOffset(Vertex v) const { /* code */ }

inline Size Graph::GetDegree(Vertex v) const { /* code */ }

inline Size Graph::GetCoreNum(Vertex v) const { /* code */ }

inline Label Graph::GetLabelFrequency(Label l) const { /* code */ }

inline Vertex Graph::GetNeighbor(Size i) const { /* code */ }

}  // namespace daf

#endif  // GRAPH_H_
