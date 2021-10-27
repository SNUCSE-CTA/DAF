#ifndef DATA_GRAPH_H_
#define DATA_GRAPH_H_

#include <algorithm>
#include <climits>
#include <string>
#include <unordered_map>
#include <utility>

#include "global/global.h"
#include "include/graph.h"

namespace daf {
class DataGraph : public Graph {
 public:
  DataGraph();
  ~DataGraph();

  DataGraph& operator=(const DataGraph&) = delete;
  DataGraph(const DataGraph&) = delete;

  void LoadAndProcessGraph(const std::string& filename);

  using Graph::GetEndOffset;
  using Graph::GetStartOffset;

  inline Size GetStartOffset(Vertex v, Label l) const;
  inline Size GetEndOffset(Vertex v, Label l) const;

  inline Label GetTransferredLabel(Label l) const;

  inline Size GetNbrBitsetSize() const;

  inline Size GetNeighborLabelFrequency(Vertex v, Label l) const;
  inline Size GetMaxLabelFrequency() const;
  inline Size GetStartOffsetByLabel(Label l) const;
  inline Size GetEndOffsetByLabel(Label l) const;
  inline Size GetVertexBySortedLabelOffset(Size i) const;
  inline Size GetMaxNbrDegree(Vertex v) const;

  inline Size GetInitCandSize(Label l, Size d) const;
  inline bool CheckAllNbrLabelExist(Vertex v, uint64_t* nbr_bitset) const;

 private:
  Label* transferred_label_;
  std::pair<Size, Size>* adj_offs_by_label_;

  Size* offs_by_label_;
  Vertex* vertices_sorted_;

  uint64_t* linear_nbr_bitset_;
  Size* max_nbr_degree_;

  Size nbr_bitset_size_;
  Size max_label_frequency_;
};

inline Size DataGraph::GetStartOffset(Vertex v, Label l) const {
  return adj_offs_by_label_[v * GetNumLabels() + l].first;
}

inline Size DataGraph::GetEndOffset(Vertex v, Label l) const { /* code */ }

inline Size DataGraph::GetNeighborLabelFrequency(Vertex v, Label l) const { /* code */ }

inline Label DataGraph::GetTransferredLabel(Label l) const { /* code */ }

inline Size DataGraph::GetStartOffsetByLabel(Label l) const { /* code */ }

inline Size DataGraph::GetEndOffsetByLabel(Label l) const { /* code */ }

inline Size DataGraph::GetVertexBySortedLabelOffset(Size i) const { /* code */ }

inline Size DataGraph::GetNbrBitsetSize() const { /* code */ }

inline Size DataGraph::GetMaxNbrDegree(Vertex v) const { /* code */ }

inline Size DataGraph::GetMaxLabelFrequency() const { /* code */ }

inline Size DataGraph::GetInitCandSize(Label l, Size d) const { /* code */ }

inline bool DataGraph::CheckAllNbrLabelExist(Vertex v,
                                             uint64_t* nbr_bitset) const { /* code */ }
}  // namespace daf
#endif  // DATA_GRAPH_H_
