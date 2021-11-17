#ifndef BACKTRACK_HELPER_H_
#define BACKTRACK_HELPER_H_

#include <boost/dynamic_bitset.hpp>

#include "global/global.h"

namespace daf {
enum MappingState { UNMAPPED, MAPPED, RESERVED };

struct SearchTreeNode {
  bool initialized;
  Vertex u;
  Vertex v;
  Size v_idx;
  boost::dynamic_bitset<> failing_set;
  bool embedding_founded;
};

class BacktrackHelper {
 public:
  inline BacktrackHelper();
  inline ~BacktrackHelper();

  BacktrackHelper &operator=(const BacktrackHelper &) = delete;
  BacktrackHelper(const BacktrackHelper &) = delete;

  inline void Initialize(Size num_query_vtx, Size degree, Size cs_size,
                         Vertex u);

  void AddMapping(Vertex u) {
    num_mapped_neighbors_ += 1;
    num_extendable_indices_[num_mapped_neighbors_] = 0;
    ancestors_[num_mapped_neighbors_] = ancestors_[num_mapped_neighbors_ - 1];
    lastly_mapped_neighbor_[num_mapped_neighbors_] = u;
  }
  void RemoveMapping() { num_mapped_neighbors_ -= 1; }

  inline boost::dynamic_bitset<> &GetAncestor() {
    return ancestors_[num_mapped_neighbors_];
  }
  inline Size &GetNumExtendable() {
    return num_extendable_indices_[num_mapped_neighbors_];
  }
  inline Size *GetExtendableIndices() {
    return extendable_indices_[num_mapped_neighbors_ - 1];
  }
  inline Size GetNumPrevExtendable() {
    return num_extendable_indices_[num_mapped_neighbors_ - 1];
  }
  inline Size *GetPrevExtendableIndices() {
    return extendable_indices_[num_mapped_neighbors_ - 2];
  }
  inline Size GetExtendableIndex(Size i) {
    if (num_mapped_neighbors_ == 0)
      return i;
    else
      return extendable_indices_[num_mapped_neighbors_ - 1][i];
  }
  inline MappingState &GetMappingState() { return mapping_state_; }
  inline Size GetNumMappedNeighbors() { return num_mapped_neighbors_; }
  inline Vertex GetLastlyMappedNeighbor() {
    return lastly_mapped_neighbor_[num_mapped_neighbors_];
  }

 private:
  Size num_mapped_neighbors_;
  boost::dynamic_bitset<> *ancestors_;
  Size **extendable_indices_;
  Size *num_extendable_indices_;
  Vertex *lastly_mapped_neighbor_;
  MappingState mapping_state_;

  Size degree_ = -1;
};
inline BacktrackHelper::BacktrackHelper() { /* code */ }

inline BacktrackHelper::~BacktrackHelper() { /* code */ }

inline void BacktrackHelper::Initialize(Size num_query_vtx, Size degree,
                                        Size cs_size, Vertex u) { /* code */ }

}  // namespace daf

#endif  // BACKTRACK_HELPER_H_
