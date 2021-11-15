#include "include/candidate_space.h"

namespace daf {
CandidateSpace::CandidateSpace(const DataGraph &data, const QueryGraph &query,
                               const DAG &dag)
    : data_(data), query_(query), dag_(dag) {
  candidate_set_size_ = new Size[query_.GetNumVertices()];
  candidate_set_ = new Vertex *[query_.GetNumVertices()];

  for (Vertex v = 0; v < query_.GetNumVertices(); ++v) {
    if (query_.IsInNEC(v) && !query_.IsNECRepresentation(v))
      candidate_set_[v] = nullptr;
    else
      candidate_set_[v] = new Vertex[dag_.GetInitCandSize(v)];
  }

  candidate_offsets_ =
      new Size *[query_.GetNumVertices() * query_.GetMaxDegree()];
  linear_cs_adj_list_ = nullptr;

  num_visit_cs_ = new Size[data_.GetNumVertices()];
  visited_candidates_ =
      new Vertex[data_.GetNumVertices() * query_.GetNumVertices()];
  cand_to_cs_idx_ = new Size[data_.GetNumVertices()];

  num_visitied_candidates_ = 0;

  std::fill(num_visit_cs_, num_visit_cs_ + data_.GetNumVertices(), 0);
  std::fill(candidate_set_size_, candidate_set_size_ + query_.GetNumVertices(),
            0);
  std::fill(
      candidate_offsets_,
      candidate_offsets_ + query_.GetNumVertices() * query_.GetMaxDegree(),
      nullptr);
  std::fill(cand_to_cs_idx_, cand_to_cs_idx_ + data_.GetNumVertices(), -1);
}

CandidateSpace::~CandidateSpace() {
  delete[] candidate_set_size_;

  for (Vertex v = 0; v < query_.GetNumVertices(); ++v) {
    if (candidate_set_[v] != nullptr) delete[] candidate_set_[v];
  }
  delete[] candidate_set_;

  for (Size i = 0; i < query_.GetNumVertices() * query_.GetMaxDegree(); ++i) {
    if (candidate_offsets_[i]) delete[] candidate_offsets_[i];
  }
  delete[] candidate_offsets_;
  if (linear_cs_adj_list_ != nullptr) delete[] linear_cs_adj_list_;

  delete[] num_visit_cs_;
  delete[] visited_candidates_;
  delete[] cand_to_cs_idx_;
}

bool CandidateSpace::BuildCS() {
  if (!FilterByTopDownWithInit()) return false;
  if (!FilterByBottomUp()) return false;
  if (!FilterByTopDown()) return false;

  ConstructCS();

  return true;
}

bool CandidateSpace::FilterByTopDownWithInit() {
  bool result = true;

  uint64_t *nbr_label_bitset = new uint64_t[data_.GetNbrBitsetSize()];
  Size max_nbr_degree;

  if (!InitRootCandidates()) {
    result = false;
  } else {
    for (Size i = 1; i < query_.GetNumVertices(); ++i) {
      Vertex cur = dag_.GetVertexOrderedByBFS(i);

      if (query_.IsInNEC(cur) && !query_.IsNECRepresentation(cur)) continue;

      Label cur_label = query_.GetLabel(cur);
      Size num_parent = 0;
      for (Size i = 0; i < dag_.GetNumParents(cur); ++i) {
        Vertex parent = dag_.GetParent(cur, i);

        for (Size i = 0; i < candidate_set_size_[parent]; ++i) {
          Vertex parent_cand = candidate_set_[parent][i];

          for (Size i = data_.GetStartOffset(parent_cand, cur_label);
               i < data_.GetEndOffset(parent_cand, cur_label); ++i) {
            Vertex cand = data_.GetNeighbor(i);

            if (data_.GetDegree(cand) < query_.GetDegree(cur)) break;

            if (num_visit_cs_[cand] == num_parent) {
              num_visit_cs_[cand] += 1;
              if (num_parent == 0) {
                visited_candidates_[num_visitied_candidates_] = cand;
                num_visitied_candidates_ += 1;
              }
            }
          }
        }
        num_parent += 1;
      }

      ComputeNbrInformation(cur, &max_nbr_degree, nbr_label_bitset);

      for (Size i = 0; i < num_visitied_candidates_; ++i) {
        Vertex cand = visited_candidates_[i];
        if (num_visit_cs_[cand] == num_parent &&
            data_.GetCoreNum(cand) >= query_.GetCoreNum(cur) &&
            data_.CheckAllNbrLabelExist(cand, nbr_label_bitset) &&
            data_.GetMaxNbrDegree(cand) >= max_nbr_degree) {
          candidate_set_[cur][candidate_set_size_[cur]] = cand;
          candidate_set_size_[cur] += 1;
        }
      }

      if (candidate_set_size_[cur] == 0) {
        result = false;
        break;
      }

      while (num_visitied_candidates_ > 0) {
        num_visitied_candidates_ -= 1;
        num_visit_cs_[visited_candidates_[num_visitied_candidates_]] = 0;
      }
    }
  }

  delete[] nbr_label_bitset;
  return result;
}

bool CandidateSpace::FilterByBottomUp() {
  bool result = true;
  for (Size i = 0; i < query_.GetNumVertices(); ++i) {
    Vertex cur = dag_.GetVertexOrderedByBFS(query_.GetNumVertices() - i - 1);

    if (dag_.GetNumChildren(cur) == 0) continue;

    Label cur_label = query_.GetLabel(cur);

    Size num_child = 0;
    for (Size i = 0; i < dag_.GetNumChildren(cur); ++i) {
      Vertex child = dag_.GetChild(cur, i);

      if (query_.IsInNEC(child) && !query_.IsNECRepresentation(child)) continue;

      for (Size i = 0; i < candidate_set_size_[child]; ++i) {
        Vertex child_cand = candidate_set_[child][i];

        for (Size i = data_.GetStartOffset(child_cand, cur_label);
             i < data_.GetEndOffset(child_cand, cur_label); ++i) {
          Vertex cand = data_.GetNeighbor(i);
          if (data_.GetDegree(cand) < query_.GetDegree(cur)) break;

          if (num_visit_cs_[cand] == num_child) {
            num_visit_cs_[cand] += 1;
            if (num_child == 0) {
              visited_candidates_[num_visitied_candidates_] = cand;
              num_visitied_candidates_ += 1;
            }
          }
        }
      }
      num_child += 1;
    }

    for (Size i = 0; i < candidate_set_size_[cur]; ++i) {
      Vertex cand = candidate_set_[cur][i];
      if (num_visit_cs_[cand] != num_child) {
        std::swap(candidate_set_[cur][i],
                  candidate_set_[cur][candidate_set_size_[cur] - 1]);
        candidate_set_size_[cur] -= 1;
        --i;
      }
    }

    if (candidate_set_size_[cur] == 0) {
      result = false;
      break;
    }

    while (num_visitied_candidates_ > 0) {
      num_visitied_candidates_ -= 1;
      num_visit_cs_[visited_candidates_[num_visitied_candidates_]] = 0;
    }
  }

  return result;
}

bool CandidateSpace::FilterByTopDown() {
  bool result = true;
  for (Size i = 1; i < query_.GetNumVertices(); ++i) {
    Vertex cur = dag_.GetVertexOrderedByBFS(i);
    if (query_.IsInNEC(cur) && !query_.IsNECRepresentation(cur)) continue;

    Label cur_label = query_.GetLabel(cur);

    Size num_parent = 0;
    for (Size i = 0; i < dag_.GetNumParents(cur); ++i) {
      Vertex parent = dag_.GetParent(cur, i);

      for (Size i = 0; i < candidate_set_size_[parent]; ++i) {
        Vertex parent_cand = candidate_set_[parent][i];

        for (Size i = data_.GetStartOffset(parent_cand, cur_label);
             i < data_.GetEndOffset(parent_cand, cur_label); ++i) {
          Vertex cand = data_.GetNeighbor(i);
          if (data_.GetDegree(cand) < query_.GetDegree(cur)) break;

          if (num_visit_cs_[cand] == num_parent) {
            num_visit_cs_[cand] += 1;
            if (num_parent == 0) {
              visited_candidates_[num_visitied_candidates_] = cand;
              num_visitied_candidates_ += 1;
            }
          }
        }
      }
      num_parent += 1;
    }

    for (Size i = 0; i < candidate_set_size_[cur]; ++i) {
      Vertex cand = candidate_set_[cur][i];
      if (num_visit_cs_[cand] != num_parent) {
        std::swap(candidate_set_[cur][i],
                  candidate_set_[cur][candidate_set_size_[cur] - 1]);
        candidate_set_size_[cur] -= 1;
        --i;
      }
    }

    if (candidate_set_size_[cur] == 0) {
      result = false;
      break;
    }

    while (num_visitied_candidates_ > 0) {
      num_visitied_candidates_ -= 1;
      num_visit_cs_[visited_candidates_[num_visitied_candidates_]] = 0;
    }
  }

  return result;
}

void CandidateSpace::ConstructCS() { /* code */ }

bool CandidateSpace::InitRootCandidates() { /* code */ }

void CandidateSpace::AllocateSpaceForCS() { /* code */ }

void CandidateSpace::ComputeNbrInformation(Vertex u, Size *max_nbr_degree,
                                           uint64_t *nbr_label_bitset) { /* code */ }
}  // namespace daf
