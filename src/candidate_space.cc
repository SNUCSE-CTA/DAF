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

void CandidateSpace::ConstructCS() {
  AllocateSpaceForCS();

  for (Size i = 0; i < query_.GetNumVertices(); ++i) {
    Vertex u = dag_.GetVertexOrderedByBFS(i);

    if (query_.IsInNEC(u) && !query_.IsNECRepresentation(u)) continue;

    Label u_label = query_.GetLabel(u);
    Size u_degree = query_.GetDegree(u);

    for (Size i = 0; i < candidate_set_size_[u]; ++i) {
      cand_to_cs_idx_[candidate_set_[u][i]] = i;
    }

    for (Size i = query_.GetStartOffset(u); i < query_.GetEndOffset(u); ++i) {
      Vertex u_adj = query_.GetNeighbor(i);

      if (query_.IsInNEC(u_adj) && !query_.IsNECRepresentation(u_adj)) continue;

      Size u_adj_idx = i - query_.GetStartOffset(u);

      Size start_offset =
          candidate_offsets_[u * query_.GetMaxDegree() + u_adj_idx][0];

      for (Size v_adj_idx = 0; v_adj_idx < candidate_set_size_[u_adj];
           ++v_adj_idx) {
        Vertex v_adj = candidate_set_[u_adj][v_adj_idx];

        for (Size i = data_.GetStartOffset(v_adj, u_label);
             i < data_.GetEndOffset(v_adj, u_label); ++i) {
          Vertex v = data_.GetNeighbor(i);

          if (data_.GetDegree(v) < u_degree) break;

          Size v_idx = cand_to_cs_idx_[v];

          if (v_idx != -1) {
            linear_cs_adj_list_[candidate_offsets_[u * query_.GetMaxDegree() +
                                                   u_adj_idx][v_idx]] =
                v_adj_idx;
            candidate_offsets_[u * query_.GetMaxDegree() + u_adj_idx][v_idx] +=
                1;
          }
        }
      }

      for (Size i = GetCandidateSetSize(u) - 1; i > 0; --i) {
        candidate_offsets_[u * query_.GetMaxDegree() + u_adj_idx][i] =
            candidate_offsets_[u * query_.GetMaxDegree() + u_adj_idx][i - 1];
      }

      candidate_offsets_[u * query_.GetMaxDegree() + u_adj_idx][0] =
          start_offset;
    }

    for (Size i = 0; i < candidate_set_size_[u]; ++i) {
      cand_to_cs_idx_[candidate_set_[u][i]] = -1;
    }
  }
}

bool CandidateSpace::InitRootCandidates() {
  Vertex root = dag_.GetRoot();
  Label root_label = query_.GetLabel(root);

  uint64_t *nbr_label_bitset = new uint64_t[data_.GetNbrBitsetSize()];
  Size max_nbr_degree;

  ComputeNbrInformation(root, &max_nbr_degree, nbr_label_bitset);

  for (Size i = data_.GetStartOffsetByLabel(root_label);
       i < data_.GetEndOffsetByLabel(root_label); ++i) {
    Vertex cand = data_.GetVertexBySortedLabelOffset(i);

    if (data_.GetDegree(cand) < query_.GetDegree(root)) break;

    if (data_.GetCoreNum(cand) >= query_.GetCoreNum(root) &&
        data_.CheckAllNbrLabelExist(cand, nbr_label_bitset) &&
        data_.GetMaxNbrDegree(cand) >= max_nbr_degree) {
      candidate_set_[root][candidate_set_size_[root]] = cand;
      candidate_set_size_[root] += 1;
    }
  }

  delete[] nbr_label_bitset;
  return candidate_set_size_[root] > 0;
}

void CandidateSpace::AllocateSpaceForCS() {
  Size cur_cand_offset = 0;

  for (Size i = 0; i < query_.GetNumVertices(); ++i) {
    Vertex u = dag_.GetVertexOrderedByBFS(i);

    if (query_.IsInNEC(u) && !query_.IsNECRepresentation(u)) continue;

    Label u_label = query_.GetLabel(u);
    Size u_degree = query_.GetDegree(u);

    for (Size i = 0; i < candidate_set_size_[u]; ++i) {
      cand_to_cs_idx_[candidate_set_[u][i]] = i;
    }

    for (Size i = query_.GetStartOffset(u); i < query_.GetEndOffset(u); ++i) {
      Vertex u_adj = query_.GetNeighbor(i);

      if (query_.IsInNEC(u_adj) && !query_.IsNECRepresentation(u_adj)) continue;

      Size u_adj_idx = i - query_.GetStartOffset(u);

      candidate_offsets_[u * query_.GetMaxDegree() + u_adj_idx] =
          new Size[GetCandidateSetSize(u) + 1];

      std::fill(candidate_offsets_[u * query_.GetMaxDegree() + u_adj_idx],
                candidate_offsets_[u * query_.GetMaxDegree() + u_adj_idx] +
                    GetCandidateSetSize(u) + 1,
                cur_cand_offset);

      for (Size v_adj_idx = 0; v_adj_idx < candidate_set_size_[u_adj];
           ++v_adj_idx) {
        Vertex v_adj = candidate_set_[u_adj][v_adj_idx];

        for (Size i = data_.GetStartOffset(v_adj, u_label);
             i < data_.GetEndOffset(v_adj, u_label); ++i) {
          Vertex v = data_.GetNeighbor(i);
          if (data_.GetDegree(v) < u_degree) break;

          if (cand_to_cs_idx_[v] != -1) {
            candidate_offsets_[u * query_.GetMaxDegree() + u_adj_idx]
                              [cand_to_cs_idx_[v] + 1] += 1;
          }
        }
      }

      for (Size i = 2; i < GetCandidateSetSize(u) + 1; ++i) {
        candidate_offsets_[u * query_.GetMaxDegree() + u_adj_idx][i] +=
            candidate_offsets_[u * query_.GetMaxDegree() + u_adj_idx][i - 1] -
            cur_cand_offset;
      }

      cur_cand_offset = candidate_offsets_[u * query_.GetMaxDegree() +
                                           u_adj_idx][GetCandidateSetSize(u)];
    }

    for (Size i = 0; i < candidate_set_size_[u]; ++i) {
      cand_to_cs_idx_[candidate_set_[u][i]] = -1;
    }
  }
  linear_cs_adj_list_ = new Vertex[cur_cand_offset];
}

void CandidateSpace::ComputeNbrInformation(Vertex u, Size *max_nbr_degree,
                                           uint64_t *nbr_label_bitset) { /* code */ }
}  // namespace daf
