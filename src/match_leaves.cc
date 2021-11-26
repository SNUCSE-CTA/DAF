#include "include/match_leaves.h"

namespace daf {
namespace {
constexpr uint64_t Factorial(Size n) {
  return n == 1 ? 1 : n * Factorial(n - 1);
}
}  // namespace

MatchLeaves::MatchLeaves(const DataGraph &data, const QueryGraph &query,
                         const CandidateSpace &cs, Vertex *mapped_query_vtx,
                         BacktrackHelper *helpers,
                         SearchTreeNode **mapped_nodes)
    : data_(data),
      query_(query),
      cs_(cs),
      backtrack_mapped_query_vtx(mapped_query_vtx),
      backtrack_helpers_(helpers),
      backtrack_mapped_nodes_(mapped_nodes) {
  nec_distinct_cands_ = new Vertex[data_.GetNumVertices()];
  cand_to_nec_idx_ = new std::vector<Size>[data_.GetNumVertices()];
  nec_ranking_ = new Size[query_.GetNumNECLabel()];

  sum_nec_cands_size_ = new Size[query_.GetNumNECLabel()];
  sum_nec_size_ = new Size[query_.GetNumNECLabel()];

  num_remaining_cands_ = new Size[query_.GetNumVertices()];
  num_remaining_nec_vertices_ = new Size[query_.GetNumVertices()];

  num_nec_distinct_cands_ = 0;

  maximum_matching_ =
      new MaximumMatching(data_, query_, cs_, backtrack_helpers_);
}

MatchLeaves::~MatchLeaves() {
  delete[] cand_to_nec_idx_;

  delete[] nec_distinct_cands_;
  delete[] nec_ranking_;
  delete[] sum_nec_cands_size_;
  delete[] sum_nec_size_;

  delete[] num_remaining_cands_;
  delete[] num_remaining_nec_vertices_;

  delete maximum_matching_;
}

uint64_t MatchLeaves::Match(uint64_t limit) {
  uint64_t result = 1;

  // first scan
  for (Size i = 0; i < query_.GetNumNECLabel(); ++i) {
    Size start_offset = query_.GetNECStartOffset(i);
    Size end_offset = query_.GetNECEndOffset(i);

    for (Size j = start_offset; j < end_offset; ++j) {
      const NECElement &nec = query_.GetNECElement(j);

      Vertex represent = nec.represent;
      Size size = nec.size;

      Vertex repr_nbr = nec.adjacent;

      BacktrackHelper *repr_helper = backtrack_helpers_ + represent;
      BacktrackHelper *repr_nbr_helper = backtrack_helpers_ + repr_nbr;

      repr_helper->AddMapping(repr_nbr);

      Size *extendable_indices = repr_helper->GetExtendableIndices();
      Size &num_extendable = repr_helper->GetNumExtendable();

      SearchTreeNode *nbr_node = backtrack_mapped_nodes_[repr_nbr];

      Size cs_idx = repr_nbr_helper->GetExtendableIndex(nbr_node->v_idx);

      for (Size i =
               cs_.GetCandidateStartOffset(repr_nbr, nec.represent_idx, cs_idx);
           i < cs_.GetCandidateEndOffset(repr_nbr, nec.represent_idx, cs_idx);
           ++i) {
        Size v_idx = cs_.GetCandidateIndex(i);
        Vertex v = cs_.GetCandidate(represent, v_idx);
        if (backtrack_mapped_query_vtx[v] == -1) {
          extendable_indices[num_extendable] = v_idx;
          num_extendable += 1;
        }
      }

      if (num_extendable < size) {
        result = 0;
        goto EXIT;
      } else {
        if (end_offset - start_offset == 1) {
          for (Size j = 0; j < nec.size; ++j) {
            result *= num_extendable - j;
          }
        } else if (num_extendable == size) {
          ReserveVertex(represent, repr_helper);
          result *= Factorial(size);
        }
      }
    }
  }

  // second scan : perform maximum matching
  for (Size i = 0; i < query_.GetNumNECLabel(); ++i) {
    Size start_offset = query_.GetNECStartOffset(i);
    Size end_offset = query_.GetNECEndOffset(i);

    sum_nec_cands_size_[i] = 0;
    sum_nec_size_[i] = 0;

    nec_ranking_[i] = i;

    if (end_offset - start_offset == 1) continue;

    Label label = query_.GetNECElement(start_offset).label;

    for (Size j = start_offset; j < end_offset; ++j) {
      const NECElement &nec = query_.GetNECElement(j);

      // assert(label == nec.label);

      Vertex represent = nec.represent;
      Size size = nec.size;

      BacktrackHelper *repr_helper = backtrack_helpers_ + represent;

      if (repr_helper->GetMappingState() == RESERVED) continue;

      for (Size k = 0; k < size; ++k) {
        maximum_matching_->AddToU(j);
      }

      sum_nec_size_[i] += size;

      for (Size k = 0; k < repr_helper->GetNumExtendable(); ++k) {
        Vertex cand =
            cs_.GetCandidate(represent, repr_helper->GetExtendableIndex(k));

        // assert(cand >= 0 && cand < data_.GetNumVertices());

        if (backtrack_mapped_query_vtx[cand] == -1) {
          if (maximum_matching_->IsAddedToV(cand) == false) {
            maximum_matching_->AddToV(cand);
            nec_distinct_cands_[num_nec_distinct_cands_] = cand;

            num_nec_distinct_cands_ += 1;
          }
        } else {
          std::swap(
              repr_helper->GetExtendableIndices()[k],
              repr_helper
                  ->GetExtendableIndices()[repr_helper->GetNumExtendable() -
                                           1]);
          repr_helper->GetNumExtendable() -= 1;
          k -= 1;
        }
      }

      sum_nec_cands_size_[i] += repr_helper->GetNumExtendable();

      if (repr_helper->GetNumExtendable() < size) {
        maximum_matching_->Clear(nec_distinct_cands_, &num_nec_distinct_cands_);
        result = 0;
        goto EXIT;
      } else if (num_nec_distinct_cands_ < sum_nec_size_[i]) {
        maximum_matching_->Clear(nec_distinct_cands_, &num_nec_distinct_cands_);
        result = 0;
        goto EXIT;
      } else if (repr_helper->GetNumExtendable() == size) {
        ReserveVertex(represent, repr_helper);
        result *= Factorial(size);
      }
    }

    Size size_mm = maximum_matching_->ComputeMaximumMatching(label);
    maximum_matching_->Clear(nec_distinct_cands_, &num_nec_distinct_cands_);

    if (size_mm < sum_nec_size_[i]) {
      result = 0;
      goto EXIT;
    }
  }

  if (result >= limit) goto EXIT;

  // count maximum matching

  std::sort(nec_ranking_, nec_ranking_ + query_.GetNumNECLabel(),
            [this](const Size &a, const Size &b) -> bool {
              if (this->sum_nec_size_[a] == this->sum_nec_size_[b]) {
                return this->sum_nec_cands_size_[a] <
                       this->sum_nec_cands_size_[b];
              } else {
                return this->sum_nec_size_[a] < this->sum_nec_size_[b];
              }
            });

  for (Size i = 0; i < query_.GetNumNECLabel(); ++i) {
    cur_label_idx = nec_ranking_[i];
    Size start_offset = query_.GetNECStartOffset(cur_label_idx);
    Size end_offset = query_.GetNECEndOffset(cur_label_idx);

    if (end_offset - start_offset == 1) continue;

    bool nec_cand_size_reduced;
    do {
      nec_cand_size_reduced = false;
      for (Size j = start_offset; j < end_offset; ++j) {
        const NECElement &nec = query_.GetNECElement(j);

        Vertex represent = nec.represent;
        Size size = nec.size;

        BacktrackHelper *repr_helper = backtrack_helpers_ + represent;

        if (repr_helper->GetMappingState() == RESERVED) continue;

        // assert(repr_helper->GetNumExtendable() >= nec.size);

        for (Size k = 0; k < repr_helper->GetNumExtendable(); ++k) {
          Vertex cand =
              cs_.GetCandidate(represent, repr_helper->GetExtendableIndex(k));

          if (backtrack_mapped_query_vtx[cand] != -1) {
            std::swap(
                repr_helper->GetExtendableIndices()[k],
                repr_helper
                    ->GetExtendableIndices()[repr_helper->GetNumExtendable() -
                                             1]);
            repr_helper->GetNumExtendable() -= 1;
            k -= 1;
          }
        }

        if (repr_helper->GetNumExtendable() == size) {
          ReserveVertex(represent, repr_helper);
          nec_cand_size_reduced = true;
          result *= Factorial(size);
        }
      }
    } while (nec_cand_size_reduced);

    for (Size j = start_offset; j < end_offset; ++j) {
      const NECElement &nec = query_.GetNECElement(j);

      Vertex represent = nec.represent;
      Size size = nec.size;

      BacktrackHelper *repr_helper = backtrack_helpers_ + represent;

      if (repr_helper->GetMappingState() == RESERVED) {
        num_remaining_nec_vertices_[j] = 0;
        continue;
      }

      num_remaining_cands_[j] = repr_helper->GetNumExtendable();
      num_remaining_nec_vertices_[j] = size;

      for (Size k = 0; k < repr_helper->GetNumExtendable(); ++k) {
        Vertex cand =
            cs_.GetCandidate(represent, repr_helper->GetExtendableIndex(k));

        if (cand_to_nec_idx_[cand].size() == 0) {
          nec_distinct_cands_[num_nec_distinct_cands_] = cand;
          num_nec_distinct_cands_ += 1;
        }
        cand_to_nec_idx_[cand].push_back(j);
      }
    }

    result *= Combine();
    // assert(result != 0);

    for (Size k = 0; k < num_nec_distinct_cands_; ++k) {
      Vertex cand = nec_distinct_cands_[k];

      cand_to_nec_idx_[cand].clear();
    }

    num_nec_distinct_cands_ = 0;

    if (result >= limit) goto EXIT;
  }

EXIT:
  ClearMemoryForBacktrack();

  return result;
}

uint64_t MatchLeaves::Combine() { /* code */ }

void MatchLeaves::ReserveVertex(Vertex represent,
                                BacktrackHelper *repr_helper) { /* code */ }

void MatchLeaves::ClearMemoryForBacktrack() { /* code */ }
}  // namespace daf
