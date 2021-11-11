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

bool CandidateSpace::FilterByTopDownWithInit() { /* code */ }

bool CandidateSpace::FilterByBottomUp() { /* code */ }

bool CandidateSpace::FilterByTopDown() { /* code */ }

void CandidateSpace::ConstructCS() { /* code */ }

bool CandidateSpace::InitRootCandidates() { /* code */ }

void CandidateSpace::AllocateSpaceForCS() { /* code */ }

void CandidateSpace::ComputeNbrInformation(Vertex u, Size *max_nbr_degree,
                                           uint64_t *nbr_label_bitset) { /* code */ }
}  // namespace daf
