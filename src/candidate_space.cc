#include "include/candidate_space.h"

namespace daf {
CandidateSpace::CandidateSpace(const DataGraph &data, const QueryGraph &query,
                               const DAG &dag)
    : data_(data), query_(query), dag_(dag) { /* code */ }

CandidateSpace::~CandidateSpace() { /* code */ }

bool CandidateSpace::BuildCS() { /* code */ }

bool CandidateSpace::FilterByTopDownWithInit() { /* code */ }

bool CandidateSpace::FilterByBottomUp() { /* code */ }

bool CandidateSpace::FilterByTopDown() { /* code */ }

void CandidateSpace::ConstructCS() { /* code */ }

bool CandidateSpace::InitRootCandidates() { /* code */ }

void CandidateSpace::AllocateSpaceForCS() { /* code */ }

void CandidateSpace::ComputeNbrInformation(Vertex u, Size *max_nbr_degree,
                                           uint64_t *nbr_label_bitset) { /* code */ }
}  // namespace daf
