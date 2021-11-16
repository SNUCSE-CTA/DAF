#include "include/ordering.h"

namespace daf {
Ordering::Ordering(Size num_query_vertices)
    : num_query_vertices_(num_query_vertices) { /* code */ }

Ordering::~Ordering() { /* code */ }

void Ordering::Insert(Vertex u, Size weight) { /* code */ }

void Ordering::UpdateWeight(Vertex u, Size weight) { /* code */ }

void Ordering::Remove(Vertex u) { /* code */ }

bool Ordering::Exists(Vertex u) { /* code */ }

Vertex Ordering::PopMinWeight() { /* code */ }
}  // namespace daf
