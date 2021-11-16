#include "include/ordering.h"

namespace daf {
Ordering::Ordering(Size num_query_vertices)
    : num_query_vertices_(num_query_vertices) {
  extendable_queue_ = new Vertex[num_query_vertices_];
  extendable_queue_size_ = 0;
  weights_ = new Size[num_query_vertices_];
  std::fill(weights_, weights_ + num_query_vertices_,
            std::numeric_limits<Size>::max());
}

Ordering::~Ordering() {
  delete[] extendable_queue_;
  delete[] weights_;
}

void Ordering::Insert(Vertex u, Size weight) { /* code */ }

void Ordering::UpdateWeight(Vertex u, Size weight) { /* code */ }

void Ordering::Remove(Vertex u) { /* code */ }

bool Ordering::Exists(Vertex u) { /* code */ }

Vertex Ordering::PopMinWeight() { /* code */ }
}  // namespace daf
