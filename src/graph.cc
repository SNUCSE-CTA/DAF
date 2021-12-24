#include "include/graph.h"

namespace daf {
Graph::Graph(const std::string &filename)
    : filename_(filename), fin_(filename) {}

Graph::~Graph() {
  delete[] start_off_;
  delete[] linear_adj_list_;
  delete[] label_;
  delete[] label_frequency_;
  delete[] core_num_;
}

void Graph::LoadRoughGraph(std::vector<std::vector<Vertex>> *graph) {
  if (!fin_.is_open()) {
    std::cerr << "Graph file " << filename_ << " not found!\n";
    exit(EXIT_FAILURE);
  }

  Size v, e;
  char type;

  fin_ >> type >> v >> e;

  num_vertex_ = v;
  num_edge_ = e;
  label_ = new Label[v];

  graph->resize(v);

  // preprocessing
  while (fin_ >> type) {
    if (type == 'v') {
      Vertex id;
      Label l;
      fin_ >> id >> l;

      label_[id] = l;
    } else if (type == 'e') {
      Vertex v1, v2;
      fin_ >> v1 >> v2;

      (*graph)[v1].push_back(v2);
      (*graph)[v2].push_back(v1);
    }
  }

  fin_.close();
}

void Graph::computeCoreNum() {
  Size *bin = new Size[max_degree_ + 1];
  Size *pos = new Size[GetNumVertices()];
  Vertex *vert = new Vertex[GetNumVertices()];

  std::fill(bin, bin + (max_degree_ + 1), 0);

  for (Vertex v = 0; v < GetNumVertices(); ++v) {
    bin[core_num_[v]] += 1;
  }

  Size start = 0;
  Size num;

  for (Size d = 0; d <= max_degree_; ++d) {
    num = bin[d];
    bin[d] = start;
    start += num;
  }

  for (Vertex v = 0; v < GetNumVertices(); ++v) {
    pos[v] = bin[core_num_[v]];
    vert[pos[v]] = v;
    bin[core_num_[v]] += 1;
  }

  for (Size d = max_degree_; d--;) bin[d + 1] = bin[d];
  bin[0] = 0;

  for (Size i = 0; i < GetNumVertices(); ++i) {
    Vertex v = vert[i];

    for (Size j = GetStartOffset(v); j < GetEndOffset(v); j++) {
      Vertex u = GetNeighbor(j);

      if (core_num_[u] > core_num_[v]) {
        Size du = core_num_[u];
        Size pu = pos[u];

        Size pw = bin[du];
        Vertex w = vert[pw];

        if (u != w) {
          pos[u] = pw;
          pos[w] = pu;
          vert[pu] = w;
          vert[pw] = u;
        }

        bin[du]++;
        core_num_[u]--;
      }
    }
  }

  delete[] bin;
  delete[] pos;
  delete[] vert;
}
}  // namespace daf
