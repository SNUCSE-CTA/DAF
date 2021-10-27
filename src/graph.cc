#include "include/graph.h"

#include <fstream>
#include <iostream>

namespace daf {
Graph::Graph() {}

Graph::~Graph() {
  delete[] start_off_;
  delete[] linear_adj_list_;
  delete[] label_;
  delete[] label_frequency_;
  delete[] core_num_;
}

void Graph::LoadRoughGraph(const std::string &filename,
                           std::vector<std::vector<Vertex>> *graph) {
  std::ifstream fin(filename);

  if (!fin.is_open()) {
    std::cerr << "Graph file " << filename << " not found!\n";
    exit(EXIT_FAILURE);
  }

  Size v, e;
  char type;

  fin >> type >> v >> e;

  num_vertex_ = v;
  num_edge_ = e;
  label_ = new Label[v];

  graph->resize(v);

  // preprocessing
  while (fin >> type) {
    if (type == 'v') {
      Vertex id;
      Label l;
      Size d;
      fin >> id >> l >> d;

      label_[id] = l;
    } else if (type == 'e') {
      Vertex v1, v2;
      fin >> v1 >> v2;

      (*graph)[v1].push_back(v2);
      (*graph)[v2].push_back(v1);
    }
  }

  fin.close();
}

void Graph::computeCoreNum() { /* code */ }
}  // namespace daf
