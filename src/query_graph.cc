#include "include/query_graph.h"

namespace daf {
QueryGraph::QueryGraph(const std::string &filename) : Graph(filename) {}
QueryGraph::~QueryGraph() {
  delete[] NEC_map_;
  delete[] NEC_elems_;

  if (NEC_start_offs_) delete[] NEC_start_offs_;
  delete[] NEC_size_;
}

bool QueryGraph::LoadAndProcessGraph(const DataGraph &data) {
  std::vector<std::vector<Vertex>> adj_list;

  LoadRoughGraph(&adj_list);

  max_degree_ = 0;
  num_label_ = 0;
  max_label_ = 0;

  label_frequency_ = new Size[data.GetNumLabels()];
  start_off_ = new Size[GetNumVertices() + 1];
  linear_adj_list_ = new Vertex[GetNumEdges() * 2];
  core_num_ = new Size[GetNumVertices()];
  std::fill(label_frequency_, label_frequency_ + data.GetNumLabels(), 0);

  Size cur_idx = 0;

  // transfer label & construct adj list and label frequency
  for (Vertex v = 0; v < GetNumVertices(); ++v) {
    Label l = data.GetTransferredLabel(label_[v]);
    if (l == INVALID_LB) return false;
    if (label_frequency_[l] == 0) num_label_ += 1;
    label_[v] = l;
    max_label_ = std::max(max_label_, l);
    label_frequency_[l] += 1;
    if (adj_list[v].size() > max_degree_) max_degree_ = adj_list[v].size();

    start_off_[v] = cur_idx;
    core_num_[v] = adj_list[v].size();

    std::copy(adj_list[v].begin(), adj_list[v].end(),
              linear_adj_list_ + cur_idx);

    cur_idx += adj_list[v].size();
  }
  start_off_[GetNumVertices()] = num_edge_ * 2;

  // preprocess for query graph
  computeCoreNum();

  is_tree_ = true;
  for (Vertex v = 0; v < GetNumVertices(); ++v) {
    if (GetCoreNum(v) > 1) {
      is_tree_ = false;
      break;
    }
  }

  ExtractResidualStructure();

  return true;
}

namespace {
struct NECInfo {
  bool visit = false;
  Vertex representative;
  Size NEC_elems_idx;
};
}  // namespace

void QueryGraph::ExtractResidualStructure() {
  NECInfo *NEC_infos_temp = new NECInfo[GetNumVertices() * (max_label_ + 1)];
  NEC_elems_ = new NECElement[GetNumVertices()];
  NEC_map_ = new Vertex[GetNumVertices()];
  NEC_size_ = new Size[GetNumVertices()];

  Size num_NEC_elems_ = 0;

  std::fill(NEC_map_, NEC_map_ + GetNumVertices(), INVALID_VTX);
  std::fill(NEC_size_, NEC_size_ + GetNumVertices(), 0);

  num_non_leaf_vertices_ = GetNumVertices();

  // construct NEC map
  for (Vertex v = 0; v < GetNumVertices(); ++v) {
    if (GetDegree(v) == 1) {
      Vertex p = GetNeighbor(GetStartOffset(v));
      Label l = GetLabel(v);

      NECInfo &info = NEC_infos_temp[GetNumVertices() * l + p];
      if (!info.visit) {
        info = {true, v, num_NEC_elems_};
        NEC_map_[v] = v;

        for (Size nbr_idx = GetStartOffset(p); nbr_idx < GetEndOffset(p);
             ++nbr_idx) {
          Vertex nbr = GetNeighbor(nbr_idx);
          if (nbr == v) {
            NEC_elems_[num_NEC_elems_] = {l, p, v, 0,
                                          nbr_idx - GetStartOffset(p)};
            break;
          }
        }
        num_NEC_elems_ += 1;
      } else {
        NEC_map_[v] = info.representative;
      }
      NEC_size_[info.representative] += 1;
      NEC_elems_[info.NEC_elems_idx].size += 1;
      num_non_leaf_vertices_ -= 1;
    } else {
      NEC_size_[v] += 1;
    }
  }

  for (Vertex v = 0; v < GetNumVertices(); ++v) {
    NEC_size_[v] = NEC_size_[GetNECRepresentative(v)];
  }

  num_NEC_label_ = 0;
  NEC_start_offs_ = nullptr;
  if (num_NEC_elems_ > 0) {
    // sort NEC elems by label
    std::sort(NEC_elems_, NEC_elems_ + num_NEC_elems_,
              [](const NECElement &a, const NECElement &b) -> bool {
                return a.label < b.label;
              });

    // construct start offsets of NEC elems for same label
    NEC_start_offs_ = new Size[GetNumVertices() + 1];
    NEC_start_offs_[0] = 0;
    num_NEC_label_ += 1;

    Label prev_label = NEC_elems_[0].label;
    for (Size i = 1; i < num_NEC_elems_; ++i) {
      if (NEC_elems_[i].label != prev_label) {
        prev_label = NEC_elems_[i].label;
        NEC_start_offs_[num_NEC_label_] = i;
        num_NEC_label_ += 1;
      }
    }
    NEC_start_offs_[num_NEC_label_] = num_NEC_elems_;
  }

  delete[] NEC_infos_temp;
}
}  // namespace daf
