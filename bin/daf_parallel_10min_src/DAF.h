#ifndef GLOBALVARIABLES_H_
#define GLOBALVARIABLES_H_
#define PARALLEL

//#define WINDOWS

//#define OUTPUT_EXTRA_INFO
//#define DEBUG_MODE
//#define MAPPING_FUNCTION_LOG
//#define REDUCE_COUNT_READING
//#define PARALLEL_TEST_1

#include <algorithm>
#include <cfloat>
#include <chrono>
#include <climits>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>
//#include "bitvector.h"
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#ifdef WINDOWS
//#include "int_vector.hpp"
#include <thread>
#else
#include <pthread.h>
#include <unistd.h>

#include <sdsl/bit_vectors.hpp>  //bit_vector
#endif
#include <limits.h>  //for LLONG_MAX which is 9223372036854775807
#ifdef PARALLEL
#include <omp.h>
#endif

using namespace std;
using namespace sdsl;
typedef long long weight_type;

// forward declaration
struct NEC_element;
struct NEC_set_array_element;
// The index unit for the global index built when explore candidate region
struct NodeIndexUnit {
  int size;         // the size for both path and candidates
  int* candidates;  // candidate set

  // for DAG. "backtrack_" is prefixed to the pre-exiting variables.
  int*** backtrack_index = NULL;  // backtrack_index[i][j] = candidates of this
                                  // unit when the i-th parent, regarding DAG,
                                  // mapped to the j-th candidate of the parent.
  int** backtrack_size_of_index =
      NULL;  // backtrack_size_of_index[i][j] = size of back_trak_index[i][j].
             // That is, the number of candidates of this unit when the i-th
             // parent mapped to the j-th candidate of the parent.
  int* backtrack_parent_cand_size =
      NULL;  // backtrack_parent_cand_size[i] = the number of candidates of the
             // i-th parent.
  long long* weight = NULL;
};

// consts
const int MAX_QUERY_NODE = 200;
int SIZEOF_INT = sizeof(int) * 8;

// global flags
bool use_path_size_order = true;
bool use_candidate_size_order = false;
bool order_flag =
    false;  // set use_candidate_size_order in readAndProcessDataGraph when
            // order_flag turns false. (order_flag turns true when program
            // argument specifies matching order)
bool use_failing_set = true;
bool isTree;

// variables for data graph
int cnt_node = 0;           // The number of nodes in the data graph
int count_edge_number = 0;  // the number of edges in the data graph
int* nodes_label = NULL;
int cnt_unique_label;  // the number of unique labels
int* label_freqency = NULL;
int* transferredLabel = NULL;
int* nodes_info = NULL;  // this array stores the start position of each node's
                         // adjacent list in the "nodes_data"
pair<int, int>* nodes_to_label_info =
    NULL;  // given a node, and a label, it can extract the children of the
           // given node with the given label
pair<int, int>* label_deg_label_pos;
vector<int> degree_array;
int* core_number_data = NULL;
int* degree_data = NULL;
int MAX_DEGREE;
int* nodes_data = NULL;
int NLF_size = -1;
int* NLF_check = NULL;  // for data
int* NLF_array = NULL;  // for query
int* MAX_NB_degree_data;

// variables for query graph
ifstream fin_query;
int MAX_DEGREE_QUERY = 0;
int label_cur;
int degree_cur;
vector<int> query_nodes_array;
vector<int> core_tree_node_child_array;
vector<int> core_tree_node_nte_array;
int nodes_label_query[MAX_QUERY_NODE];
int node_degree_query[MAX_QUERY_NODE];  // store the degree of each node in the
                                        // query graph
int core_number_query[MAX_QUERY_NODE];
int query_nodes_array_info[MAX_QUERY_NODE + 1];
// coreDecomposition
int bin_query[MAX_QUERY_NODE];
int pos_query[MAX_QUERY_NODE];
int vert_query[MAX_QUERY_NODE];
// extractResidualStructures
int residual_tree_match_seq_index = 0;
int residual_tree_leaf_node_index = 0;
int NEC_mapping_pair_index = 0;
int NEC_map[MAX_QUERY_NODE];
int* NEC_mapping = NULL;
int* NEC_mapping_Actual;
NEC_element* NEC_mapping_pair;
NEC_set_array_element* NEC_set_array;
int NEC_set_index = 0;
char visited_for_query[MAX_QUERY_NODE];  // need memset set before use
int dfs_stack_query[MAX_QUERY_NODE];
int residual_tree_match_seq[MAX_QUERY_NODE];
pair<int, double> residual_tree_leaf_node[MAX_QUERY_NODE];
int* tree_node_parent;
// start_node_selection
int root_node_id;  // the selected root node of current query graph
// BFS_DAG
int core_tree_node_child_array_index = 0;
int core_tree_node_nte_array_index = 0;
int bfs_sequence_index = 0;  // 20170503
int exploreCRSequence_indx = 0;
int visited_int_for_query[MAX_QUERY_NODE];  // need memset set before use
int queue_array_query[MAX_QUERY_NODE];
int BFS_level_query[MAX_QUERY_NODE];
int BFS_parent_query[MAX_QUERY_NODE];
int* DAG_parent_query[MAX_QUERY_NODE];  // DAG_parent_query[i] contains parent
                                        // nodes, regarding DAG of query graph,
                                        // of i-node.
int* DAG_child_query[MAX_QUERY_NODE];   // DAG_child_query[i] contains child
                                        // nodes, regarding DAG of query graph,
                                        // of i-node.
int* DAG_child_query_parent_index
    [MAX_QUERY_NODE];  // DAG_child_query_parent_index[i][j] contains parent,
                       // which is i-th node, index of the j-th child node of
                       // i-node. It is used for ajacency list prunning.
int DAG_parent_query_size
    [MAX_QUERY_NODE];  // DAG_parent_query_size[i] contains the number of parent
                       // nodes, regarding DAG of query graph, of i-node.
int DAG_child_query_size
    [MAX_QUERY_NODE];  // DAG_child_query_size[i] contains the number of child
                       // nodes, regarding DAG of query graph, of i-node.
bit_vector DAG_ancestor[MAX_QUERY_NODE];
int* label_frequency_rank;  // label_frequency_rank[i] contains label i's rank
                            // in terms of the label_frequency in increasing
                            // order; i.e., label_frequency_rank[0] has the
                            // lowest label_frequency.
int simulation_sequence[MAX_QUERY_NODE];
int simulation_sequence_index = 0;
vector<pair<int, int>> level_index;
// for failing set
// sdsl::bit_vector zero_vector;
// sdsl::bit_vector one_vector;
bit_vector zero_vector;
bit_vector one_vector;
int prev_cnt_node_query = -1;

// Dynamic Programming
NodeIndexUnit indexSet[MAX_QUERY_NODE];

// other global variables
#ifndef PARALLEL
int* global_temp_array_1;
int* global_temp_array_2;
#endif
int count_global_temp_array_1;
int count_global_temp_array_2;

struct NEC_Node {
  int node;
  NEC_Node* nextAddress;
};

struct SearchUnit {
  int* address;
  int address_size;
  int address_pos;
};

struct PointerDagSearchUnit {
  int* address = NULL;
  int address_size;
  int address_pos;
  int vertex;

  // for failing set
  bool firstly_visited = true;
  // sdsl::bit_vector return_vector;
  bit_vector return_vector;
};

// consts
const int INF = numeric_limits<int>::max();
const int MAX_QUERY_DEGREE = 100;

int max_label_counter = 0;

#ifndef PARALLEL
double LIMIT_REMAINED = 0;
int curretn_LLC_label = -1;
//=== for Maximum Matching =====
int* nec_mappiing_MM;
int nodes_MM_idx = 0;
#endif

const int NIL = 0;
#ifndef PARALLEL
int* pair_U;
int* pair_V;
int* u_cand_MM;
int dist[MAX_QUERY_NODE + 1];
#endif

// for the one time MM optimization
int Leaf_cands_size = 0;
#ifndef PARALLEL
int Leaf_cands_idx = 0;
int idx_sum_nec_cands;
vector<int> Leaf_cands;
pair<int, int>* Leaf_cands_info;
int** sum_nec_cands_array;  // added for the optimization of performing all MM
                            // first
int* size_sum_nec_cands;    // added for the optimization of performing all MM
                            // first
int* idx_v_cands;  // length: cnt_node. its each entry stores the number of
                   // query nodes contain the corresponding candidate.
vector<pair<int, int>>*
    v_cands_info;  // a combination of v_cands and v_cands_pos
pair<int, int>*
    u_cand_info;  // two-d array:
                  //[i,j] means the j-th candidate pair of the i-th nec node(in
                  // the nec-node categorized by label) each candidate pair
                  // stores a candidate and the position of "the i-th nec node"
                  // in this candidate's u_cand set.
int* flag_sum_nec_cands;  // indicate whether or not a candidate has been added
                          // into the candidate already.
#endif
inline void swap_value(int& v1, int& v2) {
  int temp = v1;
  v1 = v2;
  v2 = temp;
}

inline void swap_value(pair<int, int>& v1, pair<int, int>& v2) {
  pair<int, int> temp = v1;
  v1 = v2;
  v2 = temp;
}

//==============
#ifndef PARALLEL
int nec_count_set_size = -1;
#endif

long long recursive_call_count = 0;
long long recursive_call_count_sum = 0;
#ifndef PARALLEL
pair<int, pair<int, double>> NEC_set_ranking[MAX_QUERY_NODE];
int nec_region_size[MAX_QUERY_NODE];
char local_flag_query[MAX_QUERY_NODE];
pair<int, int>
    v_nec_count[MAX_QUERY_NODE];  // used in the map special tree function
PointerDagSearchUnit ptr_dag_su[MAX_QUERY_NODE];
int self_pos[MAX_QUERY_NODE];
int actual_mapping[MAX_QUERY_NODE];
int nec_count_set[MAX_QUERY_NODE];
int mapped_parent[MAX_QUERY_NODE];
int* iec[MAX_QUERY_NODE][MAX_QUERY_DEGREE - 1];
int iecSize[MAX_QUERY_NODE][MAX_QUERY_DEGREE - 1];
// for failing set
bool exist_u[MAX_QUERY_NODE][MAX_QUERY_NODE];
int ancestor_set[MAX_QUERY_NODE][MAX_QUERY_NODE];
int ancestor_set_index[MAX_QUERY_NODE];
#endif

weight_type WEIGHT_MAX = LLONG_MAX;
#ifndef PARALLEL
int frontier_node_idx = 0;
int frontier_min_index = -1;
long long frontier_min_weight = LLONG_MAX;
double PRE_COMPUTED_PERMUTATION;
#endif

//========== some flag arrays =====================
char* flag_prelin_char;  // data count
int* flag_child_cand;    // data count

//========== query graph variables =========
int cnt_node_query = 0;  // the number of nodes in the query graph
#ifdef PARALLEL
inline bool sort_by_pair_of_pair(pair<pair<int, int>, int>,
                                 pair<pair<int, int>, int>);
class Context {
 public:
  int mapped_parent[MAX_QUERY_NODE];
  bool mapped_query[MAX_QUERY_NODE];
  int nec_region_size[MAX_QUERY_NODE];
  char local_flag_query[MAX_QUERY_NODE];
  pair<int, int> v_nec_count[MAX_QUERY_NODE];
  PointerDagSearchUnit ptr_dag_su[MAX_QUERY_NODE];
  int self_pos[MAX_QUERY_NODE];
  int actual_mapping[MAX_QUERY_NODE];
  int nec_count_set[MAX_QUERY_NODE];
  pair<int, pair<int, double>> NEC_set_ranking[MAX_QUERY_NODE];
  // failing set
  bool exist_u[MAX_QUERY_NODE][MAX_QUERY_NODE];
  int ancestor_set[MAX_QUERY_NODE][MAX_QUERY_NODE];
  int ancestor_set_index[MAX_QUERY_NODE];
  //
  int frontier_node_idx;
  int frontier_min_index;
  long long frontier_min_weight;
  int frontier_node[MAX_QUERY_NODE];
  weight_type weight_array[MAX_QUERY_NODE];
  weight_type cur_weight[MAX_QUERY_NODE];
  int position[MAX_QUERY_NODE];
  int num_recent_insert[MAX_QUERY_NODE];
  int dist[MAX_QUERY_NODE + 1];

  pair<pair<int, int>, int> pp_array[MAX_QUERY_NODE];

  int*** iec;
  int** iecSize;
  char* mapping_flag_data;
  int* mapping_flag_data_int;
  int* global_temp_array_1;
  int* global_temp_array_2;

  vector<int> Leaf_cands;
  pair<int, int>* Leaf_cands_info;
  int* nec_mappiing_MM;
  int* u_cand_MM;
  int* pair_U;
  int* pair_V;
  pair<int, int>* u_cand_info;  // three-d array:
  int* idx_v_cands;  // length: cnt_node. its each entry stores the number of
                     // query nodes contain the corresponding candidate.
  int* flag_sum_nec_cands;  // indicate whether or not a candidate has been
                            // added into the candidate already.
  vector<pair<int, int>>*
      v_cands_info;           // a combination of v_cands and v_cands_pos
  int** sum_nec_cands_array;  // added for the optimization of performing all MM
                              // first
  int* size_sum_nec_cands;    // added for the optimization of performing all MM
                              // first

  int nec_count_set_size;
  int Leaf_cands_idx;
  int idx_sum_nec_cands;
  double PRE_COMPUTED_PERMUTATION;
  double LIMIT_REMAINED;
  int nodes_MM_idx;
  int curretn_LLC_label;

  Context() {
    frontier_node_idx = 0;
    frontier_min_index = -1;
    frontier_min_weight = LLONG_MAX;

    iec = new int**[MAX_QUERY_NODE];
    for (int j = 0; j < MAX_QUERY_NODE; j++) {
      iec[j] = new int*[MAX_QUERY_DEGREE - 1];
      for (int k = 0; k < MAX_QUERY_DEGREE - 1; k++)
        iec[j][k] = new int[max_label_counter];
    }
    iecSize = new int*[MAX_QUERY_NODE];
    for (int j = 0; j < MAX_QUERY_NODE; j++)
      iecSize[j] = new int[MAX_QUERY_DEGREE - 1];

    mapping_flag_data = new char[cnt_node];
    memset(mapping_flag_data, 0, sizeof(char) * cnt_node);
    if (use_failing_set) {
      mapping_flag_data_int = new int[cnt_node];
      memset(mapping_flag_data_int, -1, sizeof(int) * cnt_node);
      memset(ancestor_set_index, 0, sizeof(int) * MAX_QUERY_NODE);
      for (int i = 0; i < MAX_QUERY_NODE; ++i) {
        memset(exist_u[i], false, sizeof(bool) * MAX_QUERY_NODE);
      }
    }
    global_temp_array_1 = new int[cnt_node];

    global_temp_array_2 = new int[cnt_node];

    Leaf_cands.resize(MAX_QUERY_NODE * max_label_counter);  // tight upper bound
    Leaf_cands_info =
        new pair<int, int>[(cnt_unique_label + 1) * MAX_QUERY_NODE];
    nec_mappiing_MM = new int[MAX_QUERY_NODE + 1];
    u_cand_MM = new int[MAX_QUERY_NODE * max_label_counter];

    pair_U = new int[MAX_QUERY_NODE + 1];
    memset(pair_U, NIL, sizeof(int) * NIL);  // must reset before using

    pair_V = new int[cnt_node];
    memset(pair_V, NIL, sizeof(int) * cnt_node);  // must reset before using

    u_cand_info =
        new pair<int, int>[MAX_QUERY_NODE *
                           max_label_counter];  // definite upper bound]

    memset(local_flag_query, 0, sizeof(char) * MAX_QUERY_NODE);

    idx_v_cands = new int[cnt_node];
    memset(idx_v_cands, 0, sizeof(int) * cnt_node);

    flag_sum_nec_cands = new int[cnt_node];
    memset(flag_sum_nec_cands, -1, sizeof(int) * cnt_node);

    v_cands_info = new vector<pair<int, int>>[cnt_node];
    for (int i = 0; i < cnt_node; i++) v_cands_info[i].resize(MAX_QUERY_NODE);

    sum_nec_cands_array = new int*[(cnt_unique_label + 1)];
    for (int i = 0; i < (cnt_unique_label + 1); i++)
      sum_nec_cands_array[i] = new int[max_label_counter];

    size_sum_nec_cands = new int[(cnt_unique_label + 1)];

    nec_count_set_size = -1;
    Leaf_cands_idx = 0;
    LIMIT_REMAINED = 0;
    nodes_MM_idx = 0;
    curretn_LLC_label = -1;
  }

  inline void sort_v_cands_info(int max_cand, int len) {
    for (int x = 0; x < len; ++x)
      pp_array[x] = make_pair(v_cands_info[max_cand][x],
                              nec_region_size[v_cands_info[max_cand][x].first]);
    sort(pp_array, pp_array + len, sort_by_pair_of_pair);
    for (int x = 0; x < len; ++x) v_cands_info[max_cand][x] = pp_array[x].first;
  }

  inline void print_frontier();
  inline void reinsert_to_frontier(int, weight_type, int);
  inline void insert_to_frontier_first(int, weight_type);
  inline void pop_from_frontier(int&, weight_type&, int&);
  inline void remove_recently_inserted_from_frontier(int);
  inline void clear_frontier();
};
Context** ctxt;
#else
char* mapping_flag_data;
int* mapping_flag_data_int;
bool mapped_query[MAX_QUERY_NODE];
int frontier_node[MAX_QUERY_NODE];
weight_type weight_array[MAX_QUERY_NODE];
weight_type cur_weight[MAX_QUERY_NODE];
int position[MAX_QUERY_NODE];
int num_recent_insert[MAX_QUERY_NODE];
#endif

//============ variables for current query graph =======================
double mapping_found;

//======= simulation arrays =====
int* array_to_clean;
int to_clean_index = 0;

//============================

//========== DATA graph variables =========

//================ need to clean ====================================
int count_query_file;
string datagraphFile;
string querygraphFileFolder;
int largest_label = 0;

long LIMIT = 100000;  // deafult 100,000
int nThreads = 1;     // 20181211

//========== variables for the data graph
//========================================
unordered_set<int> unique_label_set;  // store the unique labels
vector<pair<int, int>>
    NEC_set_by_label_index;  // include a redunant element to set the end

//==================================================================
int sum_degree_cur;
int MAX_sum_degree_cur = 0;

// The following three are the temporary arrays for the indexSet element, used
// during the global candidate exploration
int* index_array_for_indexSet;
int count_index_array_for_indexSet;

// statistics for evaluation
class CUtility {
 private:
  std::chrono::time_point<std::chrono::high_resolution_clock> begin;

 public:
  CUtility() {}
  ~CUtility() {}
  // begin to counting
  inline void startCT() { begin = std::chrono::high_resolution_clock::now(); }
  // stop counting
  inline double endCT() {
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - begin;
    return elapsed.count();
  }
};

CUtility* cu_total = new CUtility();
CUtility* cu_querying = new CUtility();

#define TOTAL_BEGIN cu_total->startCT();
#define TOTAL_END time_total = cu_total->endCT();

#define SEARCH_BEGIN cu_querying->startCT();
#define SEARCH_END time_search = cu_querying->endCT();

// some result collecting parameters used in the main function
double sum_search_time = 0;
double sum_mapping = 0;
double sum_time_total = 0;

double time_search;
double time_total = 0;

double longest_query_time = 0.00;
int longest_query_id;

struct NEC_set_array_element {
  int parent_id;
  int represent_node;
  int sum;
  NEC_set_array_element() {}
  NEC_set_array_element(int parent_id, int represent_node, int sum) {
    this->parent_id = parent_id;
    this->represent_node = represent_node;
    this->sum = sum;
  }
};

struct NEC_element {
  int label;
  int parent_id;
  int represent_node;

  NEC_element() {}

  NEC_element(int label, int parent_id, int represent_node) {
    this->label = label;
    this->parent_id = parent_id;
    this->represent_node = represent_node;
  }
};

inline unsigned int split(const string& txt, vector<string>& strs, char ch) {
  // this is the general case
  size_t pos = txt.find(ch);
  size_t initialPos = 0;
  strs.clear();
  // Decompose statement
  while (pos != string::npos) {
    strs.push_back(txt.substr(initialPos, pos - initialPos + 1));
    initialPos = pos + 1;
    pos = txt.find(ch, initialPos);
  }
  // Add the last one
  strs.push_back(txt.substr(initialPos, min(pos, txt.size()) - initialPos + 1));
  // return the size of the vector
  return strs.size();
}

//======================= sorting functions
//========================================
inline bool sortByLabelAndDegree(int node1, int node2) {
  if (nodes_label[node1] == nodes_label[node2])
    return (degree_data[node1] > degree_data[node2]);
  else
    return (nodes_label[node1] < nodes_label[node2]);
}

inline bool sortByLabelandDegree(int node1, int node2) {
  if (nodes_label[node1] == nodes_label[node2])
    return (degree_data[node1] < degree_data[node2]);
  else
    return (nodes_label[node1] < nodes_label[node2]);
}

inline bool sortLLCnode(pair<int, pair<int, double>> a,
                        pair<int, pair<int, double>>
                            b) {  // the second pair is <node sum, cand sum>
  if (a.second.first == b.second.first)
    return a.second.second < b.second.second;
  else
    return a.second.first < b.second.first;
}
#ifdef PARALLEL
inline bool sort_by_pair_of_pair(pair<pair<int, int>, int> p1,
                                 pair<pair<int, int>, int> p2) {
  return p1.second < p2.second;
}
#else
inline bool sort_by_effective_cand_size_pair(pair<int, int> p1,
                                             pair<int, int> p2) {
  return nec_region_size[p1.first] < nec_region_size[p2.first];
}
#endif
inline bool sort_by_NEC_label(NEC_element p1, NEC_element p2) {
  return (p1.label < p2.label);
}

inline bool sort_by_second_element(pair<int, int> p1, pair<int, int> p2) {
  return (p1.second < p2.second);
}

inline bool sortByDegree_Query_dec(const int id1, const int id2) {
  return (node_degree_query[id1] > node_degree_query[id2]);
}
// gmgu
inline bool sortByLabelFrequency_Label_inc(const int id1, const int id2) {
  return (label_freqency[id1] < label_freqency[id2]);
}
inline bool sortByLabelFrequencyRank_Query_inc(const int id1, const int id2) {
  return (label_frequency_rank[nodes_label_query[id1]] <
          label_frequency_rank[nodes_label_query[id2]]);
}

inline string withCommasNoDot(double value) {
  string numWithCommas = to_string(value);
  int pos = numWithCommas.find_first_of(".");
  int insertPosition;
  if (pos == string::npos)
    insertPosition = numWithCommas.length() - 3;
  else
    insertPosition = pos - 3;
  while (insertPosition > 0) {
    numWithCommas.insert(insertPosition, ",");
    insertPosition -= 3;
  }
  pos = numWithCommas.find_first_of(".");
  return numWithCommas.substr(0, pos);
}

inline double factorization(int x) {
  //=============testing===========
  //	if (x > testing_high_factor)
  //		testing_high_factor = x;
  //===============================

  switch (x) {
    case 1:
      return 1;
    case 2:
      return 2;
    case 3:
      return 6;
    case 4:
      return 24;
    case 5:
      return 120;
    case 6:
      return 720;
    case 7:
      return 5040;
    case 8:
      return 40320;
    case 9:
      return 362880;
    case 10:
      return 3628800;
    case 11:
      return 39916800;
    case 12:
      return 479001600;
    case 13:
      return 6227020800;
    case 14:
      return 87178291200;
    case 15:
      return 1307674368000;
    case 16:
      return 20922789888000;
    case 17:
      return 355687428096000;
    case 18:
      return 6402373705728000;
    case 19:
      return 121645100408832000;
    case 20:
      return 2432902008176640000;
  }

  double result = 2432902008176640000;
  for (double i = 21; i <= x; i++) result *= i;

  return result;
}

inline void initializeStatisticParameters() {
  time_search = 0;
  time_total = 0;
}

inline void addUpStatisticParameters(int i) {
  sum_time_total += time_total;
  sum_search_time += time_search;
  sum_mapping += mapping_found;

  if (time_search > longest_query_time) {
    longest_query_time = time_search;
    longest_query_id = i;
  }
}

void getLimit_full(string str_full_limit, long& LIMIT) {
  if (str_full_limit == "1K")
    LIMIT = 1000;
  else if (str_full_limit == "10K")
    LIMIT = 10000;
  else if (str_full_limit == "100K")
    LIMIT = 100000;
  else if (str_full_limit == "100M")
    LIMIT = 100000000;
  else if (str_full_limit == "1B")
    LIMIT = 100000000000;
  else
    LIMIT = atol(str_full_limit.c_str());
}

inline void parametersInitilisingBeforeQuery() {
  for (int i = 0; i < cnt_unique_label + 1; i++)
    if (label_freqency[i] > max_label_counter)
      max_label_counter = label_freqency[i];

  max_label_counter++;  // increase by one for safety in case of overflow

  //================ MAX_CAND and MM optimization =====================
  Leaf_cands_size = MAX_QUERY_NODE * max_label_counter;
#ifdef PARALLEL
  ctxt = new Context*[nThreads];
  for (int x = 0; x < nThreads; ++x) ctxt[x] = new Context();
#else
  Leaf_cands.resize(MAX_QUERY_NODE * max_label_counter);  // tight upper bound

  Leaf_cands_info = new pair<int, int>[(cnt_unique_label + 1) * MAX_QUERY_NODE];
  nec_mappiing_MM = new int[MAX_QUERY_NODE + 1];
  u_cand_MM =
      new int[MAX_QUERY_NODE * max_label_counter];  // definite upper bound]

  pair_U = new int[MAX_QUERY_NODE + 1];
  memset(pair_U, NIL, sizeof(int) * NIL);  // must reset before using

  pair_V = new int[cnt_node];
  memset(pair_V, NIL,
         sizeof(int) * max_label_counter);  // must reset before using

  u_cand_info = new pair<int, int>[MAX_QUERY_NODE *
                                   max_label_counter];  // definite upper bound]

  memset(local_flag_query, 0, sizeof(char) * MAX_QUERY_NODE);
  idx_v_cands = new int[cnt_node];
  memset(idx_v_cands, 0, sizeof(int) * cnt_node);

  flag_sum_nec_cands = new int[cnt_node];
  memset(flag_sum_nec_cands, -1, sizeof(int) * cnt_node);

  v_cands_info = new vector<pair<int, int>>[cnt_node];
  for (int i = 0; i < cnt_node; i++) v_cands_info[i].resize(MAX_QUERY_NODE);

  sum_nec_cands_array = new int*[(cnt_unique_label + 1)];
  for (int i = 0; i < (cnt_unique_label + 1); i++)
    sum_nec_cands_array[i] = new int[max_label_counter];
  size_sum_nec_cands = new int[(cnt_unique_label + 1)];
#endif

  //============================================================
  memset(NEC_mapping, 0, sizeof(int) * (cnt_unique_label + 1) * MAX_QUERY_NODE);

  memset(visited_for_query, 0, sizeof(char) * MAX_QUERY_NODE);

  NEC_mapping_Actual = new int[(cnt_unique_label + 1) *
                               MAX_QUERY_NODE];  // label and the parent node id
  memset(NEC_mapping_Actual, 0,
         sizeof(int) * (cnt_unique_label + 1) * MAX_QUERY_NODE);

  for (int i = 0; i < MAX_QUERY_NODE; i++) {
    indexSet[i].candidates = new int[max_label_counter];
    indexSet[i].weight = new long long[max_label_counter];  // 20170414
#ifndef PARALLEL
    for (int j = 0; j < MAX_QUERY_DEGREE - 1; j++) {
      iec[i][j] = new int[max_label_counter];
    }
#endif
  }

  index_array_for_indexSet = new int[cnt_node];

  tree_node_parent = new int[MAX_QUERY_NODE];

#ifndef PARALLEL
  mapping_flag_data = new char[cnt_node];
  memset(mapping_flag_data, 0, sizeof(char) * cnt_node);
  if (use_failing_set) {
    mapping_flag_data_int = new int[cnt_node];
    memset(mapping_flag_data_int, -1, sizeof(int) * cnt_node);
    memset(ancestor_set_index, 0, sizeof(int) * MAX_QUERY_NODE);
    for (int i = 0; i < MAX_QUERY_NODE; ++i) {
      memset(exist_u[i], false, sizeof(bool) * MAX_QUERY_NODE);
    }
  }

  global_temp_array_1 = new int[cnt_node];
  global_temp_array_2 = new int[cnt_node];
#endif

  memset(NEC_map, 0, sizeof(int) * MAX_QUERY_NODE);

  //===simulation parmateters
  array_to_clean = new int[cnt_node * MAX_QUERY_NODE];
  //==========================================

  flag_prelin_char = new char[cnt_node];
  memset(flag_prelin_char, 0, sizeof(char) * cnt_node);

  flag_child_cand = new int[cnt_node];
  memset(flag_child_cand, -1, sizeof(int) * cnt_node);
}

#endif /* GLOBALVARIABLES_H_ */
