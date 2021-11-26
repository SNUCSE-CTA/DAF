#include <algorithm>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "global/timer.h"
#include "include/backtrack.h"
#include "include/candidate_space.h"
#include "include/dag.h"
#include "include/data_graph.h"
#include "include/query_graph.h"

int main(int argc, char* argv[]) {
  daf::Timer total_timer, backtrack_timer;

  std::string dir_name;
  std::string data_name;
  std::string query_name;
  uint64_t limit = std::numeric_limits<uint64_t>::max();

  for (int i = 1; i < argc; ++i) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
        case 't':
          dir_name = argv[i + 1];
          break;
        case 'd':
          data_name = argv[i + 1];
          break;
        case 'q':
          query_name = argv[i + 1];
          break;
        case 'm':
          limit = std::atoi(argv[i + 1]);
      }
    }
  }

  daf::DataGraph data;
  data.LoadAndProcessGraph(data_name);
  daf::QueryGraph query;

  total_timer.Start();
  query.LoadAndProcessGraph(query_name, data);
  total_timer.Stop();

  daf::DAG dag(data, query);

  total_timer.Start();
  dag.BuildDAG();
  total_timer.Stop();

  daf::CandidateSpace cs(data, query, dag);

  total_timer.Start();
  bool cs_constructed = cs.BuildCS();
  total_timer.Stop();

  uint64_t num_embeddings = 0;
  if (cs_constructed) {
    daf::Backtrack backtrack(data, query, cs);

    backtrack_timer.Start();
    num_embeddings = backtrack.FindMatches(limit);
    backtrack_timer.Stop();
  }

  total_timer.Add(backtrack_timer);

  std::cout << "query vertices : " << query.GetNumVertices() << "\n";
  std::cout << "query non leaf vertices : " << query.GetNumNonLeafVertices()
            << "\n";
  std::cout << "#Embeddings : " << num_embeddings << "\n";

  std::cout << "Total time: " << total_timer.GetTime() << "\n";
  std::cout << "Backtrack time: " << backtrack_timer.GetTime() << "\n";
}
