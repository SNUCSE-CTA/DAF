#include <algorithm>
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

  std::string data_name;
  std::string query_name;
  uint64_t limit = std::numeric_limits<uint64_t>::max();

  for (int i = 1; i < argc; ++i) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
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

  std::cout << "Loading data graph...\n";
  daf::DataGraph data(data_name);
  data.LoadAndProcessGraph();

  std::cout << "Loading query graph...\n";
  daf::QueryGraph query(query_name);

  total_timer.Start();
  query.LoadAndProcessGraph(data);
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
  uint64_t num_backtrack_calls = 0;
  if (cs_constructed) {
    std::cout << "Enumerating...\n";
    daf::Backtrack backtrack(data, query, cs);

    backtrack_timer.Start();
    num_embeddings = backtrack.FindMatches(limit);
    backtrack_timer.Stop();
    num_backtrack_calls = backtrack.GetNumBacktrackCalls();
  }

  total_timer.Add(backtrack_timer);

  std::cout << "#Matches: " << num_embeddings << "\n";
  std::cout << "#Recursive calls: " << num_backtrack_calls << "\n";

  std::cout << "Total time: " << total_timer.GetTime() << " ms\n";
  std::cout << "Search time: " << backtrack_timer.GetTime() << " ms\n";
}
