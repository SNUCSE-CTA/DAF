# DAF [SIGMOD 2019]
Efficient Subgraph Matching: Harmonizing Dynamic Programming, Adaptive Matching Order, and Failing Set Together

## Build
```
mkdir build
cd build
cmake ..
make
```
## Run
```
Options:
-d,   specify data graph file name
-q,   specify query graph file name
-m,   specify the number of matches to find (if not specified, find all matches)
```

```
Example:
./build/main/DAF -d dataset_example/data.graph -q dataset_example/query.graph -m 10
```

## Input File Format

The graph file format for data graphs and query graphs is a text format to store an vertex-labeled undirected graph. 
- The first line of the file should be "t #vertices #edges".
- Following lines of "v vertex-ID vertex-label" indicate the vertices in the graph.
- The vertices should be written in the file in ascending order of their IDs, and a vertex ID should be in [0, #vertices - 1].
- Following lines of "e vertex-ID1 vertex-ID2" after the vertices indicate the undirected edges in the graph.

You can find example data graph and query graph files under the dataset_example directory.

For example:
```
t 4 5
v 0 0
v 1 1
v 2 2
v 3 3
e 0 1
e 0 2
e 1 2
e 1 3
e 2 3
```
