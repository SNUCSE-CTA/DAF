# DAF [SIGMOD 2019]
Efficient Subgraph Matching: Harmonizing Dynamic Programming, Adaptive Matching Order, and Failing Set Together

## Binary Files
Binary files of DAF are available for linux
- daf_10min: DAF for subgraph matching, which sets time limit of 10 minutes for each query graph in a query set.
- daf_twitter_10min: daf_10min for a large data graph
- daf_parallel_10min: parallel version of DAF using openMP.
- daf_parallel: parallel version of DAF without time limit.
- daf_twitter: daf for a large data graph where the time limit can be given as an argument. 
- daf_parallel_twitter: parallel version of DAF for a large data graph where the time limit can be given as an argument.

The maximum number of query vertices is set to 400 for the binary file.
## Run
```
Options:
-d,   specify data graph file name
-q,   specify query graph file name
-n,   specify the number of queries in a query file
-m,   specify the number of matches to find
-p,   use path-size order. Path-size order is used by default.
-c,   use candidate-size order.
-f,   'not' using failing set. Failing set is used by default unless specified.
-h,   #threads (for daf_parallel_10min)
-l,   specify the time limit in seconds (default: 600 seconds) (for daf_twitter and daf_parallel_twitter)
```

```
Usages:
./daf_10min -d datagraphFile -q querygraphFile -n #queries -m #matches // path-size order and failing set are used by default.
./daf_10min -d datagraphFile -q querygraphFile -n #queries -m #matches -f // path-size order is used, and failing set is not used.
./daf_10min -d datagraphFile -q querygraphFile -n #queries -m #matches -c // candidate-size order and failing set are used
./daf_10min -d datagraphFile -q querygraphFile -n #queries -m #matches -c -f // candidate-size order is used, and failing set is not used.

./daf_parallel_10min -d datagraphFile -q querygraphFile -n #queries -m #matches -h #threads // using #threads
./daf_parallel -d datagraphFile -q querygraphFile -n #queries -m #matches -h #threads // without time limit

./daf_twitter -d datagraphFile -q querygraphFile -n #queries -m #matches -l timeLimit 
./daf_parallel_twitter -d datagraphFile -q querygraphFile -n #queries -m #matches -h #threads -l timeLimit 
```

```
Examples:
./daf_10min -d yeast -q yeast_50n -n 100 -m 100000
./daf_10min -d yeast -q yeast_50n -n 100 -m 100000 -f
./daf_10min -d yeast -q yeast_50n -n 100 -m 100000 -c
./daf_10min -d yeast -q yeast_50n -n 100 -m 100000 -c -f

./daf_parallel_10min -d yeast -q yeast_50n -n 100 -m 100000 -h 16
./daf_parallel -d yeast -q yeast_50n -n 100 -m 100000 -h 16
```

## Results
DAF outputs following line per a query graph in a query set.
```
id  search_time  total_time  #recursive_calls  #found_matches  (0:solved, 1:unsolved)
```

```
Example Result: $ ./daf_10min -d yeast -q yeast_50n -n 100 -m 100000
...
0 15.5061 16.6469 37313 100005 0
1 0.052802 0.646388 38 302400 0
2 0.080433 3.41642 82 101700 0
3 1.89253 3.18895 1431 100100 0
4 0.808205 2.20496 1596 103488 0
5 4.8174 5.60534 8230 17640 0
6 0.083195 0.7845 114 108864 0
7 2.26529 3.03662 2309 100016 0
8 1.34417 2.13961 285 100061 0
9 2.24164 3.32384 721 100080 0
...
```

The number of found matches can be larger than the argument due to NEC technique.


## Input File Format

Data graph file format is a text format to store an undirected graph. 
- The first line of the file should be "t ID #vertices"
- Following lines of "v vertex-ID vertex-label" indicate the vertices in the graph.
- The vertices should be written in the file in ascending order of their IDs, and a vertex ID should be in [0, #vertices - 1].
- Following lines of "e vertex-ID1 vertex-ID2 edge-label" after the vertices indicate the undirected edges in the graph.

For example:
```
Line "t 1 3112" means the start of a graph with ID=1 and #vertices=3112.
Line "v 0 1" means there is a vertex with ID=0 and vertex-label=1 in the graph.
Line "e 1 133 0" meas there is an undirected edge between vertices with IDs 1 and 133, where edge label is 0.
```

Query graph file format is a text format to store undirected graphs.
- The first line of a graph should be "t ID #vertices #edges*2"
- Each of the following lines indicates "vertex ID, label of the vertex, degree of the vertex, adjacent vertices"
- The lines should be written in the file in ascending order of vertex IDs, and a vertex ID should be in [0, #vertices - 1].

For example:
```
Line "t 0 50 156" means the start of a graph with ID=0, #vertices=50, and #edges*2=156.
Line "1 6 5 0 2 3 6 11"  means v1 (vertex with ID=1) is labeled 6, and v1 is adjacent to five verteices; v0, v2, v3, v6, v11.
```
