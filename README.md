# DAF [SIGMOD 2019]
Efficient Subgraph Matching: Harmonizing Dynamic Programming, Adaptive Matching Order, and Failing Set Together

## Binary Files
Binary files of DAF are available for linux
- daf_10min: DAF for subgraph matching, which sets time limit of 10 minutes for each query graph in a query set.
- daf_big_10min: DAF especially for large data graph such as twitter, which sets the same time limit as `daf_10min`
- daf_parallel_unlimited: parallel version of DAF, which finds all embeddings without time limit.

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
```

```
Usages:
./daf_10min -d datagraphFile -q querygraphFile -n #queries -m #matches // path-size order and failing set are used by default.
./daf_10min -d datagraphFile -q querygraphFile -n #queries -m #matches -f // path-size order is used, and failing set is not used.
./daf_10min -d datagraphFile -q querygraphFile -n #queries -m #matches -c // candidate-size order and failing set are used
./daf_10min -d datagraphFile -q querygraphFile -n #queries -m #matches -c -f // candidate-size order is used, and failing set is not used.
```

```
Examples:
./daf_10min -d yeast -q yeast_50n -n 100 -m 100000
./daf_10min -d yeast -q yeast_50n -n 100 -m 100000 -f
./daf_10min -d yeast -q yeast_50n -n 100 -m 100000 -c
./daf_10min -d yeast -q yeast_50n -n 100 -m 100000 -c -f
```

## Results
DAF outputs following line per a query graph in a query set.
```
id  search_time  total_time  #recursive_calls  #found_matches  (0:solved, 1:unsolved)
```

```
Example Result:
0 0.093706 0.590738 92 15 0
1 0.038522 0.322034 37 2 0
2 0.074682 1.84184 82 2 0
3 0.321389 0.708628 771 11 0
4 0.035832 0.433856 49 2 0
5 0.043896 0.296473 93 2 0
6 0.035493 0.292103 70 2 0
7 0.030768 0.277881 66 2 0
8 0.037769 0.314199 51 7 0
9 0.075046 0.400043 168 2 0
...
```
