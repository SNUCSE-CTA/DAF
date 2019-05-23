# DAF [SIGMOD 2019]
Efficient Subgraph Matching: Harmonizing Dynamic Programming, Adaptive Matching Order, and Failing Set Together

## Compile
To compile the query DAG construction program:
- Linux: type `./compile.sh` or `g++ -std=c++11 main.cpp dag.cpp -o program`

Binary files of DAF are available
- Linux: daf_1min
- Windows: daf_1min.exe

## Run

To run the query DAG construction program:
```
Usage: ./program datagraphFile querygraphFile #queries > outputFile
```

```
Example: ./program human human_10n 100 > human_10n.dag
```


To run DAF (time limit is set to 60 sec and #matches is set to 100,000)
```
Usages:

./daf_1min -d datagraphFile -q querygraphFile -a dagFile -n #queries // with input dagFile
or
./daf_1min -d datagraphFile -q querygraphFile -n # queries // query DAG is constructed by DAF
```

```
Examples:

./daf_1min -d human -q human_10n -a human_10n.dag -n 100
./daf_1min -d human -q human_10n -n 100
```

## Compare Results

To evaluate two results for a query set, use `sort_result.py` script as follows.

```
./daf_1min -d human -q human_10n -a human_10n.dag -n 100 > result1
./daf_1min -d human -q human_10n -n 100 > result2
python sort_result.py result1 result2
```
