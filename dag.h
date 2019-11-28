// Name     : dag.h
// Author   : Geonmo Gu
// Version  : 1.000
// Copyright: Apache License 2.0

#ifndef DAG_H
#define DAG_H

#include <iostream>
#include <fstream>
#include <algorithm> //sort
#include <sstream>
#include <string>
#include <cstring> //memset
#include <set>
#include <cfloat> //DBL_MAX

using namespace std;

void readDataGraph(string);
void readQueryGraph(ifstream&, int);
int selectRoot();
void buildDAG();
void printDAG();
void clearMemory();
bool sortByDegreeData(int, int);
bool sortByLabel(int, int);
bool sortByDegreeQuery(int, int);
bool sortByLabelFreqQuery(int, int);
int binaryLowerBound(int, int, int);

extern int numQueryNode;

#endif
