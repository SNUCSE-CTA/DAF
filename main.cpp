// Name     : main.cpp
// Author   : Geonmo Gu
// Version  : 1.000
// Copyright: Apache License 2.0

#include "dag.h"

using namespace std;

//Usage: ./program dataGraph queryGraph numQuery
int main(int argc, char *argv[]) 
{

    if(argc != 4) {
        cout << "Usage: ./program dataFile queryFile numQuery" << endl;
        return -1;
    }
    //cout << "Data File : " << argv[1] << endl;
    //cout << "Query file: " << argv[2] << endl;
    //cout << "Num Query: " << argv[3] << endl;

    int numQuery = atoi(argv[3]);

    readDataGraph(argv[1]);

    ifstream queryFile(argv[2]);
    for(int i = 0; i < numQuery; ++i) {
        char tag;
        int id;
        int num;
        int sumDegree;
        queryFile >> tag >> id >> num >> sumDegree;
        numQueryNode = num;

        readQueryGraph(queryFile, sumDegree);
        buildDAG();
    }
    queryFile.close();

    clearMemory();
    return 0;
}
