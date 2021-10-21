# Name     : sort_result.py
# Author   : Geonmo Gu
# Version  : 1.000
# Copyright: Apache License 2.0


import sys

#Usage: python program.py file1 file2

if len(sys.argv) != 3:
	print('Usage: python program.py result1 result2')
	sys.exit(-1)

infile1 = open( sys.argv[1], 'r')
infile2 = open( sys.argv[2], 'r')

totalTimeFile1 = {}
totalTimeFile2 = {}
totalRecFile1 = {}
totalRecFile2 = {}
numSolvedFile1 = 0
numSolvedFile2 = 0

################## Read File 1
for line in infile1.readlines():
	args = line.split()
	tag = args[0]

	if str.isdigit(tag):
		queryId = int(tag)
		totalTime = float(args[2])
		recursion = float(args[3])
		solved = int(args[5])

		if solved == 0: #solved
			totalTimeFile1[queryId] = totalTime
			totalRecFile1[queryId] = recursion
			numSolvedFile1 += 1

infile1.close()

################## Read File 2
for line in infile2.readlines():
	args = line.split()
	tag = args[0]
	if str.isdigit(tag):
		queryId = int(tag)
		totalTime = float(args[2])
		recursion = float(args[3])
		solved = int(args[5])

		if solved == 0: #solved
			totalTimeFile2[queryId] = totalTime
			totalRecFile2[queryId] = recursion
			numSolvedFile2 += 1

infile2.close()

minNumSolved = min( numSolvedFile1, numSolvedFile2 )

################## Compute Average Values for File 1
numAdded = 0
avgTime = 0
avgRec = 0
for queryId, totalTime in sorted( totalTimeFile1.iteritems(), key = lambda (k, v): (v, k) ):
	if numAdded >= minNumSolved:
		break
	
	numAdded += 1
	avgTime += totalTimeFile1[queryId]
	avgRec += totalRecFile1[queryId]

if numAdded != 0:
	avgTime /= numAdded
	avgRec /= numAdded

print(sys.argv[1])
print('Average elapsed time: ' + str(avgTime))
print('Average #recursive calls: ' + str(avgRec))
print('#Solved queries: ' + str(numSolvedFile1))

print('')
################## Compute Average Values for File 2
numAdded = 0
avgTime = 0
avgRec = 0
for queryId, totalTime in sorted( totalTimeFile2.iteritems(), key = lambda (k, v): (v, k) ):
	if numAdded >= minNumSolved:
		break
	
	numAdded += 1
	avgTime += totalTimeFile2[queryId]
	avgRec += totalRecFile2[queryId]

if numAdded != 0:
	avgTime /= numAdded
	avgRec /= numAdded

print(sys.argv[2])
print('Average elapsed time: ' + str(avgTime))
print('Average #recursive calls: ' + str(avgRec))
print('#Solved queries: ' + str(numSolvedFile2))
