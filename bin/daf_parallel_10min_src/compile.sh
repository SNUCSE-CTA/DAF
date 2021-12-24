sdsl_inc=/home/gmgu/include
sdsl_lib=/home/gmgu/lib
# Sequential version
#g++ -pthread -w -std=c++11 -O3 -DNDEBUG -I /home/hjkim/include -L /home/hjkim/lib DAF.cpp -o daf_10min -lsdsl -ldivsufsort -ldivsufsort64

# Parallelized version
#g++ -pthread -w -std=c++11 -O3 -DNDEBUG -I /home/hjkim/include -L /home/hjkim/lib DAF.cpp -fopenmp -o daf_parallel_10min -lsdsl -ldivsufsort -ldivsufsort64 -DPARALLEL

# Parallelized version (finding all embeddings)
#g++ -pthread -w -std=c++11 -O3 -DNDEBUG -I /home/hjkim/include -L /home/hjkim/lib DAF.cpp -fopenmp -o daf_parallel_unlimited -lsdsl -ldivsufsort -ldivsufsort64 -DPARALLEL -DUNLIMITED

# Parallelized version (+ reduce reading the number of matches)
#g++ -pthread -w -std=c++11 -O3 -DNDEBUG -I /home/hjkim/include -L /home/hjkim/lib DAF.cpp -fopenmp -o daf_10min_parallel_opt -lsdsl -ldivsufsort -ldivsufsort64 -DPARALLEL -DREDUCE_COUNT_READING

# Parallelized version for the experiment
#g++ -pthread -w -std=c++11 -O3 -DNDEBUG -I /home/hjkim/include -L /home/hjkim/lib DAF.cpp -fopenmp -o daf_10min_parallel_exp -lsdsl -ldivsufsort -ldivsufsort64 -DPARALLEL -DPARALLEL_EXPERIMENT

# Parallelized version for the experiment with local limit per region
#g++ -pthread -w -std=c++11 -O3 -DNDEBUG -I /home/hjkim/include -L /home/hjkim/lib DAF.cpp -fopenmp -o daf_10min_parallel_exp_local_limit -lsdsl -ldivsufsort -ldivsufsort64 -DPARALLEL -DPARALLEL_EXPERIMENT -DLOCAL_LIMIT_PER_REGION


#g++ -pthread -w -std=c++11 -O3 -DNDEBUG DAF.cpp -fopenmp -o daf_parallel_10min -DPARALLEL

g++ -pthread -w -std=c++11 -O3 -DNDEBUG -I ${sdsl_inc} -L ${sdsl_lib} DAF.cpp -fopenmp -o daf_parallel_10min -DPARALLEL -lsdsl -ldivsufsort -ldivsufsort64
