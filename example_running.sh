#compute linear orders for 100 query graphs
echo -ne 'computing linear orders for 100 query graphs...'
./program yeast yeast_400n 100 > yeast_400n.dag
echo 'done'

#run daf using the computed DAGs
echo -ne 'running DAF using the computed DAGs...'
./daf_2min -d yeast -q yeast_400n -a yeast_400n.dag -n 100 -m 100000 > result_dag
echo 'done'

#run original daf
echo -ne 'running the original DAF...'
./daf_2min -d yeast -q yeast_400n -n 100 -m 100000 > result_daf
echo 'done'

#compare two results
echo '*comparing result*'
python sort_result.py result_dag result_daf
