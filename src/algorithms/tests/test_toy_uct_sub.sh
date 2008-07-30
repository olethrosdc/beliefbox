# shell for the job:
#PBS -S /bin/bash
# use one node with 8 cores:
#PBS -lnodes=1:ppn=8:cores8
# job requires at most 10 hours, 0 minutes
#     and 0 seconds wallclock time
#PBS -lwalltime=76:00:00
# cd to the directory where the program is to be called:


cd $HOME/projects/beliefbox/src/algorithms/tests

rm -f $results
for gamma in  0.1 0.5 0.9 0.99 0.999 0.9999 0.99999 0.999999
do
  time ./bin/toy_uct_stopping_problem $gamma $iter 0 10000 >>$results &
done

wait
