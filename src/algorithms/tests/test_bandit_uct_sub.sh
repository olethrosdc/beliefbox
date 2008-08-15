# shell for the job:
#PBS -S /bin/bash
# use one node with 8 cores:
#PBS -lnodes=1:ppn=8:cores8
# job requires at most 10 hours, 0 minutes
#     and 0 seconds wallclock time
#PBS -lwalltime=96:00:00
# cd to the directory where the program is to be called:


cd $HOME/projects/beliefbox/src/algorithms/tests

rm -f $results
for gamma in  0.5 0.9 0.99 0.999 0.9999 0.99999
do
    time ./bin/bandit_uct $method $gamma $iter $actions 0 $experiments  >>$results &
done

wait
