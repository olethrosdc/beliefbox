# shell for the job:
#PBS -S /bin/bash
# use one node with 2 cores:
#PBS -lnodes=1:ppn=2:cores2
# job requires at most 10 hours, 0 minutes
#     and 0 seconds wallclock time
#PBS -lwalltime=01:00:00
# cd to the directory where the program is to be called:

cd $HOME/projects/beliefbox/src/algorithms/tests

echo "method: $method gamma: $gamma actions: $actions experiments $experiments"

# rm -f $results
for actions in 2 4
do
    time ./bin/bandit_uct $method $gamma $iter $actions 0 $experiments  >>a${actions}_$results &
done

wait
