# shell for the job:
#PBS -S /bin/bash
# use one node with 2 cores:
#PBS -lnodes=1:ppn=2:cores2
# job requires at most 10 hours, 0 minutes
#     and 0 seconds wallclock time
#PBS -lwalltime=72:00:00
# cd to the directory where the program is to be called:

cd $HOME/experiments

script_dir=$HOME/projects/beliefbox/src/algorithms/tests/scripts
bin_dir=$HOME/projects/beliefbox/src/algorithms/tests/bin

echo "results dir: $results method: $method gamma: $gamma actions: $actions experiments $experiments"

# rm -f $results
for actions in 2 #4
do
    time $bin_dir/bandit_uct_adjusted $method $gamma $iter $actions 0 $experiments  >>$results/a${actions}_g${gamma}_m${method}_i${iter}_uct.out &
done

wait
