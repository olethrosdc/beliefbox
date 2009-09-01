# shell for the job:
#PBS -S /bin/bash
# use one node with 8 cores:
#PBS -lnodes=1:ppn=8:cores8
# job requires at most 10 hours, 0 minutes
#     and 0 seconds wallclock time
#PBS -lwalltime=96:00:00
# cd to the directory where the program is to be called:

cd $HOME/experiments

script_dir=$HOME/projects/beliefbox/src/algorithms/tests/scripts
bin_dir=$HOME/projects/beliefbox/src/algorithms/tests/bin

echo "results dir: $results method: $method gamma: $gamma actions: $actions experiments $experiments"

# rm -f $results
for method in 0 1 3 4 5 13 14 15
do
    time $bin_dir/bandit_uct $method $gamma $iter $actions 0 $experiments  >>$results/a${actions}_g${gamma}_m${method}_i${iter}_uct.out &
done

wait
