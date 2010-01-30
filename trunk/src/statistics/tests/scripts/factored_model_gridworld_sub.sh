# shell for the job:
#PBS -S /bin/bash
# use one node with 2 cores:
#PBS -lnodes=1:ppn=4:cores8
# job requires at most 1 hours, 0 minutes
#     and 0 seconds wallclock time
#PBS -lwalltime=1:00:00
# cd to the directory where the program is to be called:

cd $HOME/experiments

script_dir=$HOME/projects/beliefbox/src/statistics/tests/scripts
bin_dir=$HOME/projects/beliefbox/src/statistics/tests/bin


for model in FMC BFMC BVMM CTW
do
    echo ./bin/factored_model Gridworld $model $iter $T $depth $maze_rand $action_rand $maze_dir/$maze_name $n_obs >${outdir}/${states}s_${depth}d_${model}.out >${out_dir}/${states}s_${depth}d_${model}.params
    time ./bin/factored_model Gridworld $model $iter $T $depth $maze_rand $action_rand $maze_dir/$maze_name $n_obs >${outdir}/${states}s_${depth}d_${model}.out

done

wait
