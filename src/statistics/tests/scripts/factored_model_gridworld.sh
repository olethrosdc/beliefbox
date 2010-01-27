#! /bin/sh

outdir=~/experiments/bvmdp/gridworld/maze7/p001_policy/pow_prior_05
mkdir -p $outdir
T=100000
iter=100
states=0
echo "factored_model $iter $T" >${outdir}/run.params

for action_rand in 0 0.01 0.1 0.5 1.0
  do
  for maze_rand in 0 0.01 0.1
    do
    for depth in 0 1 2 4 8
      do
      for model in FMC BFMC BVMM CTW
        do
        echo $model $depth
        time ./bin/factored_model $iter $T $states $depth $model >${outdir}/${states}s_${depth}d_${model}.out
      done
    done
  done
done
