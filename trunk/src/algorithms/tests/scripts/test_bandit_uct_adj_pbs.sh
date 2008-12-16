#! /bin/bash
#

mkdir -p $HOME/experiments
cd $HOME/experiments

script_dir=$HOME/projects/beliefbox/src/algorithms/tests/scripts

res_dir=$HOME/results/belief_uct_adj/bandit
mkdir -p $res_dir

experiments=1000

for gamma in  0.999 0.9999 #0.5 0.9 0.99 #0.999 0.9999 0.99999
  do
  for iter in  0 1 2 3 4 5 6 7 8 10 12 14 16 18 22 28 32 64
    do
    for method in 0 1 2 14 8 10 11 13 #0 1 2 3 4 5 6 7 8 9 10 11 12 13 14
      do
      results=$res_dir
      qsub -v"method=$method","iter=$iter","gamma=$gamma","experiments=$experiments","results=$results"  $script_dir/test_bandit_uct_adj_sub.sh
    done
  done
done