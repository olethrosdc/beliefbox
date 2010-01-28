#! /bin/bash

mkdir -p $HOME/experiments
cd $HOME/experiments

script_dir=$HOME/projects/beliefbox/src/statistics/tests/scripts
environment=1DMaze

for p in 0.0 #0.01 0.1
do
    for ap in 1.0 #0.5 0.1
    do 
	out_dir=$HOME/results/bvmdp/${environment}/${p}p_${ap}ap/
	for states in 4 #8 10 12 16 32
	do
            for depth in 1 #2 4 8 10 12
            do
		echo "./bin/factored_model 1DMaze $model $iter $T $depth $p 1.0 $states >${outdir}/${states}s_${depth}d_${model}.out" >>$out_dir/params.txt
		qsub -v"environment=$environment","iter=$iter","T=$T","depth=$depth","p=$p","pa=$pa","states=$states","out_dir=$out_dir" $script_dir/factored_model_sub.sh
            done
        done
    done
done
