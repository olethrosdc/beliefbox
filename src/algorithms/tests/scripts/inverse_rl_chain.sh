baseres=$HOME/results/inverse_rl
n_runs=100

for epsilon in 0.01 0.1
for steps in 100 200 400 800 1600 3200 6400 12800 25600 51200 102400
do

    environment=Chain
    resdir=$baseres/$environment
    mkdir -p $resdir
    for n_states in 4 5 6 8
    do
        echo $environment $n_states $steps
        fname=$resdir/${n_states}s_${steps}T
        /usr/bin/time -o $fname.time ./bin/inverse_rl --algorithm Model --n_states $n_states --environment $environment --n_runs $n_runs --n_steps $steps | tee $fname.out | grep DV >${fname}_DV.out &
    done
    wait;
done

