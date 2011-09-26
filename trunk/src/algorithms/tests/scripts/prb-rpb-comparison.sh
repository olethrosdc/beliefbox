outdir=$HOME/results/inverse_rl/Chain/
runs=100
for samples in 2000 5000 10000
do
    for n_steps in 10 100 1000 10000
    do
        inverse_rl --algorithm Model --n_steps ${n_steps} --environment Chain --n_states 5 --n_chains 1 --iterations $samples --n_chain_samples $samples --MonteCarlo --n_runs $runs | tee $outdir/out_${n_steps}steps_${samples}samples &
    done
    wait;
done
