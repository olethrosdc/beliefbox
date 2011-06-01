for samples in 1 10 20 50 100 200 500 1000 2000 5000 10000
do
    for n_steps in 10 100 1000 10000
    do
        ./bin/inverse_rl --algorithm Model --n_steps ${n_steps} --environment Chain --n_states 5 --n_chains 1 --n_chain_samples $samples --MonteCarlo --n_runs 10 | tee out_${n_steps}steps_${samples}samples &;
    done
    wait;
done
