#! /bin/bash



EXEC=/projects/beliefbox/src/algorithms/tests/bin/context_bandit_uct
method=0
gamma=0.999;
n_iter=0
n_actions=2
verbose=0
n_experiments=1000
for T in 1 2 4 8 16 32 64 128 256 512 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 12000 14000 16000 18000 20000
do
    ${EXEC} $method $gamma $n_iter $n_actions $verbose $n_experiments $T >> context_m{$method}_i${n_iter}
done
