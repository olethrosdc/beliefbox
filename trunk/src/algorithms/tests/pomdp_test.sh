n_states=4
gamma=0.9
lambda=0.7
epsilon=0.01
n_runs=1000
T=10000000
Environment=POMDPGridworld
for model in QLearning HQLearning Sarsa Model;
do
    ./bin/pomdp_algorithms $n_states 3 $gamma $lambda $epsilon $n_runs 1 $T $model $Environment 4  | grep PAYOFF >$model.payoff
done

for depth in 1 2 3 4 8
do
    ./bin/pomdp_algorithms $n_states 3 $gamma $lambda $epsilon $n_runs 1 $T BVMM $Environment $depth  | grep PAYOFF >BVMM${depth}.payoff
done

octave -q ./pomdp_test_plot.m