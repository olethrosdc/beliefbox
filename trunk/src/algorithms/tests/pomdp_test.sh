n_states=8
n_actions=3
gamma=0.99
lambda=0.7
epsilon=1.0
n_runs=10
n_episodes=10
T=1000000
Environment=MountainCar
outdir=$HOME/experiments/pomdp/$Environment
mkdir -p $outdir
for model in QLearning Sarsa USampling;
do
    ./bin/pomdp_algorithms $n_states $n_actions $gamma $lambda $epsilon $n_runs $n_episodes $T $model $Environment 4 > out
    grep PAYOFF out >$outdir/$model.payoff
    grep REWARD out >$outdir/$model.episode
done

for depth in 1 2 4
do
    ./bin/pomdp_algorithms $n_states $n_actions $gamma $lambda $epsilon $n_runs $n_episodes $T BVMM $Environment $depth >out
    grep PAYOFF out >$outdir/BVMM${depth}.payoff
    grep REWARD out >$outdir/BVMM${depth}.payoff
done

cp pomdp_test_plot.m $outdir/
cd $outdir
octave -q ./pomdp_test_plot.m
