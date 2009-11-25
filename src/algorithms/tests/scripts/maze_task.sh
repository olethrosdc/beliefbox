n_runs=100
n_episodes=10
n_steps=1000

gamma=0.9
lambda=0.9
random=0.01

n_actions=4
n_states=16

environment=Gridworld
resdir=~/experiments/aggregate/results/maze4_g${gamma}
mkdir -p $resdir
cd $resdir

exc=~/projects/beliefbox/src/algorithms/tests/bin/online_algorithms

echo "Run parameters" >run.params
echo  "exc n_states n_actions gamma lambda random n_runs n_episodes n_steps" >>run.params
echo  ${exc} $n_states $n_actions $gamma $lambda $random $n_runs $n_episodes $n_steps >>run.params

for model in Sarsa QLearning Model Collection
do
    #gdb --args ${exc} $n_states $n_actions $gamma $lambda $random $n_runs $n_episodes $n_steps ${model} $environment
    time ${exc} $n_states $n_actions $gamma $lambda $random $n_runs $n_episodes $n_steps ${model} $environment>${model}.out
    grep MSE ${model}.out >${model}.mse
    grep REWARD ${model}.out >${model}.reward
    grep PAYOFF ${model}.out >${model}.payoff
done
