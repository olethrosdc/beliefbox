n_runs=10
n_episodes=10
n_steps=1000

gamma=1.0
lambda=0.9
random=0.01

n_actions=4
n_states=16


resdir=~/experiments/aggregate/s${n_states}_a${n_actions}_g${gamma}
mkdir -p $resdir
cd $resdir

exc=~/projects/beliefbox/src/algorithms/tests/bin/online_algorithms

echo "Run parameters" >run.params
echo  "exc n_states n_actions gamma lambda random n_runs n_episodes n_steps" >>run.params
echo  ${exc} $n_states $n_actions $gamma $lambda $random $n_runs $n_episodes $n_steps >>run.params

time ${exc} $n_states $n_actions $gamma $lambda $random $n_runs $n_episodes $n_steps Sarsa >sarsa.out
grep MSE sarsa.out >sarsa.mse
grep REWARD sarsa.out >sarsa.reward
grep PAYOFF sarsa.out >sarsa.payoff

time ${exc} $n_states $n_actions $gamma $lambda $random $n_runs $n_episodes $n_steps QLearning >qlearning.out
grep MSE qlearning.out >qlearning.mse
grep REWARD qlearning.out >qlearning.reward
grep PAYOFF qlearning.out >qlearning.payoff

time ${exc} $n_states $n_actions $gamma $lambda $random $n_runs $n_episodes $n_steps Model >model.out
grep MSE model.out >model.mse
grep REWARD model.out >model.reward
grep PAYOFF model.out >model.payoff

time ${exc} $n_states $n_actions $gamma $lambda $random $n_runs $n_episodes $n_steps Collection >collection.out
grep MSE collection.out >collection.mse
grep REWARD collection.out >collection.reward
grep PAYOFF collection.out >collection.payoff
