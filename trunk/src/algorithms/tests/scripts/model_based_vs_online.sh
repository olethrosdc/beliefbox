n_runs=10000
n_episodes=10
n_steps=10000

gamma=0.999
lambda=0.9
random=0.01

n_actions=4
n_states=16

#environment=Gridworld
#resdir=~/experiments/aggregate/results/maze2_g${gamma}
environment=ContextBandit
resdir=~/experiments/aggregate/results/context
mkdir -p $resdir
cd $resdir

exc=~/projects/beliefbox/src/algorithms/tests/bin/online_algorithms

echo "Run parameters" >run.params
echo  "exc n_states n_actions gamma lambda random n_runs n_episodes n_steps" >>run.params
echo  ${exc} $n_states $n_actions $gamma $lambda $random $n_runs $n_episodes $n_steps >>run.params

time ${exc} $n_states $n_actions $gamma $lambda $random $n_runs $n_episodes $n_steps Sarsa $environment >sarsa.out
grep MSE sarsa.out >sarsa.mse
grep REWARD sarsa.out >sarsa.reward
grep PAYOFF sarsa.out >sarsa.payoff

time ${exc} $n_states $n_actions $gamma $lambda $random $n_runs $n_episodes $n_steps QLearning $environment >qlearning.out
grep MSE qlearning.out >qlearning.mse
grep REWARD qlearning.out >qlearning.reward
grep PAYOFF qlearning.out >qlearning.payoff

#time ${exc} $n_states $n_actions $gamma $lambda $random $n_runs $n_episodes $n_steps Model $environment >model.out
#grep MSE model.out >model.mse
#grep REWARD model.out >model.reward
#grep PAYOFF model.out >model.payoff

time ${exc} $n_states $n_actions $gamma $lambda $random $n_runs $n_episodes $n_steps ContextBanditGaussian $environment >model.out
grep MSE model.out >model.mse
grep REWARD model.out >model.reward
grep PAYOFF model.out >model.payoff

time ${exc} $n_states $n_actions $gamma $lambda $random $n_runs $n_episodes $n_#steps Aggregate $environment >aggregate.out
grep MSE aggregate.out >aggregate.mse
grep REWARD aggregate.out >aggregate.reward
grep PAYOFF aggregate.out >aggregate.payoff



time ${exc} $n_states $n_actions $gamma $lambda $random $n_runs $n_episodes $n_steps ContextBanditCollection $environment >collection.out
grep MSE collection.out >collection.mse
grep REWARD collection.out >collection.reward
grep PAYOFF collection.out >collection.payoff
