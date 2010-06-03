n_runs=10000
n_episodes=1
n_steps=100000

gamma=0.75
lambda=0.5
random=0.01

n_actions=2
for n_states in 4 8 6
do
#environment=Gridworld
#resdir=~/experiments/aggregate/results/maze2_g${gamma}
    environment=OneDMaze
    resdir=~/experiments/pomdp/results/1dmaze_${n_states}_g${gamma}
#environment=ContextBandit
#resdir=~/experiments/aggregate/results/context
    mkdir -p $resdir
    cd $resdir

    exc=~/projects/beliefbox/src/algorithms/tests/bin/pomdp_algorithms

    echo "Run parameters" >run.params
    echo  "exc n_states n_actions gamma lambda random n_runs n_episodes n_steps" >>run.params
    echo  ${exc} $n_states $n_actions $gamma $lambda $random $n_runs $n_episodes $n_steps >>run.params

    for model in Sarsa QLearning Model
    do
        time ${exc} $n_states $n_actions $gamma $lambda $random $n_runs $n_episodes $n_steps ${model} $environment >${model}.out
        grep REWARD ${model}.out >${model}.reward
        grep PAYOFF ${model}.out >${model}.payoff
    done
done