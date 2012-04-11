n_runs=10
n_episodes=100
n_steps=100000

gamma=0.99
lambda=0.9
random=0.01

n_actions=4
n_states=16

environment=Gridworld
maze=$HOME/projects/beliefbox/dat/maze2
resdir=~/results/maze2_g${gamma}
mkdir -p $resdir
cd $resdir

exc=~/projects/beliefbox/src/algorithms/tests/bin/online_algorithms

echo "Run parameters" >run.params
echo  "exc n_states n_actions gamma lambda random n_runs n_episodes n_steps" >>run.params

for model in Sarsa QLearning #Model 
do
    #gdb --args ${exc} $n_states $n_actions $gamma $lambda $random $n_runs $n_episodes $n_steps ${model} $environment
    time ${exc} --gamma $gamma --lambda $lambda --randomness $random --n_runs $n_runs --n_episodes $n_episodes --n_steps $n_steps --algorithm ${model} --environment $environment --maze_name $maze >${model}.out
    grep MSE ${model}.out >${model}.mse
    grep REWARD ${model}.out >${model}.reward
    grep PAYOFF ${model}.out >${model}.payoff
done


for model in Sampling
do
    for s in 1 2 4 8
    do
        for bound in lower upper
        do
            time ${exc} --gamma $gamma --lambda $lambda --randomness $random --n_runs $n_runs --n_episodes $n_episodes --n_steps $n_steps --algorithm ${model} --environment $environment --reward_prior Fixed --maze_name $maze >${model}${s}_${bound}.out &
            
            time ${exc} --gamma $gamma --lambda $lambda --randomness $random --n_runs $n_runs --n_episodes $n_episodes --n_steps $n_steps --algorithm ${model} --environment $environment --maze_name --maze>${model}${s}_${bound}.out &
        done
        wait
        
        for bound in lower upper
        do
            grep MSE ${model}${s}_${bound}.out >${model}${s}_${bound}.mse
            grep REWARD ${model}${s}_${bound}.out >${model}${s}_${bound}.reward
            grep PAYOFF ${model}${s}_${bound}.out >${model}${s}_${bound}.payoff
        done
    done
done


