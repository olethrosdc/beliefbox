maze=$HOME/projects/beliefbox/dat/maze1
width=8
height=8
T=1000
n_runs=100
n_episodes=1000
rand=0.01
lambda=0.9
gamma=0.95
model=QLearning
S=8
A=3
Environment=Gridworld

model=QLearning
echo $model
time ./bin/online_algorithms --gamma $gamma --lambda $lambda --n_runs $n_runs --n_episodes $n_episodes --n_steps $T --maze_name=$maze --maze_height=$height --maze_width=$width --algorithm $model --environment $Environment > $model.out

model=Sarsa
echo $model
time ./bin/online_algorithms --gamma $gamma --lambda $lambda --n_runs $n_runs --n_episodes $n_episodes --n_steps $T --maze_name=$maze --maze_height=$height --maze_width=$width --algorithm $model --environment $Environment > $model.out



model=Model
echo $model
time ./bin/online_algorithms --gamma $gamma --lambda $lambda --n_runs $n_runs --n_episodes $n_episodes --n_steps $T --maze_name=$maze --maze_height=$height --maze_width=$width --algorithm $model --environment $Environment > $model.out


model=Sampling
echo $model
time ./bin/online_algorithms --gamma $gamma --lambda $lambda --n_runs $n_runs --n_episodes $n_episodes --n_steps $T --maze_name=$maze --maze_height=$height --maze_width=$width --algorithm $model --environment $Environment > $model.out

wait;

for model in QLearning Sarsa Model Sampling
do
    grep REWARD $model.out >$model.reward
    grep PAYOFF $model.out >$model.payoff
done
