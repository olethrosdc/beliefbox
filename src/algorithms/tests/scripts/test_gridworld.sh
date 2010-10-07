maze=$HOME/projects/beliefbox/dat/maze1
width=8
height=8
T=10000
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
time ./bin/online_algorithms $S $A $gamma $lambda $rand $n_runs $n_episodes $T $model $Environment --maze_name=$maze --maze_height=$height --maze_width=$width  > $model.out&


model=Sarsa
echo $model
time ./bin/online_algorithms $S $A $gamma $lambda $rand $n_runs $n_episodes $T $model $Environment --maze_name=$maze --maze_height=$height --maze_width=$width  > $model.out&


model=Model
echo $model
time ./bin/online_algorithms $S $A $gamma $lambda $rand $n_runs $n_episodes $T $model $Environment --maze_name=$maze --maze_height=$height --maze_width=$width > $model.out&


model=Sampling
echo $model
time ./bin/online_algorithms $S $A $gamma $lambda $rand $n_runs $n_episodes $T $model $Environment --maze_name=$maze --maze_height=$height --maze_width=$width  > $model.out &

wait;

for model in QLearning Sarsa Model Sampling
do
    grep REWARD $model.out >$model.reward
    grep PAYOFF $model.out >$model.payoff
done
