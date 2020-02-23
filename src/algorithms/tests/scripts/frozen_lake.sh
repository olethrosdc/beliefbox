runs=100
outdir=frozen_lake_0
n_steps=100000
episode_steps=100
n_episodes=10000

mkdir -p $outdir

 for alg in QLearning Model 
 do
 	time ./bin/online_algorithms --n_episodes $n_episodes --episode_steps $episode_steps --n_steps $n_steps --n_runs $runs --algorithm $alg --environment Gridworld --maze_name ~/projects/beliefbox/dat/frozen_lake  --step_value 0 --pit_value 0 --goal_value 1 --randomness 0.33 >$outdir/$alg.out
 done

for samples in 1 2 4 8
do
	for alg in USampling LSampling
	do
		./bin/online_algorithms --n_episodes $n_episodes --episode_steps $episode_steps --n_steps $n_steps --n_runs $runs --algorithm $alg max_samples $samples --environment Gridworld --maze_name ~/projects/beliefbox/dat/frozen_lake  --step_value 0 --pit_value 0 --goal_value 1 --randomness 0.33 >$outdir/${alg}_${samples}s.out
	done
done

