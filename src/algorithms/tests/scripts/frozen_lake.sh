runs=1
outdir=frozen_lake_minus
n_steps=100000
episode_steps=1000
n_episodes=1000
samples=2
mkdir -p $outdir

for alg in PGAC #QLearning WQLearning Model Sarsa LSampling USampling UCRL TdBma LGBRL UGBRL TBRL PGAC
do
	time ./bin/online_algorithms --n_episodes $n_episodes --episode_steps $episode_steps --n_steps $n_steps --n_runs $runs --algorithm $alg --max_samples 2 --environment Gridworld --maze_name ~/projects/beliefbox/dat/frozen_lake  --step_value -0.1 --pit_value -1 --goal_value 1 --randomness 0.33 >$outdir/$alg.out
done


# for samples in 1 2 4 8
# do
# 	for alg in USampling LSampling
# 	do
# 		time ./bin/online_algorithms --n_episodes $n_episodes --episode_steps $episode_steps --n_steps $n_steps --n_runs $runs --algorithm $alg --max_samples $samples--environment Gridworld --maze_name ~/projects/beliefbox/dat/frozen_lake --step_value 0 --pit_value -0 --goal_value 1 --randomness 0.33 >$outdir/${alg}_${samples}s.out
# 	done
# done
