#! /bin/bash

discount=0.99
rollouts=100
depth=1000
grid=20
horizon=1000

for environment in Pendulum PuddleWorld Bike CartPole Acrobot
do
	./bin/UCTMC_test --environment $environment --discount $discount --depth $depth --grids $grid --horizon $horizon --n_rollouts $rollouts --n_evaluations 10 --n_episodes 100 --lambda 1.0 #>${environment}_UCTMC.out 
done
