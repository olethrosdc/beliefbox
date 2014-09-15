#! /bin/bash

discount=0.99
rollouts=1
depth=100
grid=10
horizon=1000

for environment in MountainCar Pendumlum Puddle Bicycle CartPole Acrobot
do
	./bin/UCTMC_test --environment $environment --discount $discount --depth $depth --grids $grid --horizon $horizon --n_rollouts $rollouts --n_evaluations 10 --n_episodes 100 >${environment}_UCTMC.out 
done
