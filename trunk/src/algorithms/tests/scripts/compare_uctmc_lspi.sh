#! /bin/bash

discount=0.99
rollouts=100
depth=1000
grid=20
horizon=1000

for environment in MountainCar  Pendulum PuddleWorld Bike CartPole Acrobot
do
	gdb --args ./bin/abc_lspi --environment $environment --discount $discount --depth_uct $depth --grid_uct $grid --n_rollouts_uct $rollouts --n_evaluations 10 --lambda_uct 1.0 >${environment}_ABC_UCT.out 
done
