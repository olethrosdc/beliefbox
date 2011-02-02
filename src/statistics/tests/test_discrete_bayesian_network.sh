#! /bin/bash

./bin/discrete_bayesian_network >out
grep DATA out >samples
grep JOINT out >joint

octave ./test_discrete_bayesian_network.m
