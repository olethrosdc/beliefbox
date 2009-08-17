#! /bin/bash
alpha=0.5
knn=5
rbf=.1
../bin/knn_test $alpha $knn $rbf | tee out
grep V out >Vout
grep Q0 out >Q0out
grep Q1 out >Q1out
grep Q2 out >Q2out

octave ./plot_pendulum.m
gv V_pendulum.eps&
gv Q0_pendulum.eps&
gv Q1_pendulum.eps&
gv Q2_pendulum.eps&
gv DQ02_pendulum.eps&


