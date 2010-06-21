## -*- Mode: octave -*-

load QLearning.payoff;
load HQLearning.payoff;
load Sarsa.payoff;
load Model.payoff;
load BVMM1.payoff;
load BVMM2.payoff;
load BVMM3.payoff;
load BVMM4.payoff;
load BVMM8.payoff;

subplot(2,1,1);
K=1000;
plot(moving_average(QLearning, K), ";Q-Learning;",
#     moving_average(Sarsa, K), ";Sarsa;",
     moving_average(HQLearning, K), ";HQ-Learning;",
     moving_average(Model, K), ";Model;",
#     moving_average(BVMM1, K), ";BVMM1;",
     moving_average(BVMM2, K), ";BVMM2;",
#     moving_average(BVMM3, K), ";BVMM3;",
     moving_average(BVMM4, K), ";BVMM4;",
     moving_average(BVMM8, K), ";BVMM8;");

subplot(2,1,2);
plot(cumsum(QLearning), ";Q-Learning;",
#     cumsum(Sarsa), ";Sarsa;",
     cumsum(HQLearning), ";HQ-Learning;",
     cumsum(Model), ";Model;",
#     cumsum(BVMM1), ";BVMM1;",
     cumsum(BVMM2), ";BVMM2;",
#     cumsum(BVMM3), ";BVMM3;",
     cumsum(BVMM4), ";BVMM4;",
     cumsum(BVMM8), ";BVMM8;");



     