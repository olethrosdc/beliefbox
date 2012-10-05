## -*- Mode: octave -*-

load QLearning.payoff;
##load HQLearning.payoff;
load Sarsa.payoff;
load Model.payoff;
load BVMM1.payoff;
load BVMM2.payoff;
load BVMM3.payoff;
load BVMM4.payoff;
load BVMM8.payoff;
load USampling.payoff;


K=1000;
figure(1);
plot(moving_average(QLearning, K), ";Q-Learning;",
#     moving_average(Sarsa, K), ";Sarsa;",
#     moving_average(HQLearning, K), ";HQ-Learning;",
     moving_average(USampling, K), ";USampling;",
     moving_average(Model, K), ";Model;",
#     moving_average(BVMM1, K), ";BVMM1;",
     moving_average(BVMM2, K), ";BVMM2;",
#     moving_average(BVMM3, K), ";BVMM3;",
     moving_average(BVMM4, K), ";BVMM4;");
#     moving_average(BVMM8, K), ";BVMM8;");
print("pomdp_test_reward.eps", "-dashed", "-color");

cS = cumsum(Sarsa);
figure(2);
plot(cS - cumsum(QLearning), ";Q-Learning;",
#     cS - cumsum(HQLearning), ";HQ-Learning;",
     cS - cumsum(USampling), ";USampling;",
     cS - cumsum(Model), ";Model;",
#     cS - cumsum(BVMM1), ";BVMM1;",
     cS - cumsum(BVMM2), ";BVMM2;",
#     cS - cumsum(BVMM3), ";BVMM3;",
     cS - cumsum(BVMM4), ";BVMM4;");


print("pomdp_test_regret.eps", "-dashed", "-color");

