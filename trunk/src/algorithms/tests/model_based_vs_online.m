figure(1);
load Model.reward;
load Sarsa.reward;
load QLearning.reward; 
load Sampling1.reward;
load Sampling2.reward;
load Sampling4.reward;
load Sampling8.reward;

c=1;
K=10;
plot(moving_average(Sarsa(:,c), K), ";Sarsa;",
     moving_average(QLearning(:,c), K), ";Q-learning;",
     moving_average(Model(:,c), K), ";Model;",
     moving_average(Sampling1(:,c), K), ";Sampling1;"
,     moving_average(Sampling2(:,c), K), ";Sampling2;");
     moving_average(Sampling4(:,c), K), ";Sampling4;",
     moving_average(Sampling8(:,c), K), ";Sampling8;")
xlabel("Episode");
ylabel("U");


# figure(2);

# load Model.payoff;
# load Sarsa.payoff;
# load QLearning.payoff; 
# load Sampling.payoff;

# c=1;
# plot(cumsum(Sarsa(:,c)), ";Sarsa;",
#      cumsum(QLearning(:,c)), ";Q-learning;",
#      cumsum(Model(:,c)), ";Model;",
#      cumsum(Sampling(:,c)), ";Sampling;");
# xlabel("T");
# ylabel("R");
