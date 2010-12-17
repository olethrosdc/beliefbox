figure(1);
load Model.reward;
load Sarsa.reward;
load QLearning.reward; 
load Sampling.reward;

c=2;
plot((Sarsa(:,c)), ";Sarsa;",
     (QLearning(:,c)), ";Q-learning;",
     (Model(:,c)), ";Model;",
     (Sampling(:,c)), ";Sampling;");
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
