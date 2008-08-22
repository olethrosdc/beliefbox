## -*- Mode: octave -*-

n_arms=2;

X=load ("results/belief_uct/bandit/baseline.out");
hold off
s=X(:,1)==0 & X(:,3)==n_arms;
plot(1./(1-X(s,2)), X(s,6) - X(s,4), '0@-;ucb;');
hold on
s=X(:,1)==2 & X(:,3)==n_arms;
plot(1./(1-X(s,2)), X(s,6) - X(s,4), '1-;beta;');

s=X(:,1)==3 & X(:,3)==n_arms;
plot(1./(1-X(s,2)), X(s,6) - X(s,4), '2-;biased beta;');

s=X(:,1)==4 & X(:,3)==n_arms;
plot(1./(1-X(s,2)), X(s,6) - X(s,4), '3-;greedy;');

s=X(:,1)==5 & X(:,3)==n_arms;
plot(1./(1-X(s,2)), X(s,6) - X(s,4), '4-;\epsilon-greedy 0.01;');

s=X(:,1)==6 & X(:,3)==n_arms;
loglog(1./(1-X(s,2)), X(s,6) - X(s,4), '5-;\epsilon-greedy 0.1;');

grid on;
title(strcat("n_{arms} = ",
             num2str(n_arms));
print("bandit_baseline.eps");
