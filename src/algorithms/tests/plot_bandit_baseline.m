## -*- Mode: octave -*-

n_arms=2;

iter = 5;

X=load ("results/belief_uct/bandit/baseline.out");
Y=load ("results/belief_uct/bandit/uct_complete.out");
hold off
s=X(:,1)==0 & X(:,3)==n_arms;
plot(2./(1-X(s,2)), X(s,6) - X(s,4), '0@-;ucb;');
hold on

s=X(:,1)==3 & X(:,3)==n_arms;
plot(2./(1-X(s,2)), X(s,6) - X(s,4), '1-;biased beta;');

s=Y(:,1)==0 & Y(:,2)==iter & Y(:,4)==n_arms;
loglog(2./(1-Y(s,3)), Y(s,7) - Y(s,5), '2@-;serial;');


s=Y(:,1)==1 & Y(:,2)==iter & Y(:,4)==n_arms;
loglog(2./(1-Y(s,3)), Y(s,7) - Y(s,5), '3@-;random;');
s=Y(:,1)==2 & Y(:,2)==iter & Y(:,4)==n_arms;
loglog(2./(1-Y(s,3)), Y(s,7) - Y(s,5), '4@-;LB1;');

s=Y(:,1)==3 & Y(:,2)==iter & Y(:,4)==n_arms;
loglog(2./(1-Y(s,3)), Y(s,7) - Y(s,5), '5@-;LB2;');
grid on;

title(strcat("n_{arms} = ", num2str(n_arms),
             ", n_{iter} = ", num2str(iter)));
print("bandit_baseline.eps");
