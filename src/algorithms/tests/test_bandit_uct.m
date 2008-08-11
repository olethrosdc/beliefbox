## -*- Mode: octave -*-

lw=4;

X=load("results/belief_uct/bandit/uct.out");
Y=load("results/belief_uct/bandit/ucb.out");
gamma = 0.99;
n_arms=2;
hold off;
idx = Y(:,2)==gamma & Y(:,1)==0 & Y(:,3)== n_arms; plot([1 22], [1 1]*(Y(idx,6)-Y(idx,4)), "0@-;ucb;", "linewidth", lw);
hold on;
uct=5;
oracle=7;
sel_exp = X(:,3)==gamma & X(:,4)==n_arms;

idx = sel_exp & X(:,1)==0; plot(X(idx,2), X(idx,oracle)-X(idx,uct), "1-;serial;", "linewidth", lw);
idx = sel_exp & X(:,1)==1; plot(X(idx,2), X(idx,oracle)-X(idx,uct), "1@-;random;", "linewidth", lw);
idx = sel_exp & X(:,1)==2; plot(X(idx,2), X(idx,oracle)-X(idx,uct), "2-;mean;", "linewidth", lw);
idx = sel_exp & X(:,1)==3; plot(X(idx,2), X(idx,oracle)-X(idx,uct), "2@-;mean gamma;", "linewidth", lw);
idx = sel_exp & X(:,1)==4; plot(X(idx,2), X(idx,oracle)-X(idx,uct), "3-;thompson;", "linewidth", lw);
idx = sel_exp & X(:,1)==5; plot(X(idx,2), X(idx,oracle)-X(idx,uct), "3@-;thompson gamma;", "linewidth", lw);
idx = sel_exp & X(:,1)==12; plot(X(idx,2), X(idx,oracle)-X(idx,uct), "5@-;high prob;", "linewidth", lw);
idx = sel_exp & X(:,1)==9; plot(X(idx,2), X(idx,oracle)-X(idx,uct), "4@-;high prob gamma;", "linewidth", lw);

xlabel("number of expansions");
ylabel("total regret");
legend("location", "northeast");
legend("boxon");
grid on
print("test_bandit_uct.eps");
