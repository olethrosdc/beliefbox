## -*- Mode: octave -*-

X=load("results/belief_uct/bandit/uct.out");
Y=load("results/belief_uct/bandit/ucb.out");
gamma = 0.999;
column = 5;
hold off;
idx = Y(:,2)==gamma & Y(:,1)==0 & Y(:,3)==2; plot([1 22], [1 1]*Y(idx,4), '0@-;ucb;');
hold on;
idx = X(:,3)==gamma & X(:,1)==0; plot(X(idx,2), X(idx, column), '1-;serial;');
idx = X(:,3)==gamma & X(:,1)==1; plot(X(idx,2), X(idx, column), '1@-;random;');
idx = X(:,3)==gamma & X(:,1)==2; plot(X(idx,2), X(idx, column), '2-;mean;');
idx = X(:,3)==gamma & X(:,1)==3; plot(X(idx,2), X(idx, column), '2@-;mean gamma;');
idx = X(:,3)==gamma & X(:,1)==4; plot(X(idx,2), X(idx, column), '3-;thompson;');
idx = X(:,3)==gamma & X(:,1)==5; plot(X(idx,2), X(idx, column), '3@-;thompson gamma;');
idx = X(:,3)==gamma & X(:,1)==6; plot(X(idx,2), X(idx, column), '4-;high prob;');
idx = X(:,3)==gamma & X(:,1)==7; plot(X(idx,2), X(idx, column), '4@-;high prob gamma;');
xlabel("number of expansions");
ylabel("return");
legend("location", "southeast");
legend("boxon");
grid on
print("test_bandit_uct.eps");
