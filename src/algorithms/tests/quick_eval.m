## -*- Mode: octave -*_
X=load("quick_eval2.out");
XU=load("quick_eval_ucb.out");
gamma = 0.999;

Y=X(X(:,2)==gamma, :);
YU=XU(XU(:,2)==gamma, :);

col=5;
for i=1:8;
  K=Y(Y(:,1)==i-1,:); Z(:,i)=K(:,6)-K(:,col);
endfor

for i=1:2
  KU=YU(YU(:,1)==i-1,:); Z(:,i+8)=KU(:,6)-KU(:,col);
endfor

plot(K(:,3), Z);
legend("Serial", "Random", "Mean", "Discounted Mean",  "Thompson", "Thompson gamma", "TBound", "Tbound gamma", "UCB1", "UCBgamma");
title("gamma = 0.99");
xlabel("depth");
ylabel("regret");
print("quick_eval999.eps");

