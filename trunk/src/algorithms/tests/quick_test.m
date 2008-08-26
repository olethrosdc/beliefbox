## -*- Mode: octave -*-

X=load("quick_test.out");

grid on;
grid("minor");

subplot(2,2,1);
iter=2;
hold off;
s=X(:,1)==0 & X(:,2)==iter & X(:,4)==2; loglog(X(s,9), X(s,7)-X(s,5),'0')
hold on
s=X(:,1)==1 & X(:,2)==iter & X(:,4)==2; loglog(X(s,9), X(s,7)-X(s,5),'1')
s=X(:,1)==2 & X(:,2)==iter & X(:,4)==2; loglog(X(s,9), X(s,7)-X(s,5),'2')
s=X(:,1)==3 & X(:,2)==iter & X(:,4)==2; loglog(X(s,9), X(s,7)-X(s,5),'3')
axis([100,100000,1,10000]); grid on;
title("n_{iter}=2");

subplot(2,2,2);
iter=3;
hold off;
s=X(:,1)==0 & X(:,2)==iter & X(:,4)==2; loglog(X(s,9), X(s,7)-X(s,5),'0')
hold on
s=X(:,1)==1 & X(:,2)==iter & X(:,4)==2; loglog(X(s,9), X(s,7)-X(s,5),'1')
s=X(:,1)==2 & X(:,2)==iter & X(:,4)==2; loglog(X(s,9), X(s,7)-X(s,5),'2')
s=X(:,1)==3 & X(:,2)==iter & X(:,4)==2; loglog(X(s,9), X(s,7)-X(s,5),'3')
axis([100,100000,1,10000]); grid on; #rid("minor");
title("n_{iter}=3");

subplot(2,2,3);
iter=4;
hold off;
s=X(:,1)==0 & X(:,2)==iter & X(:,4)==2; loglog(X(s,9), X(s,7)-X(s,5),'0')
hold on
s=X(:,1)==1 & X(:,2)==iter & X(:,4)==2; loglog(X(s,9), X(s,7)-X(s,5),'1')
s=X(:,1)==2 & X(:,2)==iter & X(:,4)==2; loglog(X(s,9), X(s,7)-X(s,5),'2')
s=X(:,1)==3 & X(:,2)==iter & X(:,4)==2; loglog(X(s,9), X(s,7)-X(s,5),'3')
axis([100,100000,1,10000]); grid on; #grid("minor");
title("n_{iter}=4");

subplot(2,2,4);
iter=5;
hold off;
s=X(:,1)==0 & X(:,2)==iter & X(:,4)==2; loglog(X(s,9), X(s,7)-X(s,5),'0')
hold on
s=X(:,1)==1 & X(:,2)==iter & X(:,4)==2; loglog(X(s,9), X(s,7)-X(s,5),'1')
s=X(:,1)==2 & X(:,2)==iter & X(:,4)==2; loglog(X(s,9), X(s,7)-X(s,5),'2')
s=X(:,1)==3 & X(:,2)==iter & X(:,4)==2; loglog(X(s,9), X(s,7)-X(s,5),'3')
axis([100,100000,1,10000]); grid on; #grid("minor");
title("n_{iter}=5");

legend("serial", "random", "LB1", "LB2")
