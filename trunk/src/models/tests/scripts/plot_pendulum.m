## -*- Mode: octave -*-

tx=linspace(-1.5,1.5,100);
ty=linspace(-6,6,100);
[xx, yy] = meshgrid (tx, ty);
load Vout;
load Q0out;
load Q1out;
load Q2out

NC=64;

figure(1);
contour(xx,yy, Vout, NC);
xlabel("angle");
ylabel("angular velocity");
title("V")
print("V_pendulum.eps");

figure(2);
contour(xx,yy, Q0out, NC);
xlabel("angle");
ylabel("angular velocity");
title("Q(0)")
print("Q0_pendulum.eps");

figure(3);
contour(xx,yy, Q1out, NC);
xlabel("angle");
ylabel("angular velocity");
title("Q(1)")
print("Q1_pendulum.eps");

figure(4);
contour(xx,yy, Q2out, NC);
xlabel("angle");
ylabel("angular velocity");
title("Q(2)")
print("Q2_pendulum.eps")

figure(5);
contour(xx,yy, Q0out-Q2out, NC);
xlabel("angle");
ylabel("angular velocity");
title("Q(0) - Q(2)")
print("DQ02_pendulum.eps")
