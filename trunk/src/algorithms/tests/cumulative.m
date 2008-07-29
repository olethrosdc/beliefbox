## -*- Mode: octave -*-
## test methods for finding
## P(V(x) > v)
## when we know the density p(x) but calculating the integral
## \int I(V(x) > v) p(x) dx,
## is very difficult.

1;

function v = Value(x)
  v = sinc(x/2 - 2);
end

function p = Density(x)
  p = normal_pdf(x, 0, 1);
end

function X = Sample(N)
  X = normal_rnd(0, 1, N, 1);
end

x = [-100:100]/10;
plot(x, Value(x), ";value;", x, Density(x), ";density;", x, Value(x).*Density(x), ";ED;");

N=1000;
X = Sample(N);
P = Density(X);
EV = X.*Value(X);

