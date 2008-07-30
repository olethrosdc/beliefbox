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

x = [-100000:100000]/10000;
Vx = Value(x);
Px = Density(x);
SPx = sum(Px);

#EVx = Vx .* Px;
figure(1);
plot(x, Vx, ";value;", x, Px, ";density;", [-10, 10], Vx*Px'*[1 1]/SPx, ";EV;");

N=10000;
X = Sample(N);
P = Density(X);
V = Value(X);
#EV = X.*Value(X);

vi = [-25:100]/100;
vn = length(vi);

cdf = zeros(vn,1);

SP = sum(P);
for i=1:vn
  v = vi(i);
  cdf(i) = sum(Px(Vx > v)) / SPx;
  cdf_b(i) = sum(V > v) / length(V);
  cdf_bw(i) = sum(P(V > v)) / SP;
endfor

figure(2);
semilogy(vi, cdf, ";cdf;", vi, cdf_b, ";MC;", vi, cdf_bw, ";MC IS;")