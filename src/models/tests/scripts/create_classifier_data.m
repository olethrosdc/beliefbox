## -*- Mode: octave -*-


c = randn(2,2);

p_class = 0.3;

n=1000;

printf("%d %d\n", n, 3);
for i=1:n
  if (rand < p_class) 
	class = 1;
  else
	class = 2;
  endif
  x = randn(2,1) + c(2, class);
  printf ("%f %f %d\n", x(1), x(2), class - 1);
endfor
