x = load out;
x0 = load out0;
x1 = load out1;
x2 = load out2;
x4 = load out4;
x8 = load out8;

K=1000;

plot(moving_average(x, K), ';Q;',
     moving_average(x0, K), ';0;',
     moving_average(x1, K), ';1;',
     moving_average(x2, K), ';2;',
     moving_average(x4, K), ';4;',
     moving_average(x8, K), ';8;')

plot(moving_average(x, K), ';Q;',
     moving_average(x0, K), ';0;',
     moving_average(x1, K), ';1;',
     moving_average(x2, K), ';2;');
