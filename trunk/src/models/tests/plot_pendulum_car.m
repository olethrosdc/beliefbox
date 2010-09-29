K = 1000;

x1=moving_average(load("MountainCar_d1.out"), K);
x2=moving_average(load("MountainCar_d2.out"), K);
x3=moving_average(load("MountainCar_d3.out"), K);
x4=moving_average(load("MountainCar_d4.out"), K);
x5=moving_average(load("MountainCar_d5.out"), K);

y1=moving_average(load("Pendulum_d1.out"), K);
y2=moving_average(load("Pendulum_d2.out"), K);
y3=moving_average(load("Pendulum_d3.out"), K);
y4=moving_average(load("Pendulum_d4.out"), K);
y5=moving_average(load("Pendulum_d5.out"), K);

subplot(2,1,1);
plot(cumsum(0.01 + x1(:,2)), ';1;',
     cumsum(0.01 + x2(:,2)), ';2;',
     cumsum(0.01 + x3(:,2)), ';3;',
     cumsum(0.01 + x4(:,2)), ';4;',
     cumsum(0.01 + x5(:,2)), ';5;')

subplot(2,1,2);
plot(y1(:,2), ';1;',
     y2(:,2), ';2;',
     y3(:,2), ';3;',
     y4(:,2), ';4;',
     y5(:,2), ';5;')
