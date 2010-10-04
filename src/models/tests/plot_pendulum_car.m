linewidth=3;
fontsize=24;
FF = "-F:24";
param_string="linewidth";
params=linewidth;
set (gca(), "defaultlinelinewidth", linewidth)
set (gca(), "defaulttextfontsize", fontsize)

K = 1000;

x1=moving_average(load("MountainCar_d1.out"), K);
x2=moving_average(load("MountainCar_d2.out"), K);
x3=moving_average(load("MountainCar_d3.out"), K);
x4=moving_average(load("MountainCar_d4.out"), K);
x5=moving_average(load("MountainCar_d5.out"), K);
x6=moving_average(load("MountainCar_d6.out"), K);
x7=moving_average(load("MountainCar_d7.out"), K);
x8=moving_average(load("MountainCar_d8.out"), K);
if (0)
z1=moving_average(load("MountainCar3D_d1.out"), K);
z2=moving_average(load("MountainCar3D_d2.out"), K);
z3=moving_average(load("MountainCar3D_d3.out"), K);
z4=moving_average(load("MountainCar3D_d4.out"), K);
z5=moving_average(load("MountainCar3D_d5.out"), K);
z6=moving_average(load("MountainCar3D_d6.out"), K);
z7=moving_average(load("MountainCar3D_d7.out"), K);
z8=moving_average(load("MountainCar3D_d8.out"), K);
endif
y1=moving_average(load("Pendulum_d1.out"), K);
y2=moving_average(load("Pendulum_d2.out"), K);
y3=moving_average(load("Pendulum_d3.out"), K);
y4=moving_average(load("Pendulum_d4.out"), K);
y5=moving_average(load("Pendulum_d5.out"), K);
y6=moving_average(load("Pendulum_d6.out"), K);
y7=moving_average(load("Pendulum_d7.out"), K);
y8=moving_average(load("Pendulum_d8.out"), K);

figure(1)
plot(cumsum(1 + x1(:,2)), ';1;',
     cumsum(1 + x2(:,2)), ';2;',
     cumsum(1 + x3(:,2)), ';3;',
     cumsum(1 + x4(:,2)), ';4;',
     cumsum(1 + x5(:,2)), ';5;',
     cumsum(1 + x6(:,2)), ';6;',
     cumsum(1 + x7(:,2)), ';7;',
     cumsum(1 + x8(:,2)), '0;8;');
title("mountain car");
if (0)
figure(2)
plot(cumsum(1 + z1(:,2)), ';1;',
     cumsum(1 + z2(:,2)), ';2;',
     cumsum(1 + z3(:,2)), ';3;',
     cumsum(1 + z4(:,2)), ';4;',
     cumsum(1 + z5(:,2)), ';5;',
     cumsum(1 + z6(:,2)), ';6;',
     cumsum(1 + z7(:,2)), ';7;',
     cumsum(1 + z8(:,2)), '0;8;');
title("3D mountain car");
endif
figure(3);
plot(y1(:,2), ';1;',
     y2(:,2), ';2;',
     y3(:,2), ';3;',
     y4(:,2), ';4;',
     y5(:,2), ';5;',
     y6(:,2), ';6;',
     y7(:,2), ';7;',
     y8(:,2), '0;8;');
title("Pendulum");