## -*- Mode: octave -*-

fontsize=32
linewidth=5;
FF = "-F:32";

set (gcf(), "defaultlinelinewidth", linewidth)
set(gcf(), "defaulttextfontsize", fontsize)


s="100000_16_16"; 
f_xy = ["P_XY_", s]; 
f_y_x = ["P_Y_X_", s];
P_XY=load (f_xy); 
P_Y_X=load (f_y_x); 

figure(1); 
imagesc(P_XY);
xlabel("x");
ylabel("y"); 
title("p(x, y)"); 
print([f_xy, ".eps"], FF, "-color"); 

figure(2); 
imagesc(P_Y_X); 
xlabel("x");
ylabel("y"); 
title("p(y | x)"); 
print([f_y_x, ".eps"], FF, "-color")
