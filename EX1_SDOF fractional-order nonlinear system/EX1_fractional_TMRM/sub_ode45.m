function yprime=sub_ode45(t,y)
global M A beta r 
F=[0;0;-M\[0.008*cos(0.2*t)*y(1);0.135*y(2)^2*y(4)]];
yprime=A*y+F;
