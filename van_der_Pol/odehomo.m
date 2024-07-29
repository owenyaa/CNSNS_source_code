function yprime=odehomo(t,y)
global A e a f1 f2 w1 w2 
F=[0;e*y(1)^2*y(2)-a*y(1)^3+f1*cos(w1*t)+f2*cos(w2*t)];
yprime=A*y+F;