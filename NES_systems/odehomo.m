function yprime=odehomo(t,y)
global A ep f1 k_n w M
F=[0;0;M\[ep*f1*cos(w*t)-ep*k_n*(y(1)-y(2))^3;-ep*k_n*(y(2)-y(1))^3]];
yprime=A*y+F;