function out = f(x)
%Volume force, right hand side of Laplace equation: - div(grad(u)) = f

global t
x0 = 1.5 + [cos(t);sin(t)];
d2 = (x(:,1)-x0(1)).^2+(x(:,2)-x0(2)).^2;
out = exp(-10*d2);
