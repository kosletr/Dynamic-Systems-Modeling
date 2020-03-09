function [m,b,k] = leastSquares(t,y,u,lambda)

phi=regressionVector(t,y,u,lambda);
phi_phi=phi*phi.';

phi_y=phi*y;

theta0 = linsolve(phi_phi,phi_y)
m=1/theta0(3);
k=m*(lambda(2)+theta0(2));
b=m*(lambda(1)+theta0(1));
end