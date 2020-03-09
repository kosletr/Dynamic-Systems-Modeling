function res = dxdt(t,x,u,theta_m,gamma)

res=[
    -2*x(1)+u(t);
    -x(2)*theta_m+x(1);
    -x(3)*theta_m+u(t);
    gamma*(x(2)*x(1)-x(2)*x(2)*x(4)-x(2)*x(3)*x(5));
    gamma*(x(3)*x(1)-x(3)*x(2)*x(4)-x(3)*x(3)*x(5))
    ];
end