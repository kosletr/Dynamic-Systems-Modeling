function res = dx_mixNoisedt(t,x,u,theta_m,eta)

res=[-2*x(1)+u(t);
    -x(3)*x(2)+x(4)*u(t)-theta_m*((x(1)+eta(t))-x(2));
    -((x(1)+eta(t))-x(2))*(x(1)+eta(t));
    ((x(1)+eta(t))-x(2))*u(t)
    ];
end