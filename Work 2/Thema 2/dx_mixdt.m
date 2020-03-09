function res = dx_mixdt(t,x,u,theta_m)

res=[-2*x(1)+u(t);
    -x(3)*x(2)+x(4)*u(t)-theta_m*(x(1)-x(2));
    -(x(1)-x(2))*x(1);
    (x(1)-x(2))*u(t)
    ];
end