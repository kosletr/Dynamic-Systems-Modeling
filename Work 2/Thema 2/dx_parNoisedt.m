function res = dx_parNoisedt(t,x,u,eta)

res=[-2*x(1)+u(t);
    -x(3)*x(2)+x(4)*u(t);
    -((x(1)+eta(t))-x(2))*x(2);
    ((x(1)+eta(t))-x(2))*u(t)
    ];
end