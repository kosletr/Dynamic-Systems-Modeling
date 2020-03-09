function res = dx_pardt(t,x,u)

res=[-2*x(1)+u(t);
    -x(3)*x(2)+x(4)*u(t);
    -(x(1)-x(2))*x(2);
    (x(1)-x(2))*u(t)
    ];
end