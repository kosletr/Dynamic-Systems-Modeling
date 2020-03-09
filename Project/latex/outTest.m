function res=out(t,u)
a1=6;
a2=11;
a3=6;

b0=1;
b1=-3;
b2=2;

res=lsim([b0 b1 b2],[1 a1 a2 a3],u,t);
end
