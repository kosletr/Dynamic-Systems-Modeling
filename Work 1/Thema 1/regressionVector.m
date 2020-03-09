function z = regressionVector(t,y,u,lambda)

z1=lsim([-1,0],[1,lambda(1),lambda(2)],y,t);
%figure('Name','z1 diff solution');
%plot(t,z1);

z2=lsim([-1],[1,lambda(1),lambda(2)],y,t);
%figure('Name','z2 diff solution');
%plot(t,z2);

z3=lsim([1],[1,lambda(1),lambda(2)],u,t);
%figure('Name','z3 diff solution');
%plot(t,z3);

z=[z1,z2,z3].';
end