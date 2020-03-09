function zVc = regressionVectorVc(t,x,u1,u2,lambda)

z1C=lsim([-1,0],[1,lambda(1),lambda(2)],x,t);
%figure('Name','z1C')
%plot(t,z1C)

z2C=lsim([-1],[1,lambda(1),lambda(2)],x,t);
%figure('Name','z2C')
%plot(t,z2C)

z3C=lsim([1],[1,lambda(1),lambda(2)],u1,t);
%figure('Name','z3C')
%plot(t,z3C)

z4C=lsim([1],[1,lambda(1),lambda(2)],u2,t);
%figure('Name','z4C')
%plot(t,z4C)

z5C=lsim([1,0],[1,lambda(1),lambda(2)],u1,t);
%figure('Name','z5C')
%plot(t,z5C)

z6C=lsim([1,0],[1,lambda(1),lambda(2)],u2,t);
%figure('Name','z6C')
%plot(t,z6C)

zVc=[z1C,z2C,z3C,z4C,z5C,z6C].';

end