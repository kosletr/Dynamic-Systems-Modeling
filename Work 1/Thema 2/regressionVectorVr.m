function zVr = regressionVectorVr(t,x,u1,u2,lambda)

z1R=lsim([-1,0],[1,lambda(1),lambda(2)],x,t);
%figure('Name','z1R')
%plot(t,z1R)

z2R=lsim([-1],[1,lambda(1),lambda(2)],x,t);
%figure('Name','z2R')
%plot(t,z2R)

z3R=lsim([1],[1,lambda(1),lambda(2)],u1,t);
%figure('Name','z3C')
%plot(t,z3C)

z4R=lsim([1],[1,lambda(1),lambda(2)],u2,t);
%figure('Name','z4R')
%plot(t,z4R)

z5R=lsim([1,0,0],[1,lambda(1),lambda(2)],u1,t);
%figure('Name','z5R')
%plot(t,z5R)

z6R=lsim([1,0,0],[1,lambda(1),lambda(2)],u2,t);
%figure('Name','z6R')
%plot(t,z6R)

zVr=[z1R,z2R,z3R,z4R,z5R,z6R].';

end