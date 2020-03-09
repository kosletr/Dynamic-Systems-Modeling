function [theta0Vr,theta0Vc] = leastSquares(t,x,lambda,u1,u2)

Vr=x(:,1);
Vc=x(:,2);

phiVr=regressionVectorVr(t,Vr,u1,u2,lambda);
phiVc=regressionVectorVc(t,Vc,u1,u2,lambda);

phi_phi_Vr=phiVr*phiVr.';
phi_phi_Vc=phiVc*phiVc.';

phi_Vr=phiVr*Vr;
phi_Vc=phiVc*Vc;

theta0Vr = (phi_phi_Vr)^(-1)*phi_Vr;
theta0Vc = (phi_phi_Vc)^(-1)*phi_Vc;
end