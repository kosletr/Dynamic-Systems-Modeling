%% Offline Least Squares

function [params_a,params_b] = leastSquaresOffline(t,y,u,lambda,n,m)

phi=regressionVectorOffline(t,y,u,lambda,n,m)';

phi_phi=phi*phi.';
phi_y=phi*y;
theta0 = linsolve(phi_phi,phi_y);

params_a=zeros(1,n);
params_b=zeros(1,m+1);

for i=1:1:n
    params_a(i)=theta0(i,1)+lambda(i); % parameters_a = theta_1 + lambda
end
for i=1:1:(m+1)
    params_b(i)=theta0(i+n,1);  % parameters_b = theta_2
end
end

function z = regressionVectorOffline(t,y,u,lambda,n,m)
z_y=cell(1,n);
z_u=cell(1,m+1);

j=1;
for i=n:(-1):1
    temp = zeros(1,i); % construct term (-s^i)
    temp(1)=-1;
    z_y{j}=lsim(temp,[1,lambda],y,t); % solve (-y*s^i)/(Lamda(s))
    j=j+1;
end

j=1;
for i=(m+1):(-1):1
    temp = zeros(1,i); % construct term (s^i)
    temp(1)=1;
    z_u{j}=lsim(temp,[1,lambda],u,t); %solve (u*s^i)/(Lamda(s))
    j=j+1;
end

z=[cell2mat(z_y),cell2mat(z_u)]; % construct regression vector
end