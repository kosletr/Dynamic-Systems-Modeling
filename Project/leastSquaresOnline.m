%% Online Least Squares

function [params_a2,params_b2] = leastSquaresOnline(t,y,step,u,lambda,beta,Q,n,m)

A=cell(n,1);
B=cell(n,1);
C=cell(n,1);
D=cell(n,1);

j=1;
for i=n:(-1):1
    temp_y = zeros(1,i); % construct term (-s^i)
    temp_y(1)=-1;
    [A{j},B{j},C{j},D{j}]=tf2ss(temp_y, [1,lambda]); % convert to state equations
    j=j+1;
end


for i=(m+1):(-1):1
    temp_u = zeros(1,i); % construct term (s^i)
    temp_u(1)=1;
    [A{j},B{j},C{j},D{j}]=tf2ss(temp_u, [1,lambda]); % convert to state equations
    j=j+1;
end

z_init=zeros((n+m+1)*n,1); % regression vector initial values
theta_init=zeros(n+m+1,1); % theta vector initial values
init_P=reshape(Q^(-1),[(n+m+1)^2,1]); % P initial values

[t,x] = ode45(@(t,x) odefun(t,x,step,u,y,A,B,C,D,n,m,beta),t,[z_init;theta_init;init_P],odeset('RelTol',1e-6,'AbsTol',1e-3));

estim_a=cell(1,n);
estim_b=cell(1,m+1);

j=1;
for i=1:1:n
    estim_a{j}=x(:,i+n*(n+m+1))+lambda(i);
    j=j+1;
end

j=1;
for i=n+1:1:n+m+1
    estim_b{j}=x(:,i+n*(n+m+1));
    j=j+1;
end

estim_a=cell2mat(estim_a);
estim_b=cell2mat(estim_b);

%%%%%%%%%%%%%%%%%%%%
Legend1=cell(1,n);
Legend2=cell(1,m+1);
if(n~=0)
    figure('Name','a Parameters Estimation')
    plot(t,estim_a)
    for i=1:1:n
        Legend1{i}=strcat('a_', num2str(i));
    end
    legend(Legend1);
end
figure('Name','b Parameters Estimation')
plot(t,estim_b)
for i=1:1:m+1
    Legend2{i}=strcat('b_', num2str(i)-1);
end
legend(Legend2);
%%%%%%%%%%%%%%%%%%%

params_a2=zeros(1,n);
params_b2=zeros(1,m+1);

for i=1:1:n
    params_a2(i)=estim_a(end,i);
end

for i=1:1:m+1
    params_b2(i)=estim_b(end,i);
    
end

end


function res = odefun(t,x,step,u,y,A,B,C,D,n,m,beta)
z_temp=cell(1,n*(n+m+1));
zeta=cell(1,n+m+1);
dot_z=cell(1,n+m+1);
theta=cell(1,n+m+1);

if(n==0) % case y=b_0*u
    phi=u(floor(t/step+1));
else
    k=1; % construct regression vector from state equations
    for i=1:1:n+m+1 % there are n+(m+1) zeta-components
        for j=1:1:n % each zeta component is produced out of n state equations
            z_temp{j}=x(k);
            k=k+1;
        end
        z=cell2mat(z_temp)';
        if i<=n % the first n components have y as input
            r=y;
        else
            r=u; % the last m+1 components have u as input
        end
        dot_z{i}=A{i}*z+B{i}*r(floor(t/step+1)); % DE to be solved by ODE - n(n+m+1) rows
        zeta{i}=C{i}*z+D{i}*r(floor(t/step+1))'; % Producing regression vector
    end
    dot_z=cell2mat(dot_z');
    phi=cell2mat(zeta');
end

for j=1:1:(n+m+1) % theta vector - n+(m+1) rows-components
    theta{j}=x(j+((n+m+1)*n));
end
theta=cell2mat(theta');

R=1e+19;
P = reshape(x((n+m+1)*(n+1)+1:end), [(n+m+1), (n+m+1)]); %use P as square matrix
if(norm(P)<=R)
    dP = beta*P-P*(phi*phi')*P;
else
    dP=zeros(n+m+1);
end

dP = reshape(dP, [(n+m+1)^2,1]); % convert dP from square matrix to column for ODE

if(n==0)
    res=[
        P*phi*(y(floor(t/step+1))-phi'*theta); % dot hat theta
        dP % dot P
        ];
else
    res=[
        dot_z; % Construct Regression Vectors
        P*phi*(y(floor(t/step+1))-phi'*theta); % dot hat theta
        dP % dot P
        ];
end
end