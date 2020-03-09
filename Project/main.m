step=0.0001;
finish=100;
t=0:step:finish;
format long;

%% Check Linearity - Time Invariance

LTICheck(t,step,finish);

%% LTI Model Selection

%Manual Selection

n=3; % s^n Y +a1*s^(n-1) Y +...+an*Y
m=2; % b0*s^m U +b1*s^(m-1) U +...+bm*U

%Automatic Selection

%Uncomment these lines to test many models 
% for i=0:4
%    for j=i:4
%        n=j
%        m=i

% Check System's Validation

if(n<m || n<0 || m<0)
    disp("Invalid Model Selection (n,m)")
%   continue; % Uncoment this line as well to test more models
    return;
end
%% METHODS - Least Squares
%% Filter Selection

%Automatic Selection

syms s
filter=(s+1)^n; %Lambda(s)=(s+1)^n
lambda = double(fliplr(coeffs(filter))); 
lambda(1)=[];
clear s;

%Manual Selection

%Uncomment these lines to test a particular filter
%lambda=[3 3 1];
%if(length(lambda)~=n)
%    disp("Invalid Filter Selection")
%    return;
%end

%% Offline Least Squares Method

u=cos(t')+7*sin(2*t')+3*cos(4*t');
y=out(t,u);

[params_a1,params_b1]=leastSquaresOffline(t,y,u,lambda,n,m)

%%% Offline Least Squares Method's Results
u=5*sin(10*t')+3*cos(3*t')+sin(1.5*t');
checkResults(t,u,params_a1,params_b1)

%% Online Least Squares Method

u=cos(t')+7*sin(2*t')+3*cos(4*t');
y=out(t,u);

Q=eye(n+m+1);
beta=0.2;

[params_a2,params_b2]=leastSquaresOnline(t,y,step,u,lambda,beta,Q,n,m)

%%% Online Least Squares Method's Results
u=5*sin(10*t')+3*cos(3*t')+sin(1.5*t');
checkResults(t,u,params_a2,params_b2)

% end % Uncoment this line as well to test more models
% end % Uncoment this line as well to test more models
%% Check Results

function checkResults(t,u,params_a,params_b)
y=out(t,u);
z=lsim(params_b,[1,params_a],u,t);
figure('Name','Least Squares Error');
plot(t,abs(y-z))
error=sum(abs(y-z))/length(t)
end
