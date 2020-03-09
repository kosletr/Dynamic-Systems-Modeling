k=2;
b=0.2;
m=15;

t=0:0.1:10;
format long;

u_func = @(x) 5*sin(2*x)+10.5;
u=feval(u_func,t).';

y=lsim([1/m],[1,b/m,k/m],u,t);
figure('Name','y diff solution + u input');
plot(t,u);
hold on;
plot(t,y)
title('Input - Output')
xlabel('Time t (sec)')
legend('Input','Output')

lambda=[0.13 0.013]; % polynomial coefficients
[m_est,b_est,k_est] = leastSquares(t,y,u,lambda)

y_est=lsim([1/m_est],[1,b_est/m_est,k_est/m_est],u,t);
err=abs(y-y_est);
figure('Name','error function');
plot(t,err);
title('Error');
xlabel('Time t (sec)');