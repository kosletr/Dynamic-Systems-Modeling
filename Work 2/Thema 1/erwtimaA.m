step=0.0001;
finish=30;
t=0:step:finish;
format long;
theta_m=5;
gamma=30;

u= @(a) 5*sin(3*a);

[t,x] = ode45(@(t,x) dxdt(t,x,u,theta_m,gamma),t,[0;0;0;theta_m;0]);
a_estim=theta_m-x(:,4);
b_estiim=x(:,5);

x_estim=x(:,4).*x(:,2)+x(:,5).*x(:,3);

figure('Name','Input-Output');
plot(t,u(t));
hold on;
plot(t,x(:,1));
title('Input - Output');
legend('u(t)','x(t)');
xlabel('Time t (sec)')


figure('Name','a Estimation');
plot(t,a_estim);
title('Parameter a Estimation');
xlabel('Time t (sec)')

figure('Name','b Estimation');
plot(t,b_estiim);
title('Parameter b Estimation');
xlabel('Time t (sec)')

figure('Name','x Estimation');
subplot(2,1,1);
plot(t,x_estim);
title('Output x Estimation');
xlabel('Time t (sec)')
subplot(2,1,2);
plot(t,x(:,1));
title('Output x Estimation');
xlabel('Time t (sec)')