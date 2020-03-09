step=0.0001;
finish=60;
t=0:step:finish;
format long;

u = @(a) 10*sin(2*a)+5*sin(7.5*a);

[t,x] = ode45(@(t,x) dx_dt(t,x,u),t,[0;0;0;0;0;0;0;0;0;0]);
x1_real=x(:,1);
x2_real=x(:,2);
x1_estim=x(:,3);
x2_estim=x(:,4);
a11_estim=x(:,5);
a12_estim=x(:,6);
a21_estim=x(:,7);
a22_estim=x(:,8);
b1_estim=x(:,9);
b2_estim=x(:,10);


figure('Name','x1 and x1 Estimation');
subplot(2,1,1);
plot(t,x1_real);
title('x_1 Real');
xlabel('Time t (sec)')
subplot(2,1,2);
plot(t,x1_estim);
title('x_1 Estimation');
xlabel('Time t (sec)')

figure('Name','x2 and x2 Estimation');
subplot(2,1,1);
plot(t,x2_real);
title('x_2 Real');
xlabel('Time t (sec)')
subplot(2,1,2);
plot(t,x2_estim);
title('x_2 Estimation');
xlabel('Time t (sec)')


figure('Name','a11 Estimation');
plot(t,a11_estim);
title('Parameter a11 Estimation');
xlabel('Time t (sec)')

figure('Name','a12 Estimation');
plot(t,a12_estim);
title('Parameter a12 Estimation');
xlabel('Time t (sec)')

figure('Name','a21 Estimation');
plot(t,a21_estim);
title('Parameter a21 Estimation');
xlabel('Time t (sec)')

figure('Name','a22 Estimation');
plot(t,a22_estim);
title('Parameter a22 Estimation');
xlabel('Time t (sec)')

figure('Name','b1 Estimation');
plot(t,b1_estim);
title('Parameter b1 Estimation');
xlabel('Time t (sec)')

figure('Name','b2 Estimation');
plot(t,b2_estim);
title('Parameter b2 Estimation');
xlabel('Time t (sec)')