step=0.0001;
finish=30;
t=0:step:finish;
format long;
theta_m=0.01;

u= @(a) 5*sin(3*a);


[t,y] = ode45(@(t,y) dx_pardt(t,y,u),t,[0;0;0;0]);
y_real=y(:,1);
y_estim=y(:,2);
theta_1_par=y(:,3);
theta_2_par=y(:,4);

figure('Name','x and x Estimation Parr');
subplot(2,1,1)
plot(t,y_real);
title('Output x');
xlabel('Time t (sec)')
subplot(2,1,2)
plot(t,y_estim);
title('Output x Estimation (Parallel)');
xlabel('Time t (sec)')

figure('Name','a Estimation Parr');
plot(t,theta_1_par);
title('Parameter a Estimation (Parralel)');
xlabel('Time t (sec)')

figure('Name','b Estimation Parr');
plot(t,theta_2_par);
title('Parameter b Estimation Parr');
xlabel('Time t (sec)')


[t,x] = ode45(@(t,x) dx_mixdt(t,x,u,theta_m),t,[0;0;0;0]);
x_real=x(:,1);
x_estim=x(:,2);
theta_1_mix=x(:,3);
theta_2_mix=x(:,4);

figure('Name','x and x Estimation Mix');
subplot(2,1,1)
plot(t,x_real);
title('Output x');
xlabel('Time t (sec)')
subplot(2,1,2)
plot(t,x_estim);
title('Output x Estimation Mix');
xlabel('Time t (sec)')

figure('Name','a Estimation Mix');
plot(t,theta_1_mix);
title('Parameter a Estimation Mix');
xlabel('Time t (sec)')

figure('Name','b Estimation Mix');
plot(t,theta_2_mix);
title('Parameter b Estimation Mix');
xlabel('Time t (sec)')

%%%%%%%%%%%%%%%%%%%%Noise%%%%%%%%%%%%%%%%%%%%%%%
eta = @(a) 0.15*sin(40*pi*a);

[t,yn] = ode45(@(t,yn) dx_parNoisedt(t,yn,u,eta),t,[0;0;0;0]);
yn_real=yn(:,1);
yn_estim=yn(:,2);
theta_1n_par=yn(:,3);
theta_2n_par=yn(:,4);

figure('Name','x Noise and x Estimation (Parallel)');
subplot(2,1,1);
plot(t,yn_real);
title('Output x');
xlabel('Time t (sec)')
subplot(2,1,2);
plot(t,yn_estim);
title('Output x Estimation (Parallel) - Noise');
xlabel('Time t (sec)')

figure('Name','a Estimation Parr - Noise');
plot(t,theta_1n_par);
title('Parameter a Estimation Parr - Noise');
xlabel('Time t (sec)')

figure('Name','b Estimation Parr - Noise');
plot(t,theta_2n_par);
title('Parameter b Estimation Parr - Noise');
xlabel('Time t (sec)')


[t,xn] = ode45(@(t,xn) dx_mixNoisedt(t,xn,u,theta_m,eta),t,[0;0;0;0]);
xn_real=xn(:,1);
xn_estim=xn(:,2);
theta_1n_mix=xn(:,3);
theta_2n_mix=xn(:,4);

figure('Name','Parameter x Estimation Mix - Noise');
subplot(2,1,1);
plot(t,xn_real);
title('Output x - Noise');
xlabel('Time t (sec)')
subplot(2,1,2);
plot(t,xn_estim);
title('Output x Estimation (Mix) - Noise');
xlabel('Time t (sec)')

figure('Name','a Estimation Mix - Noise');
plot(t,theta_1n_mix);
title('Parameter a Estimation Mix - Noise');
xlabel('Time t (sec)')

figure('Name','b Estimation Mix - Noise');
plot(t,theta_2n_mix);
title('Parameter b Estimation Mix - Noise');
xlabel('Time t (sec)')