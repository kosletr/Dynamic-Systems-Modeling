step=0.000001;
finish=10;
lambda=[10 25000];

t=0:step:finish;  
format long;

u1=2*sin(t);
u2= t>=0;

time = 0;
for i=1:length(t)
temp=v(time);
Vc(i)=temp(1);
Vr(i)=temp(2);
time = time + step;
end
x=[Vr.',Vc.'];

figure('Name','Inputs - Outputs');

subplot(3,1,1);
plot(t,u1);
hold on;
plot(t,u2);
title('Inputs')
xlabel('Time t (sec)')
legend('u_1(t)','u_2(t)');

subplot(3,1,2);
plot(t,x(:,1));
title('Output V_R')
xlabel('Time t (sec)')
legend('V_R(t)');

subplot(3,1,3);
plot(t,x(:,2));
title('Output V_C')
xlabel('Time t (sec)')
legend('V_C(t)');

[theta0VR,theta0VC]=leastSquares(t,x,lambda,u1,u2);

param=theta0VC+[lambda.';0;0;0;0];
RC_inv=param(1)
LC_inv=param(2)

y=lsim([0,LC_inv],[1,RC_inv,LC_inv],ones(length(t),1),t);
vc=y(:,1);
vr=u1.'+u2.'-y(:,1);

Vr_estim=vr;
Vc_estim=vc;
errVr=abs(x(:,1)-Vr_estim);
errVc=abs(x(:,2)-Vc_estim);

figure('Name','error functions');
subplot(2,1,1);
plot(t,errVr);
xlabel('Time t (sec)')
legend('e_R')
subplot(2,1,2);
plot(t,errVc);
suptitle('Error Functions');
xlabel('Time t (sec)')
legend('e_C')

x_noise=x;
for i=1:1:100
x_noise(randi((finish/step)+1),1)=randi(200)-100;
x_noise(randi((finish/step)+1),2)=randi(200)-100;
end

figure('Name','Inputs - Outputs (Noise)');

subplot(3,1,1);
plot(t,u1);
hold on;
plot(t,u2);
title('Inputs');
xlabel('Time t (sec)');
legend('u_1(t)','u_2(t)');

subplot(3,1,2);
plot(t,x_noise(:,1));
title('Output V_R (Noise)');
xlabel('Time t (sec)');
legend('V_R(t)');

subplot(3,1,3);
plot(t,x_noise(:,2));
title('Output V_C (Noise)');
xlabel('Time t (sec)');
legend('V_C(t)');


[theta0VRfalse,theta0VCfalse]=leastSquares(t,x_noise,lambda,u1,u2);
param_false=theta0VCfalse+[lambda.';0;0;0;0];
false_RC_inv=param_false(1)
false_LC_inv=param_false(2)

y=lsim([0,false_LC_inv],[1,false_RC_inv,false_LC_inv],ones(length(t),1),t);
vc=y(:,1);
vr=u1.'+u2.'-y(:,1);

Vr_estim=vr;
Vc_estim=vc;
errVr=abs(x(:,1)-Vr_estim);
errVc=abs(x(:,2)-Vc_estim);

figure('Name','error function (Noise)');
subplot(2,1,1);
plot(t,errVr);
xlabel('Time t (sec)')
legend('e_R (Noise)')
subplot(2,1,2);
plot(t,errVc);
suptitle('Error Functions (Noise)');
xlabel('Time t (sec)')
legend('e_C (Noise)')