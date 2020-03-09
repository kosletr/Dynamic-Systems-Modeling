%% Check if System is Linear Time Invariant

function LTICheck(t,step,finish)
%% Check Linearity

c_1=rand();
c2=rand();
u1 = sin(t)+2*cos(3*t)+5*cos(7*t)+9*sin(10.5*t);
y1=out(t,u1);
u2 = 12*cos(3.5*t)+15*sin(2.5*t) + cos(2*t);
y2=out(t,u2);
y_sum=c_1*y1+c2*y2;
u_sum=c_1*u1+c2*u2;
y3=out(t,u_sum);
error_lin=y3-y_sum;

figure('Name','Check Input Superposition');
plot(t,error_lin);

%% Check Time Invariance

t0=randperm(floor(finish/2),1);

t_d=t0:step:finish+t0;
u_func = @(a) cos(a)+3*sin(5*a)+8*cos(8*a)+3*sin(15*a);


y1=out(t,u_func(t_d)); %y(u(t+t_0))
warning off all % disable warning due to Simulation start at a nonzero initial time
y2=out(t_d,u_func(t_d)); %y(t+t_0)
warning on all % reenable warnings

error_time_inv=y1-y2;
figure('Name','Check Time Invariance');
plot(t,error_time_inv);

end