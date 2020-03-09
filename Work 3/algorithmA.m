step=0.001;
finish=10;
format long;

theta_m=2;
gamma_1=20;
gamma_2=30;

u = @(a) 5*sin(3*a);
init_values=[[theta_m-1.999;2],[theta_m-0.001;2],[theta_m-1.999;1.1],[theta_m-0.001;1.1],...
    [theta_m-0.5;1.2],[theta_m-0.75;1.6],[theta_m-1.8;1.7],[theta_m-1.3;1.2]];
fig = @(a) 3*floor(a/4.001);

for i=1:1:size(init_values,2)
z{i}=algoA(init_values(:,i),u,theta_m,gamma_1,gamma_2,step,finish);
end


for i=1:1:size(init_values,2)

t=0:step:(size(z{i}(:,1))-1)*step;
    
a_estim=z{i}(:,4)-theta_m;
b_estim=z{i}(:,5);
figure(1+fig(i));
subplot(2,2,mod(i+3,4)+1)
plotRect('r');
hold on;
plot(a_estim,b_estim,'Color','b')
title('a_{est} - b_{est} trajectory');
xlabel('a estimation')
ylabel('b estimation')
axis([-2.5 0.5 0.5 2.5]);


figure(2+fig(i));
subplot(2,2,mod(i+3,4)+1)
plot(t,a_estim);
hold on;
plot(t,-2*ones(size(t)),'-.','color','r')
hold on;
plot(t,zeros(size(t)),'-.','color','r')
title('Parameter a Estimation');
xlabel('Time t (sec)')
axis(gca,'tight')
ylim([-2.2 0.2])


figure(3+fig(i));
subplot(2,2,mod(i+3,4)+1)
plot(t,b_estim);
hold on;
plot(t,1.1*ones(size(t)),'-.','color','r')
hold on;
plot(t,2*ones(size(t)),'-.','color','r')
title('Parameter b Estimation');
xlabel('Time t (sec)')
axis(gca,'tight')
ylim([0.9 2.2])

end

figure('Name','Input-Output');
plot(t,u(t));
hold on;
plot(t,z{end}(:,1));
title('Input - Output');
legend('u(t)','x(t)');
xlabel('Time t (sec)')
axis(gca,'tight')


x_estim=z{end}(:,4).*z{end}(:,2)+z{end}(:,5).*z{end}(:,3);

figure('Name','Output x Error');
plot(t,abs(z{i}(:,1)-x_estim));
title('Output x Error');
xlabel('Time t (sec)')
axis(gca,'tight')

function x=algoA(init_values,u,theta_m,gamma_1,gamma_2,step,finish)
t=0:step:finish;

options = odeset('Events', @(t, x) border_reached(t, x, theta_m,gamma_1,gamma_2));
init_val=[0;0;0;init_values];

time = t;
x=[];
te=-1;
while (isempty(te)==false)
[time,y,te,xe,ie] = ode45(@(time,y) odefun(time,y,u,theta_m,gamma_1,gamma_2),time,init_val,options);
init_val=[y(end,1);y(end,2);y(end,3);y(end,4);y(end,5)];
x=[x;y(:,1),y(:,2),y(:,3),y(:,4),y(:,5)];
time=te:step:finish;
end
end

%----------------Event function --------------------
function [value,isterminal,direction] = border_reached(time,x,theta_m,gamma_1,gamma_2)

g1 =@(a) -a+theta_m-2;
g2 =@(a) a-theta_m;
g3 =@(a) -a+1.1;
g4 =@(a) a-2;
Gamma=[gamma_1,0;0,gamma_2];
negGGradK=Gamma*[x(2);x(3)]*(x(1)-[x(2),x(3)]*[x(4);x(5)]);

value=false;

if(( g1(x(4))>=0 && (negGGradK.' * [-1;0])>=0 )) 
     value=true;
 elseif(( g2(x(4))>=0 && (negGGradK.' * [1;0])>=0 ) )
     value=true;
 elseif(( g3(x(5))>=0 && (negGGradK.' * [0;-1])>=0 ))
     value=true;
 elseif(( g4(x(5))>=0 && (negGGradK.' * [0;1])>=0 ) )
    value=true;   
end

isterminal = 1; 
direction = 1;
end

function res = odefun(t,x,u,theta_m,gamma_1,gamma_2)

g1 =@(a) -a+theta_m-2;
g2 =@(a) a-theta_m;
g3 =@(a) -a+1.1;
g4 =@(a) a-2;

Gamma=[gamma_1,0;0,gamma_2];
negGGradK=Gamma*[x(2);x(3)]*(x(1)-[x(2),x(3)]*[x(4);x(5)]);
 
x45=negGGradK;

if(( g1(x(4))>=0 && (negGGradK.' * [-1;0])>=0 )) 
     gradG=[-1;0];
     x45=negGGradK-Gamma*(gradG*(gradG.')*negGGradK)/((gradG.')*Gamma*gradG);
 elseif(( g2(x(4))>=0 && (negGGradK.' * [1;0])>=0 ) )
     gradG=[1;0];
     x45=negGGradK-Gamma*(gradG*(gradG.')*negGGradK)/((gradG.')*Gamma*gradG);
 elseif(( g3(x(5))>=0 && (negGGradK.' * [0;-1])>=0 ))
     gradG=[0;-1];
     x45=negGGradK-Gamma*(gradG*(gradG.')*negGGradK)/((gradG.')*Gamma*gradG);
 elseif(( g4(x(5))>=0 && (negGGradK.' * [0;1])>=0 ) )
    gradG=[0;1];
    x45=negGGradK-Gamma*(gradG*(gradG.')*negGGradK)/((gradG.')*Gamma*gradG);    
end

    res=[
    -1.2*x(1)+1.55*u(t); %x
    -x(2)*theta_m+x(1); %ph1
    -x(3)*theta_m+u(t); %phi2
    x45 %theta1 - theta2
    ]; 
     
 end


function plotRect(color)
A = [-2,-2];
B = [1.1,2];
C = [0,0];
D = [-2,0];
E = [2,2];
F=[1.1,1.1];
plot(A,B,'--','Color',color);
hold on;
plot(C,B,'--','Color',color);
hold on;
plot(D,E,'-','Color',color);
hold on;
plot(D,F,'-','Color',color);
end