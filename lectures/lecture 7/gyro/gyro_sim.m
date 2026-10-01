close all
clear all

m = 0.001;
k = 10.0;
c = 0.1*sqrt(2*m*k);
a = 1.0;

r0 = [0; 0; 0];
rdot0 = [0; 0; 0];

x0 = [r0; rdot0];

x02 = zeros(4,1);

Tfinal = 5.0;

[t, xtraj] = ode45(@(t,x)mems_gyro(t,x,@omega,m,k,c,a),[0 Tfinal], x0);

[t2, xtraj2] = ode45(@(t,x)mems_gyro2(t,x,@omega,m,k,c,a),[0 Tfinal], x02);


subplot(3,1,1)
plot(t,xtraj(:,1))
hold on
plot(t2,xtraj2(:,1));
subplot(3,1,2)
plot(t,xtraj(:,2))
hold on
plot(t2,xtraj2(:,2))
subplot(3,1,3)
plot(t,omega(t))