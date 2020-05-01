clc; clear; close all

%% 精确解
M=1; K=100; U=0;
[x1,x2] = dsolve('Dx1=1/M*x2,Dx2=-K*x1+U','x1(0)=0,x2(0)=1');
h = 0.001;
time = 0: h: 2;
for i=1:length(time)
    t=time(i);
    val_x1(i)=eval(x1);
    val_x2(i)=eval(x2);
end
figure
plot(time,val_x1 ,'y--',time,val_x2, 'y-.', 'linewidth', 3)
hold on

clear
%% 欧拉法
h = 0.001;
t = 0: h: 2;
n = length(t);
x(:, 1) = [0; 1];

for k = 1: n-1
    x(:, k + 1) = x(:, k) + h * f(t(k), x(:, k));
end
x1_Euler = x(1, :);
x2_Euler = x(2, :);
plot(t,x1_Euler,'r--',t,x2_Euler,'r-.')
hold on

clear
%% 改进欧拉法
h = 0.001;
t = 0: h: 2;
n = length(t);
x(:, 1) = [0; 1];

for k = 1: n-1
    K1 = f(t(k), x(:, k));
    K2 = f(t(k+1), x(:, k) + h * K1);
    x(:, k + 1) = x(:, k) + h * (K1 + K2)/2;
end
x1_Improve_Euler = x(1, :);
x2_Improve_Euler = x(2, :);
plot(t,x1_Improve_Euler ,'g--',t,x2_Improve_Euler,'g-.')
hold on

clear
%% 4阶 Runge-Kutta 法
h = 0.001;
t = 0: h: 2;
n = length(t);
x(:, 1) = [0; 1];

for k = 1: n-1
    K1 = f(t(k), x(:, k));
    K2 = f(t(k) + h/2, x(:, k) + h/2 * K1);
    K3 = f(t(k) + h/2, x(:, k) + h/2 * K2);
    K4 = f(t(k+1) , x(:, k) + h * K3);
    x(:, k + 1) = x(:, k) + h * (K1 + 2*K2 + 2*K3 +K4)/6;
end
x1_RK = x(1, :);
x2_RK = x(2, :);
plot(t, x1_RK, 'b--', t, x2_RK, 'b-.')

ylim([-1.5 1.5]);
grid minor
legend('x_1精确解','x_2精确解','x_1欧拉法','x_2欧拉法',...
    'x_1改进欧拉法','x_2改进欧拉法','x_1Runge-Kutta 法','x_2Runge-Kutta 法')
title('解的精度对比')