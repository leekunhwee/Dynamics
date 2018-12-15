%% 通过解析法得到自由梁的两端频响，并通过子结构法（RCSA）使其一端与刚体耦合，并求其自由端完整频响

clc
clear
close all

%% 定义两端自由圆柱梁的结构和材料参数，获取梁的两端原点及跨点完整频响

w = (1:0.1:5000)*2*pi;           % frequency, rad/s
E = 200e9;                             % elastic modulus, N/m^2
d = 10e-3;                             % diameter, m
L = 125e-3;                           % length, m
I  = pi*d^4/64;                      % 2nd moment of area, m^4
rho = 7800;                           % density, kg/m^3
A = pi*d^2/4;                       % cross sectional area, m^2
eta = 0.01;                             % solid damping factor
EI = E*I*(1+1i*eta);               % complex stiffness, N-m^2
lambda = (w.^2*rho*A/EI).^0.25;

c1 = cos(lambda*L).*sinh(lambda*L) - sin(lambda*L).*cosh(lambda*L);
c2 = sin(lambda*L).*sinh(lambda*L);
c3 = sin(lambda*L) - sinh(lambda*L);
c4 = cos(lambda*L) - cosh(lambda*L);
c5 = cos(lambda*L).*sinh(lambda*L) + sin(lambda*L).*cosh(lambda*L);
c6 = sin(lambda*L) + sinh(lambda*L);
c7 = EI*(cos(lambda*L).*cosh(lambda*L)-1);
c8 = EI*(cos(lambda*L).*cosh(lambda*L)+1);

h11 = -c1./(lambda.^3.*c7);
l11  = c2./(lambda.^2.*c7);
n11 = l11;
p11 = c5./(lambda.*c7);

h2a2a = -c1./(lambda.^3.*c7);
l2a2a  = -c2./(lambda.^2.*c7);
n2a2a = l2a2a;
p2a2a = c5./(lambda.*c7);

h12a = c3./(lambda.^3.*c7);
l12a  = -c4./(lambda.^2.*c7);
n12a = c4./(lambda.^2.*c7);
p12a = c6./(lambda.*c7);

h2a1 = h12a;
l2a1  = n12a;
n2a1 = l12a;
p2a1 = p12a;

%% 定义与梁耦合的墙壁上耦合点完整频响
% 全设为0，即认为完全刚性
h2b2b = zeros(1, length(w));
l2b2b = zeros(1, length(w));
n2b2b = zeros(1, length(w));
p2b2b = zeros(1, length(w));

%% 通过子结构耦合法得到端点频响
% 初始化耦合后的端点频响矩阵
H11 = zeros(1, length(w));
L11 = zeros(1, length(w));
N11 = zeros(1, length(w));
P11 = zeros(1, length(w));

% 计算耦合后的频响
for cnt = 1:length(w)
    % 两端自由梁的两端原点与跨点频响矩阵
    R11     = [h11(cnt)     l11(cnt);    n11(cnt)     p11(cnt)];
    R12a   = [h12a(cnt)   l12a(cnt);  n12a(cnt)    p12a(cnt)];
    R2a2a = [h2a2a(cnt) l2a2a(cnt); n2a2a(cnt) p2a2a(cnt)];
    R2a1   = [h2a1(cnt)   l2a1(cnt);  n2a1(cnt)   p2a1(cnt)];
    
% 定义墙的完整频响矩阵
    R2b2b = [h2b2b(cnt) l2b2b(cnt); n2b2b(cnt) p2b2b(cnt)];
    
% 生成自由端在耦合后的频响矩阵
    G11 = R11 - R12a/(R2a2a + R2b2b)*R2a1;  % 注意此处是直接*和/，而不是.* 和./
    
% 更新耦合后自由端完整频响矩阵中的各个分量
    H11(cnt) = G11(1,1);
    L11(cnt) = G11(1,2);
    N11(cnt) = G11(2,1);
    P11(cnt) = G11(2,2);
end

%% 直接计算一端自由、一端固结的梁的自由端完整频响
% 定义一端固结一端自由的圆柱频响
H11cf  = -c1./(lambda.^3.*c8);
L11cf   = c2./(lambda.^2.*c8);
N11cf  = L11cf;
P11cf   = c5./(lambda.*c8);

%% 绘图
% 绘制圆柱梁子部件在自由状态下的某一端的位移/力原点频响
figure(1)
subplot(211)
plot(w/2/pi, real(h11),'r')  %频率f=w/(2*pi)
ylim([-1e-5 1e-5])
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
ylabel('\fontsize{10}\fontname{Times New Roman}Real / m·N^{-1}')
title('\fontsize{10}\fontname{宋体}自由状态下位移/力原点频响\fontname{Times New Roman}\ith\rm_{11}')

subplot(212)
plot(w/2/pi, imag(h11),'r')
ylim([-1.8e-5 1.8e-6])
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
xlabel('\fontsize{10}\fontname{Times New Roman}Frequency / Hz')
ylabel('\fontsize{10}\fontname{Times New Roman}Imag / m·N^{-1}')
set(gcf,'unit','centimeters','position',[28 5 13.53 9.03],'color','white');%对应word（13.5,9）

% 绘制组合体自由端部的位移/力原点频响，并与一端固结一端自由的圆柱梁的自由端部频响作对比
figure(2)
subplot(211)
plot(w/2/pi, real(H11), 'b', w/2/pi, real(H11cf), 'r:')
% axis([300 700 -3.75e-4 3.75e-4])
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
ylabel('\fontsize{10}\fontname{Times New Roman}Real / m·N^{-1}')
title('\fontsize{10}\fontname{宋体}耦合及直接计算所得位移/力原点频响')

subplot(212)
plot(w/2/pi, imag(H11), 'b', w/2/pi, imag(H11cf), 'r:')
% axis([300 700 -7e-4 7e-5])
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
xlabel('\fontsize{10}\fontname{Times New Roman}Frequency / Hz')
ylabel('\fontsize{10}\fontname{Times New Roman}Imag / m·N^{-1}')
legend('\fontsize{10}\fontname{宋体}子结构法','\fontsize{10}\fontname{宋体}直接计算法')
set(gcf,'unit','centimeters','position',[28 5 13.53 9.03],'color','white');%对应word（13.5,9）

% 半对数坐标
figure(3)
semilogy(w/2/pi, abs(H11), 'b', w/2/pi, abs(H11cf), 'r:')
% ylim([3e-9 1e-3])
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
xlabel('\fontsize{10}\fontname{Times New Roman}Frequency / Hz')
ylabel('\fontsize{10}\fontname{Times New Roman}Magnitude / m·N^{-1}')
legend('\fontsize{10}\fontname{宋体}子结构法','\fontsize{10}\fontname{宋体}直接计算法')
title('\fontsize{10}\fontname{宋体}耦合及直接计算所得位移/力原点频响')
set(gcf,'unit','centimeters','position',[28 5 13.53 9.03],'color','white');%对应word（13.5,9）

% 绘制组合体自由端部的转角/扭矩原点频响，并与一端固结一端自由的圆柱梁的自由端部频响作对比
figure(4)
subplot(211)
plot(w/2/pi, real(P11), 'b', w/2/pi, real(P11cf), 'r:')
% axis([300 700 -0.045 0.045])
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
ylabel('\fontsize{10}\fontname{Times New Roman}Real / m·N^{-1}')
title('\fontsize{10}\fontname{宋体}耦合及直接计算所得转角/力矩原点频响')

subplot(212)
plot(w/2/pi, imag(P11), 'b', w/2/pi, imag(P11cf), 'r:')
% axis([300 700 -0.09 0.009])
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
xlabel('\fontsize{10}\fontname{Times New Roman}Frequency / Hz')
ylabel('\fontsize{10}\fontname{Times New Roman}Imag / m·N^{-1}')
legend('\fontsize{10}\fontname{宋体}子结构法','\fontsize{10}\fontname{宋体}直接计算法')
set(gcf,'unit','centimeters','position',[28 5 13.53 9.03],'color','white');%对应word（13.5,9）