%% 通过解析法得到自由梁的两端频响，并通过子结构法（RCSA）使其一端与一个非刚体耦合，并求其自由端完整频响
% 一个圆柱梁与一个方梁耦合，该方梁又与墙刚性耦合
clc
clear
close all

%% 定义两端自由圆柱梁的结构和材料参数，获取圆柱梁的两端原点及跨点完整频响

w = (1:0.1:5000)*2*pi;      % frequency, rad/s
E = 200e9;                        % elastic modulus, N/m^2
d = 10e-3;                        % diameter, m
L = 100e-3;                       % length, m
I  = pi*d^4/64;                 % 2nd moment of area, m^4
rho = 7800;                       % density, kg/m^3
A = pi*d^2/4;                   % cross sectional area, m^2
eta = 0.01;                         % solid damping factor
EI = E*I*(1+1i*eta);           % complex stiffness, N-m^2
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

%% 定义一端自由一端固结方梁的结构和材料参数，获取方梁的自由端的完整频响

E = 200e9;                  % elastic modulus, N/m^2
s = 50e-3;                  % square side, m
L = 200e-3;                 % length, m
% 将方梁的长度变为250e-3，则方梁刚度降低，其一阶固有频率将接近耦合后的圆柱梁的一阶固有频率
% 导致在圆柱梁自由端表现出来刚度反而增大，此为“动力吸振”。而能量在两者之间更均匀分布
I = s^4/12;                 % 2nd moment of area, m^4
rho = 7800;                 % density, kg/m^3
A = s^2;                    % cross sectional area, m^2
eta = 0.01;                 % solid damping factor
EI = E*I*(1+1i*eta);         % complex stiffness, N-m^2
lambda = (w.^2*rho*A/EI).^0.25;

c1 = cos(lambda*L).*sinh(lambda*L) - sin(lambda*L).*cosh(lambda*L);
c2 = sin(lambda*L).*sinh(lambda*L);
c3 = sin(lambda*L) - sinh(lambda*L);
c4 = cos(lambda*L) - cosh(lambda*L);
c5 = cos(lambda*L).*sinh(lambda*L) + sin(lambda*L).*cosh(lambda*L);
c6 = sin(lambda*L) + sinh(lambda*L);
c7 = EI*(cos(lambda*L).*cosh(lambda*L)-1);
c8 = EI*(cos(lambda*L).*cosh(lambda*L)+1);

h2b2b = -c1./(lambda.^3.*c8);
l2b2b  = c2./(lambda.^2.*c8);
n2b2b = l2b2b;
p2b2b = c5./(lambda.*c8);

%% 通过子结构耦合法得到圆柱梁的端点频响
% 初始化耦合后的端点频响矩阵
H11 = zeros(1, length(w));
L11  = zeros(1, length(w));
N11 = zeros(1, length(w));
P11  = zeros(1, length(w));

for cnt = 1:length(w)
    % 圆柱梁自由状态下的两端原点及跨点完整频响
    R11     = [h11(cnt) l11(cnt); n11(cnt) p11(cnt)];
    R12a   = [h12a(cnt) l12a(cnt); n12a(cnt) p12a(cnt)];
    R2a2a = [h2a2a(cnt) l2a2a(cnt); n2a2a(cnt) p2a2a(cnt)];
    R2a1   = [h2a1(cnt) l2a1(cnt); n2a1(cnt) p2a1(cnt)];
    
    % 定义方梁自由端耦合点处的完整频响矩阵
    R2b2b = [h2b2b(cnt) l2b2b(cnt); n2b2b(cnt) p2b2b(cnt)];
    
    % 计算耦合后的频响    
    G11 = R11 - R12a/(R2a2a + R2b2b)*R2a1;
    
    % 更新耦合后自由端完整频响矩阵中的各个分量
    H11(cnt) = G11(1,1);
    L11(cnt)  = G11(1,2);
    N11(cnt) = G11(2,1);
    P11(cnt)  = G11(2,2);
end

%% 绘图
% 耦合前圆柱梁和方梁的自由端部位移/力频响
figure(1)
subplot(211)
plot(w/2/pi, real(h11), 'b', w/2/pi, real(h2b2b), 'r:')
ylim([-5e-6 5e-6])
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
ylabel('\fontsize{10}\fontname{Times New Roman}Real / m·N^{-1}')
title('\fontsize{10}\fontname{宋体}耦合前自由端部位移/力频响')

subplot(212)
plot(w/2/pi, imag(h11), 'b', w/2/pi, imag(h2b2b), 'r:')
ylim([-9e-6 9e-7])
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
xlabel('\fontsize{10}\fontname{Times New Roman}Frequency / Hz')
ylabel('\fontsize{10}\fontname{Times New Roman}Imag / m·N^{-1}')
kk=legend('\fontsize{10}\fontname{宋体}圆柱梁','\fontsize{10}\fontname{宋体}方梁');
set(kk,'Position',[0.590119045268922 0.195634918182615 0.17178571663584 0.113571431023734]);
set(gcf,'unit','centimeters','position',[28 5 13.53 9.03],'color','white');%对应word（13.5,9）

% 耦合后圆柱梁的自由端部位移/力频响
figure(2)
subplot(211)
plot(w/2/pi, real(H11))
axis([500 1200 -2.0e-4 2.0e-4])
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
ylabel('\fontsize{10}\fontname{Times New Roman}Real / m·N^{-1}')
title('\fontsize{10}\fontname{宋体}耦合后圆柱梁自由端部位移/力频响\fontname{Times New Roman}\itH\rm11')

subplot(212)
plot(w/2/pi, imag(H11))
axis([500 1200 -3.5e-4 5e-5])
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
xlabel('\fontsize{10}\fontname{Times New Roman}Frequency / Hz')
ylabel('\fontsize{10}\fontname{Times New Roman}Imag / m·N^{-1}')
set(gcf,'unit','centimeters','position',[28 5 13.53 9.03],'color','white');%对应word（13.5,9）