%% ͨ���������õ�������������Ƶ�죬��ͨ���ӽṹ����RCSA��ʹ��һ����һ���Ǹ�����ϣ����������ɶ�����Ƶ��
% һ��Բ������һ��������ϣ��÷�������ǽ�������
clc
clear
close all

%% ������������Բ�����Ľṹ�Ͳ��ϲ�������ȡԲ����������ԭ�㼰�������Ƶ��

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

%% ����һ������һ�˹᷽̽���Ľṹ�Ͳ��ϲ�������ȡ���������ɶ˵�����Ƶ��

E = 200e9;                  % elastic modulus, N/m^2
s = 50e-3;                  % square side, m
L = 200e-3;                 % length, m
% �������ĳ��ȱ�Ϊ250e-3�������նȽ��ͣ���һ�׹���Ƶ�ʽ��ӽ���Ϻ��Բ������һ�׹���Ƶ��
% ������Բ�������ɶ˱��ֳ����նȷ������󣬴�Ϊ���������񡱡�������������֮������ȷֲ�
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

%% ͨ���ӽṹ��Ϸ��õ�Բ�����Ķ˵�Ƶ��
% ��ʼ����Ϻ�Ķ˵�Ƶ�����
H11 = zeros(1, length(w));
L11  = zeros(1, length(w));
N11 = zeros(1, length(w));
P11  = zeros(1, length(w));

for cnt = 1:length(w)
    % Բ��������״̬�µ�����ԭ�㼰�������Ƶ��
    R11     = [h11(cnt) l11(cnt); n11(cnt) p11(cnt)];
    R12a   = [h12a(cnt) l12a(cnt); n12a(cnt) p12a(cnt)];
    R2a2a = [h2a2a(cnt) l2a2a(cnt); n2a2a(cnt) p2a2a(cnt)];
    R2a1   = [h2a1(cnt) l2a1(cnt); n2a1(cnt) p2a1(cnt)];
    
    % ���巽�����ɶ���ϵ㴦������Ƶ�����
    R2b2b = [h2b2b(cnt) l2b2b(cnt); n2b2b(cnt) p2b2b(cnt)];
    
    % ������Ϻ��Ƶ��    
    G11 = R11 - R12a/(R2a2a + R2b2b)*R2a1;
    
    % ������Ϻ����ɶ�����Ƶ������еĸ�������
    H11(cnt) = G11(1,1);
    L11(cnt)  = G11(1,2);
    N11(cnt) = G11(2,1);
    P11(cnt)  = G11(2,2);
end

%% ��ͼ
% ���ǰԲ�����ͷ��������ɶ˲�λ��/��Ƶ��
figure(1)
subplot(211)
plot(w/2/pi, real(h11), 'b', w/2/pi, real(h2b2b), 'r:')
ylim([-5e-6 5e-6])
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
ylabel('\fontsize{10}\fontname{Times New Roman}Real / m��N^{-1}')
title('\fontsize{10}\fontname{����}���ǰ���ɶ˲�λ��/��Ƶ��')

subplot(212)
plot(w/2/pi, imag(h11), 'b', w/2/pi, imag(h2b2b), 'r:')
ylim([-9e-6 9e-7])
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
xlabel('\fontsize{10}\fontname{Times New Roman}Frequency / Hz')
ylabel('\fontsize{10}\fontname{Times New Roman}Imag / m��N^{-1}')
kk=legend('\fontsize{10}\fontname{����}Բ����','\fontsize{10}\fontname{����}����');
set(kk,'Position',[0.590119045268922 0.195634918182615 0.17178571663584 0.113571431023734]);
set(gcf,'unit','centimeters','position',[28 5 13.53 9.03],'color','white');%��Ӧword��13.5,9��

% ��Ϻ�Բ���������ɶ˲�λ��/��Ƶ��
figure(2)
subplot(211)
plot(w/2/pi, real(H11))
axis([500 1200 -2.0e-4 2.0e-4])
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
ylabel('\fontsize{10}\fontname{Times New Roman}Real / m��N^{-1}')
title('\fontsize{10}\fontname{����}��Ϻ�Բ�������ɶ˲�λ��/��Ƶ��\fontname{Times New Roman}\itH\rm11')

subplot(212)
plot(w/2/pi, imag(H11))
axis([500 1200 -3.5e-4 5e-5])
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
xlabel('\fontsize{10}\fontname{Times New Roman}Frequency / Hz')
ylabel('\fontsize{10}\fontname{Times New Roman}Imag / m��N^{-1}')
set(gcf,'unit','centimeters','position',[28 5 13.53 9.03],'color','white');%��Ӧword��13.5,9��