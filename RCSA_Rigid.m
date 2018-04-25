%% ͨ���������õ�������������Ƶ�죬��ͨ���ӽṹ����RCSA��ʹ��һ���������ϣ����������ɶ�����Ƶ��

clc
clear
close all

%% ������������Բ�����Ľṹ�Ͳ��ϲ�������ȡ��������ԭ�㼰�������Ƶ��

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

%% ����������ϵ�ǽ������ϵ�����Ƶ��
% ȫ��Ϊ0������Ϊ��ȫ����
h2b2b = zeros(1, length(w));
l2b2b = zeros(1, length(w));
n2b2b = zeros(1, length(w));
p2b2b = zeros(1, length(w));

%% ͨ���ӽṹ��Ϸ��õ��˵�Ƶ��
% ��ʼ����Ϻ�Ķ˵�Ƶ�����
H11 = zeros(1, length(w));
L11 = zeros(1, length(w));
N11 = zeros(1, length(w));
P11 = zeros(1, length(w));

% ������Ϻ��Ƶ��
for cnt = 1:length(w)
    % ����������������ԭ������Ƶ�����
    R11     = [h11(cnt)     l11(cnt);    n11(cnt)     p11(cnt)];
    R12a   = [h12a(cnt)   l12a(cnt);  n12a(cnt)    p12a(cnt)];
    R2a2a = [h2a2a(cnt) l2a2a(cnt); n2a2a(cnt) p2a2a(cnt)];
    R2a1   = [h2a1(cnt)   l2a1(cnt);  n2a1(cnt)   p2a1(cnt)];
    
% ����ǽ������Ƶ�����
    R2b2b = [h2b2b(cnt) l2b2b(cnt); n2b2b(cnt) p2b2b(cnt)];
    
% �������ɶ�����Ϻ��Ƶ�����
    G11 = R11 - R12a/(R2a2a + R2b2b)*R2a1;  % ע��˴���ֱ��*��/��������.* ��./
    
% ������Ϻ����ɶ�����Ƶ������еĸ�������
    H11(cnt) = G11(1,1);
    L11(cnt) = G11(1,2);
    N11(cnt) = G11(2,1);
    P11(cnt) = G11(2,2);
end

%% ֱ�Ӽ���һ�����ɡ�һ�˹̽���������ɶ�����Ƶ��
% ����һ�˹̽�һ�����ɵ�Բ��Ƶ��
H11cf  = -c1./(lambda.^3.*c8);
L11cf   = c2./(lambda.^2.*c8);
N11cf  = L11cf;
P11cf   = c5./(lambda.*c8);

%% ��ͼ
% ����Բ�����Ӳ���������״̬�µ�ĳһ�˵�λ��/��ԭ��Ƶ��
figure(1)
subplot(211)
plot(w/2/pi, real(h11),'r')  %Ƶ��f=w/(2*pi)
ylim([-1e-5 1e-5])
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
ylabel('\fontsize{10}\fontname{Times New Roman}Real / m��N^{-1}')
title('\fontsize{10}\fontname{����}����״̬��λ��/��ԭ��Ƶ��\fontname{Times New Roman}\ith\rm_{11}')

subplot(212)
plot(w/2/pi, imag(h11),'r')
ylim([-1.8e-5 1.8e-6])
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
xlabel('\fontsize{10}\fontname{Times New Roman}Frequency / Hz')
ylabel('\fontsize{10}\fontname{Times New Roman}Imag / m��N^{-1}')
set(gcf,'unit','centimeters','position',[28 5 13.53 9.03],'color','white');%��Ӧword��13.5,9��

% ������������ɶ˲���λ��/��ԭ��Ƶ�죬����һ�˹̽�һ�����ɵ�Բ���������ɶ˲�Ƶ�����Ա�
figure(2)
subplot(211)
plot(w/2/pi, real(H11), 'b', w/2/pi, real(H11cf), 'r:')
% axis([300 700 -3.75e-4 3.75e-4])
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
ylabel('\fontsize{10}\fontname{Times New Roman}Real / m��N^{-1}')
title('\fontsize{10}\fontname{����}��ϼ�ֱ�Ӽ�������λ��/��ԭ��Ƶ��')

subplot(212)
plot(w/2/pi, imag(H11), 'b', w/2/pi, imag(H11cf), 'r:')
% axis([300 700 -7e-4 7e-5])
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
xlabel('\fontsize{10}\fontname{Times New Roman}Frequency / Hz')
ylabel('\fontsize{10}\fontname{Times New Roman}Imag / m��N^{-1}')
legend('\fontsize{10}\fontname{����}�ӽṹ��','\fontsize{10}\fontname{����}ֱ�Ӽ��㷨')
set(gcf,'unit','centimeters','position',[28 5 13.53 9.03],'color','white');%��Ӧword��13.5,9��

% ���������
figure(3)
semilogy(w/2/pi, abs(H11), 'b', w/2/pi, abs(H11cf), 'r:')
% ylim([3e-9 1e-3])
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
xlabel('\fontsize{10}\fontname{Times New Roman}Frequency / Hz')
ylabel('\fontsize{10}\fontname{Times New Roman}Magnitude / m��N^{-1}')
legend('\fontsize{10}\fontname{����}�ӽṹ��','\fontsize{10}\fontname{����}ֱ�Ӽ��㷨')
title('\fontsize{10}\fontname{����}��ϼ�ֱ�Ӽ�������λ��/��ԭ��Ƶ��')
set(gcf,'unit','centimeters','position',[28 5 13.53 9.03],'color','white');%��Ӧword��13.5,9��

% ������������ɶ˲���ת��/Ť��ԭ��Ƶ�죬����һ�˹̽�һ�����ɵ�Բ���������ɶ˲�Ƶ�����Ա�
figure(4)
subplot(211)
plot(w/2/pi, real(P11), 'b', w/2/pi, real(P11cf), 'r:')
% axis([300 700 -0.045 0.045])
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
ylabel('\fontsize{10}\fontname{Times New Roman}Real / m��N^{-1}')
title('\fontsize{10}\fontname{����}��ϼ�ֱ�Ӽ�������ת��/����ԭ��Ƶ��')

subplot(212)
plot(w/2/pi, imag(P11), 'b', w/2/pi, imag(P11cf), 'r:')
% axis([300 700 -0.09 0.009])
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
xlabel('\fontsize{10}\fontname{Times New Roman}Frequency / Hz')
ylabel('\fontsize{10}\fontname{Times New Roman}Imag / m��N^{-1}')
legend('\fontsize{10}\fontname{����}�ӽṹ��','\fontsize{10}\fontname{����}ֱ�Ӽ��㷨')
set(gcf,'unit','centimeters','position',[28 5 13.53 9.03],'color','white');%��Ӧword��13.5,9��