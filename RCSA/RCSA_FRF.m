%% Ƶ�����

clc
clear
close all

%% Ƶ����Ϸ��õ������ô���Ƶ��
% ���������ӽṹ�ھ��������µ��������
m1 = 3;             % kg
c1 = 200;           % N-s/m
k1 = 2e6;           % N/m

m2 = 2;             % kg
c2 = 100;           % N-s/m
k2 = 1e6;           % N/m

% ��ϲ�λ�ĸն�
kc = 5e5;           % N/m

% ���������ӽṹ���Ե�λ��/��Ƶ��
w = (0:0.1:300)'*2*pi; % frequency, rad/s
FRF_I  = 1./(-w.^2*m1 + 1i*w*c1 + k1); % �ṹI��Ƶ�캯��
FRF_II = 1./(-w.^2*m2 + 1i*w*c2 + k2); % �ṹII��Ƶ�캯��

% Ƶ����Ϸ��õ�����Ϊ������������2��
FRF_III_rc = FRF_II - FRF_II./(FRF_I + FRF_II + 1/kc).*FRF_II;

%% ģ̬�������õ�Ƶ��
% ���s^2�ĸ�
s_squared = roots([(m1*m2) (m1*(k2+kc)+m2*(k1+kc)) ((k1+kc)*(k2+kc)-kc^2)]);
s1_squared = s_squared(1);
s2_squared = s_squared(2);
% ����ֵ����Ŀ����ʹ�� wn1 < wn2
if s1_squared < s2_squared
    temp = s1_squared;
    s1_squared = s2_squared;
    s2_squared = temp;
end
wn1 = sqrt(-s1_squared);
wn2 = sqrt(-s2_squared);

p1 = kc/(m1*s1_squared + (k1+kc));
p2 = kc/(m1*s2_squared + (k1+kc));

% ��������ϵ�µĶ���ѧ��������
M = [m1 0; 0 m2];
C = [c1 0; 0 c2];
K = [k1+kc -kc; -kc k2+kc];

% ģ̬�����µ�ģ̬����
P = [p1 p2; 1 1];
mq = P'*M*P;
cq = P'*C*P;
kq = P'*K*P;

mq1 = mq(1,1);
mq2 = mq(2,2);
cq1  = cq(1,1);
cq2  = cq(2,2);
kq1  = kq(1,1);
kq2  = kq(2,2);

zetaq1 = cq1/(2*sqrt(kq1*mq1));
zetaq2 = cq2/(2*sqrt(kq2*mq2));

r1 = w/wn1;
r2 = w/wn2;

% ģ̬����
FRF_III_modal = 1/kq1*((1-r1.^2) - 1i*(2*zetaq1*r1))./((1-r1.^2).^2 + (2*zetaq1*r1).^2) + 1/kq2*((1-r2.^2) - 1i*(2*zetaq2*r2))./((1-r2.^2).^2 + (2*zetaq2*r2).^2);

%% Ƶ�캯��������
% ���������淨�õ�Ƶ�캯�����󣬼���ν��Ƶ�캯����������
FRF_III_inversion = (-w.^2*m1 + 1i*w*c1 + k1 + kc)./((-w.^2*m1 + 1i*w*c1 + k1 + kc).*(-w.^2*m2 + 1i*w*c2 + k2 + kc) - kc^2);

%% ��ͼ
% �ֱ�����ӽṹ���Ե�Ƶ�캯�����ӽṹ��Ϸ��õ���Ƶ�캯����ģ̬�������õ���Ƶ�캯���Լ�Ƶ�캯���������õ���Ƶ�캯��
figure(1)
subplot(211)
plot(w/2/pi, real(FRF_I), 'b:', w/2/pi, real(FRF_II), 'b-.', w/2/pi, real(FRF_III_rc), 'b', w/2/pi, real(FRF_III_modal), 'g', w/2/pi, real(FRF_III_inversion), 'r')
ylim([-8e-6 9e-6])
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
ylabel('\fontsize{10}\fontname{Times New Roman}Real / m��N^{-1}')
title('\fontsize{10}\fontname{����}�����õ㴦ԭ��Ƶ��')

subplot(212)
plot(w/2/pi, imag(FRF_I), 'b:', w/2/pi, imag(FRF_II), 'b-.', w/2/pi, imag(FRF_III_rc), 'b', w/2/pi, imag(FRF_III_modal), 'g', w/2/pi, imag(FRF_III_inversion), 'r')
ylim([-16e-6 16e-7])
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
xlabel('\fontsize{10}\fontname{Times New Roman}Frequency / Hz')
ylabel('\fontsize{10}\fontname{Times New Roman}Imag / m��N^{-1}')
legend('\fontsize{10}\fontname{����}�ṹIƵ��','\fontsize{10}\fontname{����}�ṹIIƵ��','\fontsize{10}\fontname{����}��Ϸ�','\fontsize{10}\fontname{����}ģ̬������','\fontsize{10}\fontname{����}Ƶ�������')
legend boxoff
% �������
% ��������ΪƵ�캯���������õ����Ǿ�ȷ���
diff_modal = FRF_III_inversion - FRF_III_modal;
diff_rc = FRF_III_inversion - FRF_III_rc;
set(gcf,'unit','centimeters','position',[28 5 13.53 9.03],'color','white');%��Ӧword��13.5,9��

figure(2)
subplot(211)
plot(w/2/pi, real(diff_modal), 'b', w/2/pi, imag(diff_modal), 'r:')
kk1=legend('\fontsize{10}\fontname{Times New Roman}Re', '\fontsize{10}\fontname{Times New Roman}Im');
set(kk1,'Position',[0.792738094982647 0.622164313899322 0.0946428576537541 0.113571431023734]);
ylim([-2e-7 2e-7])
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
ylabel('\fontsize{10}\fontname{Times New Roman}Modal errors / m��N^{-1}')
title('\fontsize{10}\fontname{����}���ַ��������Ա�')
legend boxoff
subplot(212)
plot(w/2/pi, real(diff_rc), 'b', w/2/pi, imag(diff_rc), 'r:')
kk2=legend('\fontsize{10}\fontname{Times New Roman}Re', '\fontsize{10}\fontname{Times New Roman}Im');
set(kk2,'Position',[0.794523809268361 0.15213662840652 0.0946428576537541 0.113571431023734]);
ylim([-6e-21 4e-21])
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
xlabel('\fontsize{10}\fontname{Times New Roman}Frequency / Hz')
ylabel('\fontsize{10}\fontname{Times New Roman}Receptance errors / m��N^{-1}')
legend boxoff
set(gcf,'unit','centimeters','position',[28 5 13.53 9.03],'color','white');%��Ӧword��13.5,9��