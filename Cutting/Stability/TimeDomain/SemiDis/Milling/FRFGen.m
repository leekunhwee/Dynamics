clc
close all
clear

% Define mode parameters
zetaX = [0.0316 0.0114];
wnX = [369.85 1146.36]*2*pi;               % rad/s
kX = [wnX.*wnX];                     % N/m

zetaY = [1.68*0.01 1.55*0.01];
wnY = [411.21 1069.21]*2*pi;               % rad/s
kY = [wnY.*wnY];                     % N/m
f_start=0;
f_end=3000-0.5;
df=0.5;
w = (f_start:df:f_end)'*2*pi;   % frequency, rad/s
f=w/2/pi;

% FRFs synthesis
FRFY = (wnY(1)^2/kY(1))./(wnY(1)^2 - w.^2 + 1i*2*zetaY(1)*wnY(1).*w);
for cnt = 2:length(kY)
    FRFY = FRFY + (wnY(cnt)^2/kY(cnt))./(wnY(cnt)^2 - w.^2 + 1i*2*zetaY(cnt)*wnY(cnt).*w);
end
FRFX = (wnX(1)^2/kX(1))./(wnX(1)^2 - w.^2 + 1i*2*zetaX(1)*wnX(1).*w);
for cnt = 2:length(kX)
    FRFX = FRFX + (wnX(cnt)^2/kX(cnt))./(wnX(cnt)^2 - w.^2 + 1i*2*zetaX(cnt)*wnX(cnt).*w);
end
FRFMagX=sqrt(real(FRFX).^2+imag(FRFX).^2);
FRFMagY=sqrt(real(FRFY).^2+imag(FRFY).^2);

% Plot FRFs
% figure(1)
% plot(f,FRFMagY,'b-','linewidth',2)
% hold on
% xlim([200,2200])
% grid on  
% legend('\fontsize{10}\fontname{Times New Roman}\itFRF','location','northwest')
% set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
% set(gcf,'unit','centimeters','position',[0 10 13.53 9.03],'color','white');% word?°ß13.5,9??
% title('\fontsize{10}FRF')
% xlabel('\fontsize{10}\fontname{Times New Roman}\it f\rm/ Hz')
% ylabel('\fontsize{10}\fontname{Times New Roman}\it FRF\rm/ m?°ËN^{-1}')
% hold on

h1X=real(FRFX);
h2X=imag(FRFX);
figure(1)

subplot(211)
plot(f,h1X,'b-','linewidth',2)
hold on
xlim([200,2200])
grid on  
legend('\fontsize{10}\fontname{Times New Roman}\itRea;','location','northwest')
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
set(gcf,'unit','centimeters','position',[0 10 13.53 9.03],'color','white');% word?°ß13.5,9??
title('\fontsize{10}FRF')
xlabel('\fontsize{10}\fontname{Times New Roman}\it f\rm/ Hz')
ylabel('\fontsize{10}\fontname{Times New Roman}\it FRF\rm/ m?°ËN^{-1}')

subplot(212)
plot(f,h2X,'b-','linewidth',2)
hold on
xlim([200,2200])
grid on  
legend('\fontsize{10}\fontname{Times New Roman}\itImag','location','northwest')
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
set(gcf,'unit','centimeters','position',[0 10 13.53 9.03],'color','white');% word?°ß13.5,9??
title('\fontsize{10}FRF')
xlabel('\fontsize{10}\fontname{Times New Roman}\it f\rm/ Hz')
ylabel('\fontsize{10}\fontname{Times New Roman}\it FRF\rm/ m?°ËN^{-1}')


h1Y=real(FRFY);
h2Y=imag(FRFY);

figure(2)
subplot(211)
plot(f,h1Y,'b-','linewidth',2)
hold on
xlim([200,2200])
grid on  
legend('\fontsize{10}\fontname{Times New Roman}\itRea;','location','northwest')
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
set(gcf,'unit','centimeters','position',[0 10 13.53 9.03],'color','white');% word?°ß13.5,9??
title('\fontsize{10}FRF')
xlabel('\fontsize{10}\fontname{Times New Roman}\it f\rm/ Hz')
ylabel('\fontsize{10}\fontname{Times New Roman}\it FRF\rm/ m?°ËN^{-1}')

subplot(212)
plot(f,h2Y,'b-','linewidth',2)
hold on
xlim([200,2200])
grid on  
legend('\fontsize{10}\fontname{Times New Roman}\itImag','location','northwest')
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
set(gcf,'unit','centimeters','position',[0 10 13.53 9.03],'color','white');% word?°ß13.5,9??
title('\fontsize{10}FRF')
xlabel('\fontsize{10}\fontname{Times New Roman}\it f\rm/ Hz')
ylabel('\fontsize{10}\fontname{Times New Roman}\it FRF\rm/ m?°ËN^{-1}')

% Save as a .mat files, name is H00, including H00 f_start f_end df
HX00=[h1X,h2X];
HY00=[h1Y,h2Y];

save ('H00','HX00','HY00','f_start','f_end','df')