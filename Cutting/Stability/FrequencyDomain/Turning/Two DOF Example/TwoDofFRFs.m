%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------Copyright--------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Jianhui Li                         %
% Time: 09/22/2019                           %
% University of British Columbia, BC, Canada %
% Affiliation:                               %
% Department of Mechanical Engineering       %
% Manufacturing Automation Laboratary        %
% E-mail: jianhui.li@alumni.ubc.ca           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference£º
%Altintas, Y. (2012). Manufacturing Automation: 
%Metal Cutting Mechanics, Machine Tool Vibrations, and CNC Design. 
%Applied Mechanics Reviews (Vol. 54). https://doi.org/10.1115/1.1399383
%P131
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define modal parmeters
k    = [2.26e8 2.13e8];               % Modal stiffness N/m
zeta = [0.012  0.010 ];               % Modal damping ratio
wn   = [250    150   ]*2*pi;          % Natural frequency rad/s

f_start= 0;
f_end  = 300-0.5;
df     = 0.1;
w      = (f_start:df:f_end)*2*pi;   % frequency, rad/s
f      = w/2/pi;

%% FRF Synthesis (To one direction)
FRF1 = cos( 30/180*pi)^2*(wn(1)^2/k(1))./(wn(1)^2 - w.^2 + 1i*2*zeta(1)*wn(1).*w);
FRF2 = cos(-45/180*pi)^2*(wn(2)^2/k(2))./(wn(2)^2 - w.^2 + 1i*2*zeta(2)*wn(2).*w);
FRF  = FRF1 + FRF2;

% for cnt = 2:length(k)
%     FRF = FRF + (wn(cnt)^2/k(cnt))./(wn(cnt)^2 - w.^2 + 1i*2*zeta(cnt)*wn(cnt).*w);
% end

%% Plot FRF Magnitude and Phase
figure(1)
subplot(211)
plot(f,abs(FRF),'b-','linewidth',2)
hold on
xlim([0,300])
grid on  
% legend('\fontsize{10}\fontname{Times New Roman}\itFRF','location','northwest')
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
ylabel('\fontsize{10}\fontname{Times New Roman}\it FRF\rm/ m.N^{-1}')
subplot(212)
plot(f,angle(FRF)*180/pi,'b-','linewidth',2)
xlim([0,300])
ylim([-180,0])
set(gca,'ytick',[-180 -140 -100 -60 -20 0])
grid on 
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
set(gcf,'unit','centimeters','position',[0 10 13.53 9.03],'color','white');
title('\fontsize{10}FRF')
xlabel('\fontsize{10}\fontname{Times New Roman}\it f\rm/ Hz')
ylabel('\fontsize{10}\fontname{Times New Roman} degree/ ^\circ')
hold on

%% Plot FRF Real-Image Parts
figure(2)
subplot(211)
plot(f,real(FRF),'b-','linewidth',2)
hold on
grid on  
title('\fontsize{10}FRF')
xlim([0,300])
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
ylabel('\fontsize{10}\fontname{Times New Roman}\it FRF\rm/ m.N^{-1}')
subplot(212)
plot(f,imag(FRF),'b-','linewidth',2)
xlim([0,300])
grid on  
% legend('\fontsize{10}\fontname{Times New Roman}\itFRF','location','northwest')
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
set(gcf,'unit','centimeters','position',[15 3 13.53 9.03],'color','white');
xlabel('\fontsize{10}\fontname{Times New Roman}\it f\rm/ Hz')
ylabel('\fontsize{10}\fontname{Times New Roman}\it FRF\rm/ m.N^{-1}')

save('FRF','FRF','w','f_start','f_end','df')