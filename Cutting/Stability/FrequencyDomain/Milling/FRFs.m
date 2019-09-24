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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using residues and modes superposition to %
% calculate FRFs                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference£º
%Altintas, Y. (2012). Manufacturing Automation: 
%Metal Cutting Mechanics, Machine Tool Vibrations, and CNC Design. 
%Applied Mechanics Reviews (Vol. 54). https://doi.org/10.1115/1.1399383
%P162
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
close all
clear

%% Modal parameters in X direction of cutting tool 
wnX   = [452.77 1448.53]*2*pi;               % rad/s
zetaX = [0.123718 0.01651];
resiX = [9.202966E-05-1j*1.862195E-04 -4.181562E-05-1j*3.043618E-04]; %m/N

% Modal parameters in Y direction of cutting tool
wnY   = [516.17 1407.64]*2*pi;               % rad/s
zetaY = [0.0243 0.0324];
resiY = [-2.3929E-06-1j*1.721539E-04 4.055052E-05-1j*3.618808E-04]; %m/N

%% Analysis frequency width and its resolution
f_start = 0;
f_end = 8000-0.5;
df = 1;
w = (f_start:df:f_end)'*2*pi;   % frequency, rad/s
f = w/2/pi;

%% Define cells

FRFiX={1,length(resiX)};
FRFiY={1,length(resiY)};

%% Initialization
FRFX = zeros(length(w),1);
FRFY = zeros(length(w),1);

alphaX = zeros(1,length(resiX));
alphaY = zeros(1,length(resiY));

betaX = zeros(1,length(resiX));
betaY = zeros(1,length(resiY));

wdX = zeros(1,length(resiX));
wdY = zeros(1,length(resiY));

sigmaX = real(resiX);
sigmaY = real(resiY);

nuX = imag(resiX);
nuY = imag(resiY);

%% Modes superposition
for cnt = 1:length(resiX)

    wdX(cnt)=wnX(cnt)*sqrt(1-zetaX(cnt)^2);
    alphaX(cnt)=2*(zetaX(cnt)*wnX(cnt)*sigmaX(cnt)-wdX(cnt)*nuX(cnt));
    betaX(cnt)=2*sigmaX(cnt);

    FRFiX{cnt} = (alphaX(cnt).*ones(length(w),1)+1j*betaX(cnt).*w)./(wnX(cnt)^2-w.^2 +1j* (2*zetaX(cnt)*wnX(cnt).*w));
    FRFX= FRFX + FRFiX{cnt};
end

for cnt = 1:length(resiY)

    wdY(cnt)=wnY(cnt)*sqrt(1-zetaY(cnt)^2);
    
    alphaY(cnt)=2*(zetaY(cnt)*wnY(cnt)*sigmaY(cnt)-wdY(cnt)*nuY(cnt));
    betaY(cnt)=2*sigmaY(cnt);

    FRFiY{cnt} = (alphaY(cnt).*ones(length(w),1)+1j*betaY(cnt).*w)./(wnY(cnt)^2-w.^2 +1j* (2*zetaY(cnt)*wnY(cnt).*w));
    FRFY= FRFY + FRFiY{cnt};
end

%% Plot lobes
figure(1)
subplot(211)
for cnt = 1:length(resiX)
    plot(f,real(FRFiX{cnt}))
    hold on
end
plot(f,real(FRFX),'b','linewidth',2)
xlim([0,8000])
grid on  
title('\fontsize{10}FRF X')
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
xlabel('\fontsize{10}\fontname{Times New Roman}\it f\rm/ Hz')
ylabel('\fontsize{10}\fontname{Times New Roman}\it FRF\rm/ m¡¤N^{-1}')

subplot(212)
for cnt = 1:length(resiX)
    plot(f,imag(FRFiX{cnt}))
    hold on
end
plot(f,imag(FRFX),'b','linewidth',2)
xlim([0,8000])
grid on  
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
set(gcf,'unit','centimeters','position',[0 10 13.53 9.03],'color','white');% word£¨13.5,9£©
xlabel('\fontsize{10}\fontname{Times New Roman}\it f\rm/ Hz')
ylabel('\fontsize{10}\fontname{Times New Roman}\it FRF\rm/ m¡¤N^{-1}')

figure(2)
for cnt = 1:length(resiX)
    plot(f,abs(FRFiX{cnt}))
    hold on
end
k11=plot(f,abs(FRFX),'b','linewidth',2);
xlim([0,8000])
grid on  
title('\fontsize{10}FRF X')
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
set(gcf,'unit','centimeters','position',[0 10 13.53 9.03],'color','white');% word£¨13.5,9£©
xlabel('\fontsize{10}\fontname{Times New Roman}\it f\rm/ Hz')
ylabel('\fontsize{10}\fontname{Times New Roman}\it FRF\rm/ m¡¤N^{-1}')
legend(k11,'X-Direction')

figure(3)
subplot(211)
for cnt = 1:length(resiY)
    plot(f,real(FRFiY{cnt}))
    hold on
end
plot(f,real(FRFY),'b','linewidth',2)
xlim([0,8000])
grid on  
title('\fontsize{10}FRF Y')
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
xlabel('\fontsize{10}\fontname{Times New Roman}\it f\rm/ Hz')
ylabel('\fontsize{10}\fontname{Times New Roman}\it FRF\rm/ m¡¤N^{-1}')

subplot(212)
for cnt = 1:length(resiY)
    plot(f,imag(FRFiY{cnt}))
    hold on
end
plot(f,imag(FRFY),'b','linewidth',2)
xlim([0,8000])
grid on  
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
set(gcf,'unit','centimeters','position',[15 10 13.53 9.03],'color','white');% word£¨13.5,9£©
xlabel('\fontsize{10}\fontname{Times New Roman}\it f\rm/ Hz')
ylabel('\fontsize{10}\fontname{Times New Roman}\it FRF\rm/ m¡¤N^{-1}')

figure(4)
for cnt = 1:length(resiY)
    plot(f,abs(FRFiY{cnt}))
    hold on
end
k22=plot(f,abs(FRFY),'b','linewidth',2);
xlim([0,8000])
grid on  
title('\fontsize{10}FRF Y')
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
set(gcf,'unit','centimeters','position',[15 10 13.53 9.03],'color','white');% word£¨13.5,9£©
xlabel('\fontsize{10}\fontname{Times New Roman}\it f\rm/ Hz')
ylabel('\fontsize{10}\fontname{Times New Roman}\it FRF\rm/ m¡¤N^{-1}')
legend(k22,'Y-Direction')
% 
% HX(:,1) = real(FRFX);
% HX(:,2) = imag(FRFX);
% 
% HY(:,1) = real(FRFY);
% HY(:,2) = imag(FRFY);

x_modes = length(resiX);
y_modes = length(resiY);

wnmax = max(max(wnX),max(wnY));

save('FRF','FRFX','FRFY','f_start','f_end','df','x_modes','y_modes','wnmax')