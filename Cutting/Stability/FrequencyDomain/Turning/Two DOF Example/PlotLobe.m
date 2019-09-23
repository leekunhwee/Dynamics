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

load FRF

% Define force parameters
Kf=1000e6;          % N/m^2

%% Determine valid chatter frequency range 
% The chatter may happen when the real part of FRF is negtive
FRF_real = real(FRF);
FRF_imag = imag(FRF);
index = find(FRF_real < 0);
FRF_real = FRF_real(index);
FRF_imag = FRF_imag(index);
w = w(index);

% To illustrate the effects of diffrent modes, separate the two modes 
cou=0;
for cnt=1:length(index)-1
    if index(cnt+1)-index(cnt)==1
        continue;
    else
        cou=cnt;
        break;
    end
end

FRF_real1=FRF_real(1:cou);
FRF_imag1=FRF_imag(1:cou);
FRF_real2=FRF_real(cou+1:end);
FRF_imag2=FRF_imag(cou+1:end);
w1=w(1:cou);
w2=w(cou+1:end);
 
%% Plot Lobe
% Calculate blim
blim = -1./(2*Kf*FRF_real1);  % m
blim = blim*1e3;        % convert to mm

% Calculate epsilon which retlated to the rotation speed of spindle
for cnt = 1:length(FRF_imag1)
    if FRF_imag1(cnt)<0
        epsilon(cnt) = 2*pi-2*atan(abs(FRF_real1(cnt)/FRF_imag1(cnt)));
    else
        epsilon(cnt) = pi-2*atan(abs(FRF_imag1(cnt)/FRF_real1(cnt)));
    end
end

% Calculate spindle speeds for N = 0 to 2
omega0 = w1./(0*2*pi + epsilon);   % rps
omega1 = w1./(1*2*pi + epsilon);
omega2 = w1./(2*2*pi + epsilon);   % rps
omega3 = w1./(3*2*pi + epsilon);

figure(3)
plot(omega0*60, blim, 'r-.', omega1*60, blim, 'r-.', omega2*60, blim, 'r-.', omega3*60, blim, 'r-.','linewidth',2)
axis([2000 15000 0 60])
set(gca,'FontSize', 14)
xlabel('\fontsize{10}\fontname{Times New Roman} \Omega / rpm')
ylabel('\fontsize{10}\fontname{Times New Roman}\itb_{lim} \rm/ mm')
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
set(gcf,'unit','centimeters','position',[15 15 13.53 9.03],'color','white');
title('\fontsize{10}Lobe Diagram')
hold on

% Calculate blim
blim = -1./(2*Kf*FRF_real2);  % m
blim = blim*1e3;        % convert to mm
% Calculate epsilon
for cnt = 1:length(FRF_imag2)
    if FRF_imag2(cnt) < 0
        epsilon(cnt) = 2*pi - 2*atan(abs(FRF_real2(cnt)/FRF_imag2(cnt)));
    else
        epsilon(cnt) = pi - 2*atan(abs(FRF_imag2(cnt)/FRF_real2(cnt)));
    end
end

% Calculate spindle speeds for N = 0 to 2
omega0 = w2./(0*2*pi + epsilon);   % rps
omega1 = w2./(1*2*pi + epsilon);
omega2 = w2./(2*2*pi + epsilon);   % rps
omega3 = w2./(3*2*pi + epsilon);
plot(omega0*60, blim, 'b', omega1*60, blim, 'b', omega2*60, blim, 'b', omega3*60, blim, 'b','linewidth',2)
grid on

annotation('textarrow',[.597,.57],[.75,.68],'String','150Hz','FontName', 'Times New Roman'); %[x1,x2],[y1,y2]
annotation('textarrow',[.76,.68],[.75,.68],'String','250Hz','FontName', 'Times New Roman'); %[x1,x2],[y1,y2]

%% Put modes together
clear

load FRF

%% Determine valid chatter frequency range
FRF_real=real(FRF);
FRF_imag=imag(FRF);
index = find(FRF_real < 0);
FRF_real = FRF_real(index);
FRF_imag = FRF_imag(index);
w = w(index);

% Define force parameters
Kf=1000e6;          % N/m^2

%% Plot Lobe
% Calculate blim
blim = -1./(2*Kf*FRF_real);  % m
blim = blim*1e3;        % convert to mm

%% Calculate epsilon (1)
% for cnt = 1:length(FRF_imag)
%     if FRF_imag(cnt)<0
%         epsilon(cnt)=2*pi-2*atan(abs(FRF_real(cnt)/FRF_imag(cnt)));
%     else
%         epsilon(cnt)=pi-2*atan(abs(FRF_imag(cnt)/FRF_real(cnt)));
%     end
% end

%% %% Calculate epsilon (2)
for cnt = 1:length(FRF_imag)
    if FRF_imag(cnt)<0
        phi(cnt) = -pi+atan(abs(FRF_imag(cnt)/FRF_real(cnt)));
    else
        phi(cnt) = -pi-atan(abs(FRF_imag(cnt)/FRF_real(cnt)));
    end

    epsilon(cnt) = 3*pi+2*phi(cnt);
end


% Calculate spindle speeds for N = 0 to 2
omega0 = w./(0*2*pi + epsilon);   % rps
omega1 = w./(1*2*pi + epsilon);
omega2 = w./(2*2*pi + epsilon);   % rps
omega3 = w./(3*2*pi + epsilon);

figure(4)
plot(omega0*60, blim, 'g', omega1*60, blim, 'g', omega2*60, blim, 'g', omega3*60, blim, 'g','linewidth',2)
axis([2000 15000 0 60])
set(gca,'FontSize', 14)
xlabel('\fontsize{10}\fontname{Times New Roman} \Omega / rpm')
ylabel('\fontsize{10}\fontname{Times New Roman}\itb_{lim} \rm/ mm')
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
set(gcf,'unit','centimeters','position',[15 15 13.53 9.03],'color','white');
title('\fontsize{10}Lobe Diagram')
hold on
grid on


%% Different Phi
clear

load FRF

%% Determine valid chatter frequency range
FRF_real=real(FRF);
FRF_imag=imag(FRF);
index = find(FRF_real < 0);
FRF_real = FRF_real(index);
FRF_imag = FRF_imag(index);
w = w(index);

% Define force parameters
Kf=1000e6;          % N/m^2

%% Plot Lobe
% Calculate blim
blim = -1./(2*Kf*FRF_real);  % m
blim = blim*1e3;        % convert to mm

%% Calculate epsilon (1)
% for cnt = 1:length(FRF_imag)
%     if FRF_imag(cnt)<0
%         epsilon(cnt)=2*pi-2*atan(abs(FRF_real(cnt)/FRF_imag(cnt)));
%     else
%         epsilon(cnt)=pi-2*atan(abs(FRF_imag(cnt)/FRF_real(cnt)));
%     end
% end

%% %% Calculate epsilon (2)
for cnt = 1:length(FRF_imag)
    if FRF_imag(cnt)<0
        phi(cnt) = -pi+atan(abs(FRF_imag(cnt)/FRF_real(cnt)));
    else
        phi(cnt) = -pi-atan(abs(FRF_imag(cnt)/FRF_real(cnt)));
    end

    epsilon(cnt) = 3*pi+2*phi(cnt);
end


% Calculate spindle speeds for N = 0 to 2
omega0 = w./(0*2*pi + epsilon);   % rps
omega1 = w./(1*2*pi + epsilon);
omega2 = w./(2*2*pi + epsilon);   % rps
omega3 = w./(3*2*pi + epsilon);
% omega4 = w./(4*2*pi + epsilon);

figure(5)
plot(omega0*60, blim, 'g', omega1*60, blim, 'g', omega2*60, blim, 'g', omega3*60, blim, 'g','linewidth',2)
axis([2000 15000 0 60])
set(gca,'FontSize', 14)
xlabel('\fontsize{10}\fontname{Times New Roman} \Omega / rpm')
ylabel('\fontsize{10}\fontname{Times New Roman}\itb_{lim} \rm/ mm')
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
set(gcf,'unit','centimeters','position',[15 15 13.53 9.03],'color','white');
title('\fontsize{10}Lobe Diagram')
hold on
grid on