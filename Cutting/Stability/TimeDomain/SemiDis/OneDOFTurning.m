%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------Copyright------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Jianhui Li                         %
% Time: 03/14/2019                           %
% University of British Columbia, BC, Canada %
% Affiliation:                               %
% Department of Mechanical Engineering       %
% Manufacturing Automation Laboratary        %
% E-mail: jianhui.li@alumni.ubc.ca           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Turning stability using semi-discrete time domain solution

clc
clear
close all

%% Dynamic Parameters
% One DOF
zeta = 0.012;               % Damping ratio
k = 2.26e8;                 % Stiffness N/m 
theta = 30;                 % Orientation Deg
k_y = k/(cosd(theta)^2);    % Orientation Stiffness N/m 
wn = 250*2*pi;              % Natural Frequency rad/sec
K_f = 1e9;                  % Cutting Constant N/m2

%% Simulation Parameters
deltaT = 0.1e-3;    % sec
depth_st = 0e-3;    % Starting depth of cut(m)
depth_fi = 60e-3;   % Final depth of cut(m)
speed_st = 2e3;     % Starting spindle speed(rpm)
speed_fi = 14e3;    % Final spindle speed(rpm)

%% Initial conditions
I = eye(2);
unstableA=0;
unstableN=0;

% Spindle Speed loop
for n=speed_st:100:speed_fi
    % Cutting depth loop
    da = 0.2e-3; % Axis depth step(m)
    for a=depth_st:da:depth_fi
        % Build the L and R matrix
        L = [0, 1; -(wn^2)*(1+(K_f*a/k_y)), -2*zeta*wn];
        R = [0, 0; (wn^2)*K_f*a/k_y, 0];
        % Spindle period
        T = 60/n;
        % Discrete time number
        m = round(T/deltaT);
        % Initial the B1 and B2 matrix, Attention! Column and row
        B1 = zeros(2*(m+1),2*(m+1));
        B2 = zeros(2*(m+1),2*(m+1));
        % Define the first two Column and Row of B1
        B1(1:2,1:2) = expm(L*deltaT);
        % Build the complete matrix of B1
        for i=1:2:2*m-1
            B1(i+2:i+3,i:i+1) = I;
        end      
        % Build the complete matrix of B2 
        B2((1:2),(end-1:end))   = 0.5*(expm(L*deltaT)-I)/L*R;
        B2((1:2),(end-3:end-2)) = 0.5*(expm(L*deltaT)-I)/L*R;
        % Build the B matrix
        B = B1+B2;
        % Check all the eigenvalues of B matrix
        if any(abs(eig(B))>=1)
            unstableA = [unstableA, a-da];
            unstableN = [unstableN, n];
%             plot(n,a,'. blue');
            disp([n a]);
            break;
        end
    end
end

% hold on
figure;

plot(unstableN,unstableA*1000,'.');
xlabel('Spindle speed [rev/min]');
ylabel('a_l_i_m [mm]');
title(['Stability of a 1-DOF shaping process - dt = ',num2str(deltaT),'s']);
set(gcf,'unit','centimeters','position',[18 5 13.53 9.03],'color','white');%∂‘”¶word£®13.5,9£©
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman','xtick',[2000 4000 6000 8000 10000 12000 14000 16000])
axis([2000 15000 0 60])
grid on