%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------Copyright----------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Jianhui Li                         %
% Time: 03/14/2019                           %
% University of British Columbia, BC, Canada %
% Affiliation:                               %
% Department of Mechanical Engineering       %
% Manufacturing Automation Laboratary        %
% E-mail: jianhui.li@alumni.ubc.ca           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Semi-discretization Method 
% for periodic delay-differential equations 
% with discrete delay (One DOF)
clc
clear
% close all
%% Tool Parameters 
Nt = 2;% number of teeth
Kt = 6e8; % tangential cutting force coefficient (N/m2)
Kn = 2e8; % normal cutting force coefficient (N/m2)
%% Dynamic Parameters
w0 = 922*2*pi; % angular natural frequency (rad/s) 
zeta = 0.011;  % relative damping (1)
% m_t = 0.03993; % mass (kg)
k_t=1340049.648; %N/m
%% Cutting Parameters
aD = 1;% radial depth ratio of cut £» Slot
up_or_down =-1; % 1: up-milling, -1: down-milling
if up_or_down == 1  % up-milling
    fist = 0;% start angle % exit angle
    fiex = acos(1-2*aD); 
elseif up_or_down == -1 % down-milling  
    fist = acos(2*aD-1); % start angle
    fiex = pi;% exit angle
end 
%% Simulation Parameters
step_speed = 100;% steps of spindle speed
step_depth = 50;% steps of depth of cut
depth_st = 0e-3; % starting depth of cut (m) 
depth_fi = 10e-3;  % final depth of cut (m)
speed_st = 5e3;% starting spindle speed (rpm) Int.
speed_fi = 25e3; % final spind le speed (rpm)
%% computational parameters 
k = 40;% number of discretization interval over one period T
intk = 20;% number of numerical integration steps in ti to ti+1 for Equation (37)
m = k;% since time delay = time period
wa = 1/2;% since time delay = time period
wb = 1/2;% since time delay = time period
D = zeros(m +2,m +2);% matrix D
d = ones(m +1, 1); 
d(1 : 2) = 0;
D = D+diag(d,-1); 
D(3, 1) = 1; 
dfi = 2*pi/Nt/k; % ¡÷¦Õ,if ¦Õp = 2¦Ð/N step Pitch angle
h_i = zeros(1,k);
g = zeros(1,intk);
fi = zeros(1,intk);
%% numerical integration of specific cutting force coefficient 
% according to Equation (37)
for i = 1 : k 
%     h_i(i) = 0;
    for j = 1 : Nt % loop for tooth j
        for dint = 1 : intk % loop for numerical integration of hi
            fi(dint) = (i-1)*dfi +(j-1)*2*pi/Nt +dint*dfi/intk; 
            if (fi(dint)>= fist)&&(fi(dint)<= fiex) 
                g(dint) = 1;% tooth is in the cut
            else
                g(dint) = 0;% tooth is out of cut
            end
        end
        % Add all teeth
        h_i(i) = h_i(i)+sum(g.*(Kt.* cos(fi)+Kn.* sin(fi)).* sin(fi))/intk;
    end
end
%% start of computation
for x = 1 : step_speed+1% loop for spindle speeds
    speed = speed_st +(x-1)*(speed_fi-speed_st)/step_speed;% spindle speed
    tau = 60/speed/Nt;% time delay
    dt = tau/(m);% time step
    for y = 1 : step_depth+1 % loop for depth of cuts 
        w = depth_st +(y-1)*(depth_fi-depth_st)/step_depth; % depth of cut % construct transition matrix Fi
        Fi = eye(m +2,m +2); 
        for i = 1 : m % different from turning
            A = zeros(2, 2);% matrix Ai
            A(1, 2) = 1; 
            A(2, 1) =-w0^2-h_i(i)*w*w0^2/k_t; 
            A(2, 2) =-2*zeta*w0; 
            B = zeros(2, 2);% matrix Bi 
            B(2, 1) = h_i(i)*w*w0^2/k_t; 
            P = expm(A*dt);% matrix Pi
            R = (expm(A*dt)-eye(2))/A*B; % matrix Ri 
            D(1 : 2, 1 : 2) = P; 
            D(1 : 2,m + 1) = wa*R(1 : 2, 1 : 1); 
            D(1 : 2,m + 2) = wb*R(1 : 2, 1 : 1); 
            Fi = D*Fi;% transition matrix phi
        end
        ss(x, y) = speed;% matrix of spindle speeds
        dc(x, y) = w*1000;% matrix of depth of cuts
        ei(x, y) = max(abs(eig(Fi)));% matrix of eigenvalues
    end
step_speed+1-x          % Process Display
end
%% figure 
hold on
contour(ss,dc,ei,[1, 1])
xlabel('Spindle speed [rev/min]');
ylabel('a_l_i_m [mm]');
title('Stability of a one DOF milling process');
set(gcf,'unit','centimeters','position',[18 5 13.53 9.03],'color','white');
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
axis([speed_st speed_fi depth_st*1000 depth_fi*1000])
grid on