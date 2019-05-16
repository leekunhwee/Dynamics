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

%% Full-discretization Method 
% for periodic delay-differential equations with discrete delay (One DOF)
clc
clear
close all

tic
%% Tool Parameters 
N = 2;% number of teeth
Kt = 6e8; % tangential cutting force coefficient (N/m2)
Kn = 2e8; % normal cutting force coefficient (N/m2)

%% Dynamic Parameters
w0 = 922*2*pi; % angular natural frequency (rad/s) 
zeta = 0.011;  % relative damping (1)
% m_t = 0.03993; % mass (kg)
k_t=1340049.648; %N/m
%% Cutting Parameters
aD = 1;% radial depth of cut

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
step_depth = 20;% steps of depth of cut
depth_st = 0e-3; % starting depth of cut (m) 
depth_fi = 10e-3;  % final depth of cut (m)
speed_st = 5e3;% starting spindle speed (rpm) Int.
speed_fi = 25e3; % final spind le speed (rpm)

%% computational parameters 
% k = 40;% number of discretization interval over one period T
% intk = 20;% number of numerical integration steps for Equation (37)
% wa = 1/2;% since time delay = time period
% wb = 1/2;% since time delay = time period
m = 40; % number of discretization interval over one period T
D = zeros(m + 2,m +2);% matrix D
d = ones(m +1, 1); 
d(1 : 2) = 0;
D = D+diag(d,-1); 
D(3, 1) = 1; 

%% numerical integration of specific cutting force coefficient 
% according to Equation (29)
for i = 1 : m + 1 
    dfi = 2*pi/N/m; % Delta_Phi,if Phi_p = 2*Pai/N
    h_i(i) = 0; 
    for j = 1 : N % loop for tooth j
%         for dint = 1 : intk % loop for numerical integration of hi
%             fi(dint) = i*dfi +(j-1)*2*pi/N +dint*dfi/intk; 
        fi = i*dfi + (j-1)*2*pi/N;
%         if (fi(dint)>= fist)&&(fi(dint)<= fiex) 
%             g(dint) = 1;% tooth is in the cut
        if (fi >= fist)&&(fi <= fiex) 
            g = 1;
        else
            g = 0;
%             g(dint) = 0;% tooth is out of cut
%             end
        end
        h_i(i) = h_i(i)+g*(Kt*cos(fi)+Kn*sin(fi))*sin(fi);
%         h_i(i) = h_i(i)+sum(g.*(Kt.* cos(fi)+Kn.* sin(fi)).* sin(fi))/intk;
    end
end

% Begin of the proposed method
% A0 = [-zeta*w0,1/m_t;m_t*((zeta*w0)^2-w0^2),-zeta*w0];

A0 = [0,w0^2/k_t;-k_t,-2*zeta*w0];
I = eye(size(A0));
invA0=inv(A0);
%% start of computation
for x = 1 : step_speed+1% loop for spindle speeds
    speed = speed_st +(x-1)*(speed_fi-speed_st)/step_speed;% spindle speed
    tau = 60/speed/N;% time delay
    dt = tau/(m);% time step
    
    %--------------
    Fi0 = expm(A0*dt);
    Fi1 = invA0 * (Fi0-I);
    Fi2 = invA0 * (Fi0*dt-Fi1);
    Fi3 = invA0 * (Fi0*dt*dt - 2*Fi2);
    
    %--------------
    
    for y = 1 : step_depth+1 % loop for depth of cuts 
        w = depth_st +(y-1)*(depth_fi-depth_st)/step_depth; % depth of cut % construct transition matrix Fi
        Fi = eye(m + 2,m +2); 
        for i = 1 : m 
            A0k = [0,0;-w*h_i(i+1) 0 ];
            A1k = [0,0;w*(h_i(i+1)-h_i(i))/dt 0 ];
            B0k = [0,0;w*h_i(i+1) 0 ];
            B1k = [0,0;w*(h_i(i)-h_i(i+1))/dt 0 ];
            F01 = Fi2*A0k/dt + Fi3*A1k/dt;              % F0,1
            Fkp1 = (Fi1-Fi2/dt)*A0k + (Fi2-Fi3/dt)*A1k; % Fk+1
            Fm1 = (Fi1-Fi2/dt)*B0k + (Fi2-Fi3/dt)*B1k; % Fm-1
            Fm = Fi2*B0k/dt + Fi3*B1k/dt;              % Fm
            
            inv0fImFkp1 = inv(I-Fkp1);                  % [I-Fk+1]^-1  
            
            D(1:2,1:2) = inv0fImFkp1 * (Fi0+F01);
            D(1:2,m+1) = inv0fImFkp1 * Fm1(1:2,1:1);
            D(1:2,m+2) = inv0fImFkp1 * Fm(1:2,1:1); % A()and B() are negative of each other
            
            Fi = D*Fi;
        end
        ss(x, y) = speed;% matrix of spindle speeds
        dc(x, y) = w*1000;% matrix of depth of cuts
        ei(x, y) = max(abs(eig(Fi)));% matrix of eigenvalues
    end
    step_speed+1-x          % Process Display
end
toc
%%
figure 
contour(ss,dc,ei,[1, 1],'k')

xlabel('Spindle speed [rev/min]');
ylabel('a_l_i_m [mm]');
title('Stability of a one DOF milling process');
set(gcf,'unit','centimeters','position',[18 5 13.53 9.03],'color','white');
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
axis([speed_st speed_fi depth_st*1000 depth_fi*1000])
grid on