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
% for periodic delay-differential equations with discrete delay (Two DOF)
clc
clear
close all

tic
%% Tool Parameters 
N = 2;% number of teeth
Kt = 6e8; % tangential cutting force coefficient (N/m2)
Kn = 2e8; % normal cutting force coefficient (N/m2)

%% Dynamic Parameters
w0x = 922*2*pi; % angular natural frequency x (rad/s) 
zetax = 0.011;  % relative damping x(1)
w0y = 922*2*pi; % angular natural frequency y(rad/s) 
zetay = 0.011;  % relative damping y(1)

k_tx = 1340049.648; %N/m
k_ty = 1340049.648; %N/m
%% Cutting Parameters
aD = 0.05;% radial depth of cut
up_or_down = 1; % 1: up-milling, -1: down-milling
if up_or_down == 1  % up-milling
    fist = 0;% start angle % exit angle
    fiex = acos(1-2*aD); 
elseif up_or_down == -1 % down-milling  
    fist = acos(2*aD-1); % start angle
    fiex = pi;% exit angle
end 

%% Simulation Parameters
step_speed = 200;% steps of spindle speed
step_depth = 40;% steps of depth of cut
depth_st = 0e-3; % starting depth of cut (m) 
depth_fi = 10e-3;  % final depth of cut (m)
speed_st = 5e3;% starting spindle speed (rpm) Int.
speed_fi = 25e3; % final spind le speed (rpm)

%% computational parameters 
m = 40; % number of discretization interval over one period T
D = zeros(2*m +4,2*m +4);% matrix D
d = ones(2*m +2, 1); 
d(1 : 4) = 0;
D = D+diag(d,-2); 
D(5, 1) = 1; 
D(6, 2) = 1; 

hxx = zeros(1,m); 
hxy = zeros(1,m); 
hyx = zeros(1,m);  
hyy = zeros(1,m); 

%% numerical integration of specific cutting force coefficient 
for i = 1 : m + 1 
    dfi = 2*pi/N/m; % Delta_Phi,if Phi_p = 2*Pai/N
    hxx(i) = 0; 
    hxy(i) = 0; 
    hyx(i) = 0; 
    hyy(i) = 0; 
    for j = 1 : N % loop for tooth j
        fi = i*dfi + (j-1)*2*pi/N;
        if (fi >= fist)&&(fi <= fiex) 
            g = 1;
        else
            g = 0;
        end
        hxx(i) = hxx(i)+g*(Kt* cos(fi)+Kn* sin(fi))* sin(fi); 
        hxy(i) = hxy(i)+g*(Kt* cos(fi)+Kn* sin(fi))* cos(fi); 
        hyx(i) = hyx(i)+g*(-Kt* sin(fi)+Kn* cos(fi))* sin(fi); 
        hyy(i) = hyy(i)+g*(-Kt* sin(fi)+Kn* cos(fi))* cos(fi);
    end
end

% Begin of the proposed method
A0 = [0     0       w0x^2/k_tx      0;
      0     0       0               w0y^2/k_ty;
    -k_tx   0       -2*zetax*w0x    0;
      0     -k_ty   0               -2*zetay*w0y];
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
        Fi = eye(2*m + 4,2*m +4); 
        for i = 1 : m 
            A0k = [ 0,              0,              0,              0;
                    0,              0,              0,              0;
                    -w*hxx(i+1)     -w*hxy(i+1)     0               0;
                    -w*hyx(i+1)     -w*hyy(i+1)     0               0];
                
            A1k = [ 0,                      0,                          0,              0;
                    0,                      0,                          0,              0;
                    w*(hxx(i+1)-hxx(i))/dt  w*(hxy(i+1)-hxy(i))/dt      0,              0;
                    w*(hyx(i+1)-hyx(i))/dt  w*(hyy(i+1)-hyy(i))/dt      0,              0];
                
            B0k = [ 0,              0,              0,              0;
                    0,              0,              0,              0;
                    w*hxx(i+1)      w*hxy(i+1)      0               0;
                    w*hyx(i+1)      w*hyy(i+1)      0               0];
            
            B1k = [ 0,                      0,                          0,              0;
                    0,                      0,                          0,              0;
                    w*(hxx(i)-hxx(i+1))/dt  w*(hxy(i)-hxy(i+1))/dt      0,              0;
                    w*(hyx(i)-hyx(i+1))/dt  w*(hyy(i)-hyy(i+1))/dt      0,              0];
            
            F01 = Fi2*A0k/dt + Fi3*A1k/dt;              % F0,1
            Fkp1 = (Fi1-Fi2/dt)*A0k + (Fi2-Fi3/dt)*A1k; % Fk+1
            Fm1 = (Fi1-Fi2/dt)*B0k + (Fi2-Fi3/dt)*B1k; % Fm-1
            Fm = Fi2*B0k/dt + Fi3*B1k/dt;              % Fm
            
            inv0fImFkp1 = inv(I-Fkp1);                  % [I-Fk+1]^-1  
            
            D(1:4,1:4) = inv0fImFkp1 * (Fi0+F01);
            D(1:4,2*m+1:2*m+2) = inv0fImFkp1 * Fm1(1:4,1:2);
            D(1:4,2*m+3:2*m+4) = inv0fImFkp1 * Fm(1:4,1:2); % A()and B() are negative of each other
            
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