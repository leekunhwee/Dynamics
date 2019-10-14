%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------Copyright---------------%
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
% Updated semi-discretization method for periodic delay-differential equations with discrete delay
% Two DOF 
clc
% close all
clear 
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

% m_tx = 0.03993; % mass x (kg)
% m_ty = 0.03993; % mass y (kg)

%% Cutting Parameters
aD = 0.05;% radial depth of cut
up_or_down =-1; % 1: up-milling, -1: down-milling
if up_or_down == 1  % up-milling
    fist = 0;% start angle % exit angle
    fiex = acos(1-2*aD); 
elseif up_or_down == -1 % down-milling % start angle % exit angle
    fist = acos(2*aD-1); 
    fiex = pi;
end 

%% Simulation Parameters
step_speed = 100;% steps of spindle speed
step_depth = 20;% steps of depth of cut
depth_st = 0e-3; % starting depth of cut (m) 
depth_fi = 10e-3;  % final depth of cut (m)
speed_st = 5e3;% starting spindle speed (rpm) Int.
speed_fi = 25e3; % final spindle speed (rpm)

%% computational parameters 
k = 40;% number of discretization interval over one period
intk = 20;% number of numerical integration steps for Equation (37)
m = k;% since time delay = time period
wa = 1/2;% since time delay = time period
wb = 1/2;% since time delay = time period
D = zeros(2*m +4,2*m +4);% matrix D
d = ones(2*m +2, 1); 
d(1 : 4) = 0;
D = D+diag(d,-2); 
D(5, 1) = 1; 
D(6, 2) = 1; 
% numerical integration of specific cutting force coefficient according to
% Equation (40)-(43)
hxx = zeros(1,k); 
hxy = zeros(1,k); 
hyx = zeros(1,k);  
hyy = zeros(1,k); 
fi = zeros(1,intk); 
g = zeros(1,intk);
dtr = 2*pi/N/k; % Delta_Phi,if Phi_p = 2 Pai/N
%% numerical integration of specific cutting force coefficient 
for i = 1 : k 
    for j = 1 : N   % loop for tooth j 
        for h = 1 : intk % loop for numerical integration of hi
            fi(h) = (i-1)*dtr +(j-1)*2*pi/N + h*dtr/intk; 
            if (fi(h)>= fist)&&(fi(h)<= fiex) 
                g(h) = 1;% tooth is in the cut
            else
                g(h) = 0;% tooth is out of cut
            end
        end
        hxx(i) = hxx(i)+sum(g.*(Kt.* cos(fi)+Kn.* sin(fi)).* sin(fi))/intk; 
        hxy(i) = hxy(i)+sum(g.*(Kt.* cos(fi)+Kn.* sin(fi)).* cos(fi))/intk; 
        hyx(i) = hyx(i)+sum(g.*(-Kt* sin(fi)+Kn.* cos(fi)).* sin(fi))/intk; 
        hyy(i) = hyy(i)+sum(g.*(-Kt* sin(fi)+Kn.* cos(fi)).* cos(fi))/intk;
    end
end
% start of computation 
for x = 1 : step_speed+1   % loop for spindle speeds
    speed = speed_st +(x-1)*(speed_fi-speed_st)/step_speed;% spindle speed
    tau = 60/speed/N;% time delay
    dt = tau/(m);% time step
    for y = 1 : step_depth+1 % loop for depth of cuts
        w = depth_st +(y-1)*(depth_fi-depth_st)/step_depth; % depth of cut 
        % construct transition matrix Fi
        Fi = eye(2*m + 4, 2*m + 4); 
        for i = 1 : m 
            A = zeros(4, 4);    % matrix Ai
            A(1, 3) = 1; A(2, 4) = 1; 
            A(3, 1) =-w0x^2-hxx(i)*w*w0x^2/k_tx; 
            A(3, 2) =-hxy(i)*w*w0x^2/k_tx; 
            A(3, 3) =-2*zetax*w0x; 
            A(4, 1) =-hyx(i)*w*w0y^2/k_ty; 
            A(4, 2) =-w0y^2-hyy(i)*w*w0y^2/k_ty; 
            A(4, 4) =-2*zetay*w0y; 
            B = zeros(4, 4);% matrix Bi
            B(3, 1) = hxx(i)*w*w0x^2/k_tx; 
            B(3, 2) = hxy(i)*w*w0x^2/k_tx; 
            B(4, 1) = hyx(i)*w*w0y^2/k_ty; 
            B(4, 2) = hyy(i)*w*w0y^2/k_ty; 
            P = expm(A*dt);% matrix Pi
            R = (expm(A*dt)-eye(4))/A * B; % matrix Ri 
            D(1 : 4, 1 : 4) = P; 
            D(1 : 4,(2*m + 1) : (2*m + 2)) = wa*R(1 : 4, 1 : 2); 
            D(1 : 4,(2*m + 3) : (2*m + 4)) = wb*R(1 : 4, 1 : 2); 
            Fi = D*Fi;  % transition matrix Phi
        end
        ss(x, y) = speed;   % matrix of spindle speeds
        dc(x, y) = w*1000;   % matrix of depth of cuts
        ei(x, y) = max(abs(eig(Fi)));   % matrix of eigenvalues
    end
    step_speed+1-x
end
figure
contour(ss,dc,ei,[1, 1],'k')
xlabel('Spindle speed [rev/min]');
ylabel('a_l_i_m [mm]');
title('Stability of a two DOF milling process');
set(gcf,'unit','centimeters','position',[18 5 13.53 9.03],'color','white');%¶ÔÓ¦word£¨13.5,9£©
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
axis([speed_st speed_fi depth_st*1000 depth_fi*1000])
grid on