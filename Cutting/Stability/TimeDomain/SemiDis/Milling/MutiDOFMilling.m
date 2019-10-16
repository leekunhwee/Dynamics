%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------Copyright---------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Jianhui Li                         %
% Time: 10/15/2019                           %
% University of British Columbia, BC, Canada %
% Affiliation:                               %
% Department of Mechanical Engineering       %
% Manufacturing Automation Laboratary        %
% E-mail: jianhui.li@alumni.ubc.ca           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;
% close all
%% Cutting Parameters
N  = 4;         % Number of teeth
aD = 0.5;       % radial depth of cut / Ratio
up_or_down = -1; % 1: up-milling, -1: down-milling
if up_or_down == 1  % up-milling
    fist = 0;   % start angle % exit angle
    fiex = acos(1-2*aD);
elseif up_or_down == -1 % down-milling % start angle % exit angle
    fist = acos(2*aD-1); 
    fiex = pi;
end
%% Material parameter
Kt = 750e6; 
Kn = 150e6; % MPa
%% Modal parameters
% X direction
fx = [369.85*2*pi, 1146.36*2*pi]; % Natural Frequency
% kx = 1.0481e7;
% mx = [1.428,0.163];
mx = [1, 1];
zeta_x = [0.0316, 0.0114];  %0.022;  % Damping ratio
% Modal matrix X
Px = [0.8369, 2.4757];      % mode shape matrix
% Y direction
fy = [411.21*2*pi, 1069.21*2*pi];%1289*2*pi;% Natural Frequency
% ky=1.3718e7;%58062000;
% my=[1.057,0.129];
my = [1,1];
zeta_y = [1.68*0.01,1.55*0.01];%0.014;
% Modal matrix Y
Py = [0.9727, 2.7837];
% Build Dynamic Matrix
Cq = diag([zeta_x*2.*fx,zeta_y*2.*fy]);
Kq = diag([fx.*fx      ,      fy.*fy]);
mq = diag([mx          ,          my]);
invMq = inv(mq);
UU=[Px,0,0;
    0,0,Py];
M = length(fx)+length(fy);
%% computational parameters 
k = 40;% number of discretization interval over one period
intk =10;% number of numerical integration steps for Equation (37)
m = k;% since time delay = time period
wa = 1/2;% since time delay = time period
wb = 1/2;% since time delay = time period
D = zeros(M*m +2*M,M*m +2*M);% matrix D
d = ones(M*m +M, 1); 
d(1 : 2*M) = 0;
D = D + diag(d,-M); 
D((2*M+1):(2*M+M), 1:M)= eye(M); 
%% Simulation Parameters
step_speed = 100;% steps of spindle speed
step_depth = 50;% steps of depth of cut
depth_st = 0e-3; % starting depth of cut (m) 
depth_fi = 10e-3;  % final depth of cut (m)
speed_st = 5e2;% starting spindle speed (rpm) Int.
speed_fi = 5e3; % final spindle speed (rpm)
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
        Fi = eye(M*m + 2*M, M*m + 2*M); 
        for i = 1 : m 
            KK = w * UU' * [hxx(i),hxy(i);hyx(i),hyy(i)] * UU;
            L = zeros(2*M, 2*M);      % matrix Ri
            L(1:M, (M+1):2*M) = eye(M); 
            L((M+1):2*M, 1:M) = -Kq-KK;
            L((M+1):2*M, (M+1):2*M)=-Cq;
            R = zeros(2*M, 2*M);      % matrix Li
            R((M+1):2*M, 1:M) = KK;
            P = expm(L*dt);% matrix Pi
            Q = (expm(L*dt)-eye(2*M))/L * R; % matrix Ri 
            D(1 : 2*M, 1 : 2*M) = P; 
            D(1 : 2*M,(M*m + 1) : (M*m + M)) = wa*Q(1 : 2*M, 1 : M); 
            D(1 : 2*M,(M*m + M+1) : (M*m + 2*M)) = wb*Q(1 : 2*M, 1 : M); 
            Fi = D*Fi;  % transition matrix Phi
        end
        ss(x, y) = speed;   % matrix of spindle speeds
        dc(x, y) = w * 1000;   % matrix of depth of cuts
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
% axis([speed_st speed_fi depth_st*1000 depth_fi*1000])
axis([500 5000 0.8 2.4])

grid on