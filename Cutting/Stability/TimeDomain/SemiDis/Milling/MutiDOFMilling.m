clc
clear
close all

%% Cutting Condition
N  = 2;        % Number of teeth
D  = 31.75;    % Diameter of tool mm
ae = D/2;      % Radial immersion depth
operation = 1; % Down-milling

%% Cutting angle
phis = acos((D/2 - ae)/(D/2));  % (rad) the immersion angle

if operation == 0               % Up milling 
    phist = 0;                  % Start angle rad
    phiex = phist + phis;       % Exit angle rad
else                            % Down milling
    phiex = pi;                 % Start angle rad
    phist = phiex - phis;       % Exit angle rad
end

%% Material parameter
Ktc = 750e6; %4042e6;%758e6;
Krc = 150e6; %2814e6;%367e6; %MPa
Kr  = Krc / Ktc; 

%% Modal parameters
fx = [369.85*2*pi,1146.36*2*pi];%1268.1*2*pi;
% kx = 1.0481e7;
% mx = [1.428,0.163];
mx = [1,1];
zeta_x=[0.0316,0.0114];%0.022;
%modal matrix
Px = [0.8369,2.4757];

fy = [411.21*2*pi, 1069.21*2*pi];%1289*2*pi;
% ky=1.3718e7;%58062000;
% my=[1.057,0.129];
my = [1,1];
zeta_y = [1.68*0.01,1.55*0.01];%0.014;
%modal matrix Y
Py = [0.9727,2.7837];

Cq = diag([zeta_x*2.*fx,zeta_y*2.*fy]);
Kq = diag([fx.*fx      ,      fy.*fy]);
mq = diag([mx          ,          my]);
invMq = inv(mq);

UU=[
    Px,0,0;
    0,0,Py
    ];

%% Simulation parameters
m = 5; % discrete number

rpm_array     = 500:20:5000;
a_array       = 0.001*(0:1:20);
spd_mat       = zeros(length(rpm_array),length(a_array));
depth_cut_mat = zeros(size(spd_mat));
eig_mat       = zeros(size(spd_mat));

[a_xx,a_xy,a_yx,a_yy] = genDynamicCoe1(m,N,phist,phiex,Kr);

a_min = 0;
a_max = 10e-3;
delta_a = .1e-4;

aa = zeros(size(rpm_array));
    
%count time
tic
t1=toc;
% for ii=1:length(rpm_array)
%     rpm=rpm_array(ii); %speed of the spindle RPM
%     w=rpm/60*2*pi; %speed of the spindle in rad/s
%     w_T=w*N;
%     T=2*pi/w_T; % tooth period
%     dt=T/m;
%     for jj=1:length(a_array)
%         %Depth of cut
%         a=a_array(jj);%3/1000; 
%         Delta=a*Ktc/2;
% %         [a_xx,a_xy,a_yx,a_yy]=genDynamicCoe(dt,T,N,phi_st,phi_ex,w,Kr);
%         eigMax=checkStabilityMultiMode_2(Cq,Kq,invMq,UU,Delta,a_xx,a_xy,a_yx,a_yy,m,dt);
%         spd_mat(ii,jj)=rpm;
%         depth_cut_mat(ii,jj)=a;
%         eig_mat(ii,jj)=eigMax;
%     end
% end


for ii=1:length(rpm_array)
    rpm=rpm_array(ii); %speed of the spindle RPM
    w=rpm/60*2*pi; %speed of the spindle in rad/s
    w_T=w*N;
    T=2*pi/w_T; % tooth period
    dt=T/m;
    p=a_max-a_min;
    a1=a_max;
    a2=a_min;
    last_a=a_min;
    while (p>delta_a)
        a=(a1+a2)/2;
        Delta=a*Ktc/2;
%         eigMax=checkStability(fx,kx,fy,ky,Delta,zeta_x,zeta_y,a_xx,a_xy,a_yx,a_yy,m,dt);
        eigMax = checkStabilityMultiMode_2(Cq,Kq,invMq,UU,Delta,a_xx,a_xy,a_yx,a_yy,m,dt);
        if eigMax<1
            p = abs(a-last_a);
            last_a = a;      
            a2 = a;
        else
            a1 = a;
        end
    end
    aa(ii) = a;
end

t2=toc;
running_t = t2-t1

% figure
% contour(spd_mat,depth_cut_mat,eig_mat,[1, 1],'LineWidth',2);
% title('up milling, half immersion, 1045 steel')
figure;
plot(rpm_array,aa,'LineWidth',2);

data=[rpm_array;aa]';
save('al_two_modes_2T_down','data');

% contour(spd_mat,depth_cut_mat,eig_mat,[1, 1],'LineWidth',2);
% title('up milling, half immersion, 1045 steel')