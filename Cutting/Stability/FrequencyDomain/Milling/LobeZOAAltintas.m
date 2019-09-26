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
%% Have some Problems: Nagetive frequency points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Zero-Order Solution            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all

load FRF

%% Cutting force coefficients
Kt = 1319.4e6; % N/m^2
Kn = 788.8/1319.4;
Nt = 2;

%% Cutting condition
D=31.75;       % mm
ae=D/2;
operation = 1; % down-milling

%% Cutting angle
phis = acos((D/2 - ae)/(D/2));  % (rad) the immersion angle

if operation == 0             % Up milling 
    phist = 0;                % Start angle rad
    phiex = phist + phis;     % Exit angle rad
else                          % Down milling
    phiex = pi;               % Start angle rad
    phist = phiex - phis;     % Exit angle rad
end

%% Average angle parameters
alphaxx = 0.5*(( cos(2*phiex)-2*Kn*phiex+Kn*sin(2*phiex))-( cos(2*phist)-2*Kn*phist+Kn*sin(2*phist)));
alphaxy = 0.5*((-sin(2*phiex)-2*phiex+Kn*cos(2*phiex))   -(-sin(2*phist)-2*phist+Kn*cos(2*phist)));
alphayx = 0.5*((-sin(2*phiex)+2*phiex+Kn*cos(2*phiex))   -(-sin(2*phist)+2*phist+Kn*cos(2*phist)));
alphayy = 0.5*((-cos(2*phiex)-2*Kn*phiex-Kn*sin(2*phiex))-(-cos(2*phist)-2*Kn*phist-Kn*sin(2*phist)));

%% Initialization
w = (f_start:df:f_end)'*2*pi;   % frequency, rad/s
FRFxx = FRFX;
FRFyy = FRFY;

%% Calculate Lambda
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reference 
%Altintas, Y. (2012). Manufacturing Automation: 
%Metal Cutting Mechanics, Machine Tool Vibrations, and CNC Design. 
%Applied Mechanics Reviews (Vol. 54). https://doi.org/10.1115/1.1399383
%P156
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a0 = FRFxx.*FRFyy*(alphaxx*alphayy - alphaxy*alphayx);
a1 = alphaxx*FRFxx + alphayy*FRFyy;

Lambda1 = zeros(length(w),1);
Lambda2 = zeros(length(w),1);

for cnt = 1:length(w)
    lambda	= roots([a0(cnt) a1(cnt) 1]);	% Eigenvalue 
    Lambda1(cnt) = lambda(1);
    Lambda2(cnt) = lambda(2);
end

blim1 = -(2*pi/Nt/Kt) .* (real(Lambda1) .* (1 + (imag(Lambda1)./real(Lambda1)).^2));  % m
blim2 = -(2*pi/Nt/Kt) .* (real(Lambda2) .* (1 + (imag(Lambda2)./real(Lambda2)).^2));

[index1] = find(blim1 > 0);
blim1 = blim1(index1);
blim1 = blim1 * 1e3;      % mm
w1 = w(index1);
psi1 = atan2(imag(Lambda1), real(Lambda1));
psi1 = psi1(index1);
epsilon1 = 3*pi - 2 * psi1;

[index2] = find(blim2 > 0);
blim2 = blim2(index2);
blim2 = blim2 * 1e3;
w2 = w(index2);
psi2 = atan2(imag(Lambda2), real(Lambda2));
psi2 = psi2(index2);
epsilon2 = 3*pi - 2 * psi2;

lobeNumber = 20;

omega1 = {lobeNumber,1};
for k = 1:lobeNumber
    omega1{k} = (60/Nt)*w1./(epsilon1 + 2*(k-1)*pi);
end

omega2 = {lobeNumber,1};
for k = 1:lobeNumber
    omega2{k} = (60/Nt)*w2./(epsilon2 + 2*(k-1)*pi);
end

figure(1)
hold on
for k = 1:lobeNumber
    plot(omega1{k},blim1, 'r','linewidth',2);
    plot(omega2{k},blim2, 'b','linewidth',2);
end

axis([0 16000 0 6])
xlabel('\it\Omega \rm/ rpm')
ylabel('\itb_{lim} \rm/ mm')
set(gca,'FontSize', 11 ,'FontName', 'Times New Roman')
set(gcf,'unit','centimeters','position',[28 5 13.53 9.03],'color','white');% word 13.5,9
set(gca,'xtick',[2000 4000 6000 8000 10000 12000 14000 16000])
grid on