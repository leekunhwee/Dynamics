%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %---------------Copyright------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Jianhui Li                         %
% Time: 07/18/2019                           %
% University of British Columbia, BC, Canada %
% Affiliation:                               %
% Department of Mechanical Engineering       %
% Manufacturing Automation Laboratary        %
% E-mail: jianhui.li@alumni.ubc.ca           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculating FRFs using Timoshenko Beam Theory 
clc
clear
close all

%% Parameters setting
% simulation
f1 = 50;
f2 = 25000;
r = 0.390625;
f = f1:r:f2;
N = length(f);

%  Parameters of cantilever beam
l_1 = 0.19;     % Length of Beam: 200mm; Diameter: 16mm; With clamp length 10mm
d_1 = 0.016;    %
l_2 = 0.05;     % Length to be removed: 50mm; Diameter: 16mm
d_2 = 0.016;
l_3 = 0.14;     % Length of another Beam: 150mm; Diameter: 16mm; With clamp length 10mm
d_3 = 0.016;    %
l_4 = 0.05;     % Length to be Added: 50mm; Diameter: 10mm
d_4 = 0.010;    

% Material parameters
rho = 2.81e3;   % Density of 7075 Alluminum alloy
E = 7.2e10;     % Young's modulus
nu_v = 0.33;    % Poisson's ratio
el = 1e-2;      % m; Finite Element Length

%% FRF calculation 
% Calculating the Free-Free FRF of the beam with Length: 190mm; and  Diameter: 16mm
[RA11, RA21, RA12, RA22] = Beam_FRF(l_1, d_1,0, f1, f2, r, rho,0, E,0, nu_v,0, el); %
RA = zeros(2,2,N);
% Calculate the Direct FRF of the end of a Cantilever Beam with Length: 190mm; and  Diameter: 16mm
% Using RCSA method 
% Refer to: Schmitz, T. L., & Smith, K. S. (2009). Machining Dynamics. 
% (Intergovernmental Panel on Climate Change, Ed.), springer. Boston, MA: Springer US. 
% https://doi.org/10.1007/978-0-387-09645-2

for cnt = 1:N
RA(:,:,cnt) = RA11(:,:,cnt) - RA12(:,:,cnt)*((RA22(:,:,cnt)) \ RA21(:,:,cnt));
end

% Calculating the Free-Free FRF of the beam with Length: 50mm; and  Diameter: 16mm
[RB11, RB21, RB12, RB22] = Beam_FRF(l_2, d_2,0, f1, f2, r, rho,0, E,0, nu_v,0, el);%

% Calculating the Free-Free FRF of the beam with Length: 140mm; and  Diameter: 16mm
[RC11, RC21, RC12, RC22] = Beam_FRF(l_3, d_3,0, f1, f2, r, rho,0, E,0, nu_v,0, el);%
RC=zeros(2,2,N);
% Calculate the Direct FRF of the end of a Cantilever Beam with Length: 140mm; and  Diameter: 16mm
for cnt = 1:N
RC(:,:,cnt) = RC11(:,:,cnt) - RC12(:,:,cnt)*((RC22(:,:,cnt)) \ RC21(:,:,cnt));
end

% Calculating the Free-Free FRF of the beam with Length: 50mm; and  Diameter: 10mm
[RD11, RD21, RD12, RD22] = Beam_FRF(l_4, d_4,0, f1, f2, r, rho,0, E,0, nu_v,0, el);%

%% Save
% save the Theoretical results 
save('Beam_Theoretical','RA','RB11', 'RB21', 'RB12', 'RB22','RC','RD11','RD21','RD12','RD22');
% RA is the FRFs at the end of a Cantilever Beam with length 190mm
% RC is the FRFs at the end of a Cantilever Beam with length 140mm
% RB11 RB21 RB12 RB22 are the FRFs of a Free-Free Beam with length 50mm; Diameter = 16mm 
% l_2 = 50mm is the length to be removed
% RD11 RD21 RD12 RD22 are the FRFs of a Free-Free Beam with length 50mm; Diameter = 10mm 
% l_4 = 50mm is the length to be added