clc
clear
close all
%% flag
flag = 1; % Residues type

%% Parameters
Parameter = struct;
% Dynamics data
Parameter.Dynamic.f.x = [452.77 1448.53];
Parameter.Dynamic.f.y = [516.17 1407.64];
Parameter.Dynamic.zeta.x = [0.123718 0.01651];
Parameter.Dynamic.zeta.y = [0.0243    0.0324];
Parameter.Dynamic.residues.x = [9.202966E-05-1j*1.862195E-04 -4.181562E-05-1j*3.043618E-04];
Parameter.Dynamic.residues.y = [-2.3929E-06-1j*1.721539E-04   4.055052E-05-1j*3.618808E-04];
% Simulation parameter
Parameter.Simulation.f_start = 0;       % Hz
Parameter.Simulation.f_end = 8000;  % Hz
Parameter.Simulation.df = 1;            % Hz
% FRF Generation
FRF = FRFGeneration(Parameter.Dynamic,Parameter.Simulation,flag);
% Cutting Coefficient
Parameter.CuttingCoef.Kt = 1319.4e6; % N/m^2
Parameter.CuttingCoef.Kr = 788.8e6;
% Tool geometry
Parameter.ToolGeo.Nt = 2;
Parameter.ToolGeo.D = 31.75;
% Cutting parameter
Parameter.Cutting.ae = Parameter.ToolGeo.D/2;
Parameter.Cutting.operation = -1; % -1 Downmilling ; 1 Upmilling 

%% Plot Lobe
ZOA(FRF,Parameter.CuttingCoef,Parameter.ToolGeo,Parameter.Cutting)

%% Adjustment
axis([0 16000 0 6])
xlabel('\it\Omega \rm/ rpm')
ylabel('\itb_{lim} \rm/ mm')
set(gca,'FontSize', 11 ,'FontName', 'Times New Roman')
set(gcf,'unit','centimeters','position',[28 5 13.53 9.03],'color','white');% word£¨13.5,9£©
set(gca,'xtick',[2000 4000 6000 8000 10000 12000 14000 16000])
grid on