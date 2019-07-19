%% Read Test FRF Data from .xls files
clc
clear
close all

% Start frequency 50Hz; End frequency 25000Hz; Resolution 0.390625Hz
f1 = 50;
f2 = 25000;
r = 0.390625;
f = f1:r:f2;
g = 9.8; % Gravitational Acceleration

% Mass of 352C23 acceleration sensor: 0.2g. Refer to PCB User Manual in Page 11
m_acc = 2e-4; % kg

%% Read FRF data from .xls files
% Notice!!! the Unit of acceleration FRF in File BEP0620.xls is g/N, 
% Multiply by g can transfer to (m/s^2)/N, Which should be double checked!

%% FRF of Beam B, Direct FRF
b_190_h11_real = g*(xlsread('C:\Work\RCSA\Example\BEP0620.xls', 'b_190_h11', 'B249:B64121'))';
b_190_h11_imag = g*(xlsread('C:\Work\RCSA\Example\BEP0620.xls', 'b_190_h11', 'C249:C64121'))';
b_190_h11 = b_190_h11_real + 1i * b_190_h11_imag; % Acceleration FRF 
b_190_h11n = b_190_h11./(1 - m_acc * b_190_h11);  % Eliminate the mass of transducer
b_190_h11NEW = b_190_h11n./(-(f * 2 * pi).^2);    % Transfer from acceleration to displacement
% Refer to: Cakar, O., & Sanliturk, K. Y. (2005). Elimination of transducer 
% mass loading effects from frequency response functions. 
% Mechanical Systems and Signal Processing, 19(1), 87?104. 
% https://doi.org/10.1016/s0888-3270(03)00086-4
% Aij* = Aij/(1-m*Aii) 

b_190_h12_real = g*(xlsread('C:\Work\RCSA\Example\BEP0620.xls', 'b_190_h21', 'B249:B64121'))';
b_190_h12_imag = g*(xlsread('C:\Work\RCSA\Example\BEP0620.xls', 'b_190_h21', 'C249:C64121'))';
b_190_h12 = b_190_h12_real + 1i*b_190_h12_imag;
b_190_h12n = b_190_h12./(1 - m_acc*b_190_h11);
b_190_h12NEW = b_190_h12n./(-(f * 2 * pi).^2);

%% FRF of Beam C 
c_140_h11_real = g*(xlsread('C:\Work\RCSA\Example\BEP0620.xls', 'c_140_h11', 'B249:B64121'))';
c_140_h11_imag = g*(xlsread('C:\Work\RCSA\Example\BEP0620.xls', 'c_140_h11', 'C249:C64121'))';
c_140_h11 = c_140_h11_real + 1i*c_140_h11_imag;
c_140_h11n = c_140_h11./(1 - m_acc*c_140_h11);
c_140_h11NEW = c_140_h11n./(-(f * 2 * pi).^2);

%% FRF of Beam D
d_190_h11_real = g*(xlsread('C:\Work\RCSA\Example\BEP0620.xls', 'd_190_h11', 'B249:B64121'))';
d_190_h11_imag = g*(xlsread('C:\Work\RCSA\Example\BEP0620.xls', 'd_190_h11', 'C249:C64121'))';
d_190_h11 = d_190_h11_real + 1i*d_190_h11_imag;
d_190_h11n = d_190_h11./(1 - m_acc*d_190_h11);
d_190_h11NEW = d_190_h11n./(-(f * 2 * pi).^2);

d_190_h12_real = g*(xlsread('C:\Work\RCSA\Example\BEP0620.xls', 'd_190_h21', 'B249:B64121'))';
d_190_h12_imag = g*(xlsread('C:\Work\RCSA\Example\BEP0620.xls', 'd_190_h21', 'C249:C64121'))';
d_190_h12 = d_190_h12_real + 1i*d_190_h12_imag;
d_190_h12n = d_190_h12./(1 - m_acc*d_190_h11);
d_190_h12NEW = d_190_h12n./(-(f * 2 * pi).^2);

% %% plot 
% figure(1)
% plot(f,real(b_190_h11NEW),'ro','linewidth',1.5)

%% Save .mat
save('BEPData','b_190_h11NEW','b_190_h12NEW','c_140_h11NEW','d_190_h11NEW','d_190_h12NEW')