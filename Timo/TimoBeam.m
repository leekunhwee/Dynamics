% Timoshenko Beam Tip FRF 
% Verified ¡Ì
% Jianhui.Li
% 2018.12.05
clc
clear
close all

% Simulation parameters
w = (1:5000)*2*pi;          % frequency, rad/s
f = w/2/pi;

% Define free-free cylinder parameters

E = 200e9;                  % elastic modulus, N/m^2
d = 10e-3;                   % diameter, m
L = 125e-3;                 % length, m
I = pi*d^4/64;              % 2nd moment of area, m^4
rho = 7800;                 % density, kg/m^3
A = pi*d^2/4;               % cross sectional area, m^2
eta = 0.01;                 % solid damping factor
EI = E*I*(1+1i*eta);        % complex stiffness, N-m^2
el = 1e-3;                  % element length, m
n = ceil(L/el);             % number of element
nu = 0.29;                          % Poisson's ratio
G = E/(2*(1+nu));                   % shear modulus, N/m^2
AG = G*A;                           % produce of area and shear modulus for composite beam, N
rg = (I/A).^0.5;                    % radius of gyration, m
mpl = rho*A;                        % mass per unit length, kg/m
kp = 6*(1+nu)^2/(7+12*nu+4*nu^2);   % shape factor

% free-free two end receptances
[Rt11, Rt21, Rt12, Rt22] = timo_free_free(f, EI, L, AG, kp, rg, mpl, n);

% FRF of wall
Rtw=zeros(2,2,length(w));

% Assemble tip receptance of the cantilever beam
Gt11 = zeros(2,2,length(w));

for cnt = 1:length(w)
    Gt11(:,:,cnt) = Rt11(:,:,cnt) - Rt12(:,:,cnt)*((Rt22(:,:,cnt) + Rtw(:,:,cnt))\Rt21(:,:,cnt));
end

% Assemble tip FRF of the cantilever beam
Ht11=reshape(Gt11(1,1,:),1,length(w));

% Plot tip FRF of the cantilever beam
figure(1)
subplot(211)
plot(w/2/pi,real(Ht11),'k')
% axis([300 700 -3.75e-4 3.75e-4])
set(gca,'FontSize', 14)
ylabel('Real (m/N)')
subplot(212)
plot(w/2/pi,imag(Ht11),'k')
% axis([300 700 -7e-4 7e-5])
set(gca,'FontSize', 14)
xlabel('Frequency (Hz)')
ylabel('Imag (m/N)')

