clc
clear
close all

load IRCSA
load BEPData
load Beam_Theoretical

% simulation parameters
f1 = 50;
f2 = 25000;
r = 0.390625;
f = f1:r:f2;
N = length(f);

% G11RCSA is the Full FRFs Matrix after coupling
G11RCSAFinal = RCSA(G11IRCSA, RD11, RD21, RD12, RD22, N);
H11RCSAFinal = reshape(G11RCSAFinal(1,1,:),1,N); % This is the Displacement/Force FRFs get by RCSA method 
L11RCSAFinal = reshape(G11RCSAFinal(1,2,:),1,N);
N11RCSAFinal = reshape(G11RCSAFinal(2,1,:),1,N);
P11RCSAFinal = reshape(G11RCSAFinal(2,2,:),1,N);

%% Plot
figure(1)
subplot(211)
plot(f,real(d_190_h11NEW),'r',f,real(H11RCSAFinal),'b','LineWidth',1.5);
subplot(212)
plot(f,imag(d_190_h11NEW),'r',f,imag(H11RCSAFinal),'b','LineWidth',1.5);