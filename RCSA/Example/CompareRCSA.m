clc
clear
close all

load BEPData
load Beam_Theoretical

% simulation parameters
f1 = 50;
f2 = 25000;
r = 0.390625;
f = f1:r:f2;
N = length(f);

% G11RCSA is the Full FRFs Matrix after coupling
G11RCSA = RCSA(RC, RB11, RB21, RB12, RB22, N);
H11RCSA = reshape(G11RCSA(1,1,:),1,N); % This is the Displacement/Force FRFs get by RCSA method 
L11RCSA = reshape(G11RCSA(1,2,:),1,N);
N11RCSA = reshape(G11RCSA(2,1,:),1,N);
P11RCSA = reshape(G11RCSA(2,2,:),1,N);

%% Plot
figure(1)
subplot(211)
plot(f,real(b_190_h11NEW),'r',f,real(H11RCSA),'b','LineWidth',1.5);
subplot(212)
plot(f,imag(b_190_h11NEW),'r',f,imag(H11RCSA),'b','LineWidth',1.5);