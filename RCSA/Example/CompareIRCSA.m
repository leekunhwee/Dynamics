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
L = 0.05; % The displacement between exciting and measuring points of cross FRFs 

% G11IRCSA is the Full FRFs Matrix after removing
G11IRCSA = Experimental_IRCSA(b_190_h11NEW, b_190_h12NEW, L, RB11, RB21, RB12, RB22);
H11IRCSA = reshape(G11IRCSA(1,1,:),1,N); % This is the Displacement/Force FRFs get by Inverse RCSA method 
L11IRCSA = reshape(G11IRCSA(1,2,:),1,N);
N11IRCSA = reshape(G11IRCSA(2,1,:),1,N);
P11IRCSA = reshape(G11IRCSA(2,2,:),1,N);

figure(1)
subplot(211)
plot(f,real(c_140_h11NEW),'r',f,real(H11IRCSA),'b','LineWidth',1.5);
subplot(212)
plot(f,imag(c_140_h11NEW),'r',f,imag(H11IRCSA),'b','LineWidth',1.5);

%% Save the FRFs at the end of a modified (removed some parts) cantilever beam calculated by Inverse RCSA method
save('IRCSA', 'G11IRCSA', 'H11IRCSA', 'L11IRCSA', 'N11IRCSA', 'P11IRCSA')