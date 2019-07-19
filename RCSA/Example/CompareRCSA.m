%% Compare the coupled 140mm 16mm + 50mm 16mm Beam with Experimental Results
clc
clear
close all

load BEPData
load Beam_Theoretical

% simulation parameters
f1 = 50;        % Frequency start
f2 = 25000;     % Frequency end
r = 0.390625;   % Frequency resolution
f = f1:r:f2;    % Frequency point
fplot1 = 50;    % Plot start
fplot2 = 3000;  % Plot end
N = length(f);  % Number of frequency points

% G11RCSA is the Full FRFs Matrix after coupling
% RC is the end FRFs of the 140_16mm cantilever beam 
G11RCSA = RCSA(RC, RB11, RB21, RB12, RB22, N);
H11RCSA = reshape(G11RCSA(1,1,:),1,N); % This is the Displacement/Force FRFs get by RCSA method 
L11RCSA = reshape(G11RCSA(1,2,:),1,N); % This is the Displacement/Torque FRFs
N11RCSA = reshape(G11RCSA(2,1,:),1,N); % This is the Angle/Force FRFs
P11RCSA = reshape(G11RCSA(2,2,:),1,N); % This is the Angle/Torque FRFs

%% Plot 
figure(1)

subplot(211)
plot(f,real(b_190_h11NEW),'r',f,real(H11RCSA),'b','LineWidth',1.5);
xlim([fplot1,fplot2])
ylabel('\it H_{11} \rm m / N');
title('End FRF of 190mm \phi 16mm Beam');
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
legend('Experimental', 'Experimental+RCSA')
grid minor

subplot(212)
plot(f,imag(b_190_h11NEW),'r',f,imag(H11RCSA),'b','LineWidth',1.5);
xlim([fplot1,fplot2])
ylabel('\it H_{11} \rm m / N');
xlabel('\it Frequency \rm Hz');
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
legend('Experimental', 'Experimental+RCSA')
set(gcf,'unit','centimeters','position',[18 5 13.53 9.03],'color','white'); % word 13.5 9
grid minor