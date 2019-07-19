%% Compare the removed 190mm 16mm - 50mm 16mm Beam with Experimental Results
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
L = 0.05; % The displacement between exciting and measuring points of cross FRFs 

% G11IRCSA is the Full FRFs Matrix after removing
G11IRCSA = Experimental_IRCSA(b_190_h11NEW, b_190_h12NEW, L, RB11, RB21, RB12, RB22);
H11IRCSA = reshape(G11IRCSA(1,1,:),1,N); % This is the Displacement/Force FRFs get by Inverse RCSA method 
L11IRCSA = reshape(G11IRCSA(1,2,:),1,N); % This is the Displacement/Torque FRFs
N11IRCSA = reshape(G11IRCSA(2,1,:),1,N); % This is the Angle/Force FRFs
P11IRCSA = reshape(G11IRCSA(2,2,:),1,N); % This is the Angle/Torque FRFs

%% plot
figure(1)

subplot(211)
plot(f,real(c_140_h11NEW),'r',f,real(H11IRCSA),'b','LineWidth',1.5);
xlim([fplot1,fplot2])
ylabel('\it H_{11} \rm m / N');
title('End FRF of 140mm \phi 16mm Beam');
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
legend('Experimental', 'Experimental+IRCSA')
grid minor

subplot(212)
plot(f,imag(c_140_h11NEW),'r',f,imag(H11IRCSA),'b','LineWidth',1.5);
xlim([fplot1,fplot2])
ylabel('\it H_{11} \rm m / N');
xlabel('\it Frequency \rm Hz');
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
legend('Experimental', 'Experimental+IRCSA')
set(gcf,'unit','centimeters','position',[18 5 13.53 9.03],'color','white'); % word 13.5 9
grid minor

%% Save the FRFs at the end of a modified (removed some parts) cantilever beam calculated by Inverse RCSA method
save('IRCSA', 'G11IRCSA', 'H11IRCSA', 'L11IRCSA', 'N11IRCSA', 'P11IRCSA')