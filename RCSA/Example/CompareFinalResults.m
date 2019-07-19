%% Compare the coupled 140mm 16mm + 50mm 10mm Beam with Experimental Results
clc
clear
close all

load IRCSA % The FRFs at the end of the Beam 190mm after removing 50mm
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
G11RCSAFinal = RCSA(G11IRCSA, RD11, RD21, RD12, RD22, N);
H11RCSAFinal = reshape(G11RCSAFinal(1,1,:),1,N); % This is the Displacement/Force FRFs get by RCSA method 
L11RCSAFinal = reshape(G11RCSAFinal(1,2,:),1,N); % This is the Displacement/Torque FRFs
N11RCSAFinal = reshape(G11RCSAFinal(2,1,:),1,N); % This is the Angle/Force FRFs
P11RCSAFinal = reshape(G11RCSAFinal(2,2,:),1,N); % This is the Angle/Torque FRFs

%% Plot
figure(1)
subplot(211)
plot(f,real(d_190_h11NEW),'r',f,real(H11RCSAFinal),'b','LineWidth',1.5);
xlim([fplot1,fplot2])
ylabel('\it H_{11} \rm m / N');
title('End FRF of 140mm \phi 16mm and 50mm \phi 10mm Stepped Beam');
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
legend('Experimental', 'Experimental+IRCSA+RCSA')
grid minor

subplot(212)
plot(f,imag(d_190_h11NEW),'r',f,imag(H11RCSAFinal),'b','LineWidth',1.5);
xlim([fplot1,fplot2])
ylabel('\it H_{11} \rm m / N');
xlabel('\it Frequency \rm Hz');
set(gca,'FontSize', 10 ,'FontName', 'Times New Roman')
legend('Experimental', 'Experimental+IRCSA+RCSA')
set(gcf,'unit','centimeters','position',[18 5 13.53 9.03],'color','white'); % word 13.5 9
grid minor