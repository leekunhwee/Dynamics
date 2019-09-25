%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------Copyright--------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Jianhui Li                         %
% Time: 09/22/2019                           %
% University of British Columbia, BC, Canada %
% Affiliation:                               %
% Department of Mechanical Engineering       %
% Manufacturing Automation Laboratary        %
% E-mail: jianhui.li@alumni.ubc.ca           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Zero-Order Solution            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all

load FRF

%% Cutting force coefficients
Kt = 1319.4e6; % N/m^2
Kn = 788.8/1319.4;
Nt = 2;

%% Cutting condition
D=31.75;       % mm
ae=D/2;
operation = 1; % down-milling

%% Cutting angle
phis = acos((D/2 - ae)/(D/2));  % (rad) the immersion angle

if operation == 0             % Up milling 
    phist = 0;                % Start angle rad
    phiex = phist + phis;     % Exit angle rad
else                          % Down milling
    phiex = pi;               % Start angle rad
    phist = phiex - phis;     % Exit angle rad
end

%% Average angle parameters
alphaxx = 0.5*(( cos(2*phiex)-2*Kn*phiex+Kn*sin(2*phiex))-( cos(2*phist)-2*Kn*phist+Kn*sin(2*phist)));
alphaxy = 0.5*((-sin(2*phiex)-2*phiex+Kn*cos(2*phiex))   -(-sin(2*phist)-2*phist+Kn*cos(2*phist)));
alphayx = 0.5*((-sin(2*phiex)+2*phiex+Kn*cos(2*phiex))   -(-sin(2*phist)+2*phist+Kn*cos(2*phist)));
alphayy = 0.5*((-cos(2*phiex)-2*Kn*phiex-Kn*sin(2*phiex))-(-cos(2*phist)-2*Kn*phist-Kn*sin(2*phist)));

%% Initialization
w = (f_start:df:f_end)'*2*pi;   % frequency, rad/s
FRFxx = FRFX;
FRFyy = FRFY;
lambda1=zeros(1,length(w));
lambda2=zeros(1,length(w));

%% Calculate Lambda
for cnt = 1:length(w)
    % Oriented FRF (Ignore the cross FRFs)
    FRF_or = [alphaxx*FRFxx(cnt) alphaxy*FRFyy(cnt); alphayx*FRFxx(cnt) alphayy*FRFyy(cnt)];    % m/N
    
    % Calculate two eigenvalues -- lambda
    % Note: the definations of lambda and Lambda are slightly different
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reference£º
    % T.L. Schmitz, K.S. Smith, Machining Dynamics, springer, Springer US, 
    % Boston, MA, 2009. https://doi.org/10.1007/978-0-387-09645-2.
    % P133
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Calculate eigenvalues value for the standard form by eig() function
    E = eig(FRF_or);
    temp = E(1);
    lambda1(cnt) = temp;
    temp = E(2);
    lambda2(cnt) = temp;
    % Because the eig() sorts the eigenvalue by absolute value 
    if (cnt > 1)
        dot_prod1 = real(lambda2(cnt))*real(lambda2(cnt-1)) + imag(lambda2(cnt))*imag(lambda2(cnt-1));
        dot_prod2 = real(lambda2(cnt))*real(lambda1(cnt-1)) + imag(lambda2(cnt))*imag(lambda1(cnt-1));
        % This part is really important. For conjugate pair, the value of dot product means the projection 
        % which can be used to find the nearest next point. See the picture in README.md. 
        if (dot_prod2 > dot_prod1)
            temp = lambda2(cnt);
            lambda2(cnt) = lambda1(cnt);
            lambda1(cnt) = temp;
        end
    end
end

lambda1 = lambda1';
lambda2 = lambda2';

% Refer to Machining Dynamics P133
blim1 = (2*pi/Nt/Kt)./((real(lambda1)).^2 + (imag(lambda1)).^2) .* (real(lambda1) .* (1 + (imag(lambda1)./real(lambda1)).^2));  % m
blim2 = (2*pi/Nt/Kt)./((real(lambda2)).^2 + (imag(lambda2)).^2) .* (real(lambda2) .* (1 + (imag(lambda2)./real(lambda2)).^2));

[index1] = find(blim1 > 0);
blim1 = blim1(index1);
blim1 = blim1 * 1e3;      % mm
w1 = w(index1);
psi1 = atan2(imag(lambda1), real(lambda1));
psi1 = psi1(index1);
epsilon1 = pi - 2 * psi1;

[index2] = find(blim2 > 0);
blim2 = blim2(index2);
blim2 = blim2 * 1e3;
w2 = w(index2);
psi2 = atan2(imag(lambda2), real(lambda2));
psi2 = psi2(index2);
epsilon2 = pi - 2 * psi2;

lobeNumber = 15;
omega1 = {lobeNumber,1};
for k = 1:lobeNumber
    omega1{k} = (60/Nt)*w1./(epsilon1 + 2*(k-1)*pi);
end

omega2 = {lobeNumber,1};
for k = 1:lobeNumber
    omega2{k} = (60/Nt)*w2./(epsilon2 + 2*(k-1)*pi);
end

figure(1)
hold on
for k = 1:lobeNumber
    plot(omega1{k},blim1, 'r','linewidth',2);
    plot(omega2{k},blim2, 'b','linewidth',2);
end

axis([0 16000 0 6])
xlabel('\it\Omega \rm/ rpm')
ylabel('\itb_{lim} \rm/ mm')
set(gca,'FontSize', 11 ,'FontName', 'Times New Roman')
set(gcf,'unit','centimeters','position',[28 5 13.53 9.03],'color','white');% word£¨13.5,9£©
set(gca,'xtick',[2000 4000 6000 8000 10000 12000 14000 16000])
grid on