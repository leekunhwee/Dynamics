function ZOA (FRF,CuttingCoef,ToolGeo,Cutting)
%% Parameter 
% Simulation parameter
w  = FRF.w;% frequency, rad/s
% f  = FRF.f;% frequency, Hz
% Cutting Coefficient
Kt = CuttingCoef.Kt; % N/m^2
Kn = CuttingCoef.Kr/CuttingCoef.Kt;
% Tool geometry
Nt = ToolGeo.Nt;
D  = ToolGeo.D;
% Cutting condition
ae = Cutting.ae;
operation = Cutting.operation;

%% Cutting angle
% 1: up-milling, -1: down-milling
if operation == 1  % up-milling
    phist = 0;   % start angle % exit angle
    phiex = acos(1-2*(ae/D));
elseif operation == -1 % down-milling % start angle % exit angle
    phist = acos(2*(ae/D)-1); 
    phiex = pi;
end

%% Average angle parameters
alphaxx = 0.5*(( cos(2*phiex)-2*Kn*phiex+Kn*sin(2*phiex))-( cos(2*phist)-2*Kn*phist+Kn*sin(2*phist)));
alphaxy = 0.5*((-sin(2*phiex)-2*phiex+Kn*cos(2*phiex))   -(-sin(2*phist)-2*phist+Kn*cos(2*phist)));
alphayx = 0.5*((-sin(2*phiex)+2*phiex+Kn*cos(2*phiex))   -(-sin(2*phist)+2*phist+Kn*cos(2*phist)));
alphayy = 0.5*((-cos(2*phiex)-2*Kn*phiex-Kn*sin(2*phiex))-(-cos(2*phist)-2*Kn*phist-Kn*sin(2*phist)));

%% Initialization
FRFxx = FRF.X;
FRFyy = FRF.Y;
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

figure
hold on
for k = 1:lobeNumber
    plot(omega1{k},blim1, 'r','linewidth',2);
    plot(omega2{k},blim2, 'b','linewidth',2);
end
