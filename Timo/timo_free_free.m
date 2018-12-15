% timo_free_free.m
% This program uses n Timoshenko beam elements to determine end receptances
% for free-free beam.
% Input variables are: f, frequency, Hz; EI,
% structural rigidity including hysteretic damping, N-m^2;
% L, beam length, m; density, kg/m^3; A, cross sectional area, 
% m^2; kp, shear factor; G, shear modulus, N/m^2;
% rg, radius of gyration; n, number of elements

function [R11, R21, R12, R22] = timo_free_free(f, EI, L, AG, kp, rg, mpl, n)
l = L/n;                    % length of each finite element
phi = 12*EI/(kp*AG*l^2);    % shear deformation parameter

% Single element matrices for Timoshenko beam 
% Mass matrix
Mt = mpl*l/(1+phi)^2*[(13/35+7*phi/10+phi^2/3) (11/210+11*phi/120+phi^2/24)*l (9/70+3*phi/10+phi^2/6) -(13/420+3*phi/40+phi^2/24)*l;
    (11/210+11*phi/120+phi^2/24)*l (1/105+phi/60+phi^2/120)*l^2 (13/420+3*phi/40+phi^2/24)*l -(1/140+phi/60+phi^2/120)*l^2;
    (9/70+3*phi/10+phi^2/6) (13/420+3*phi/40+phi^2/24)*l (13/35+7*phi/10+phi^2/3) -(11/210+11*phi/120+phi^2/24)*l;
    -(13/420+3*phi/40+phi^2/24)*l -(1/140+phi/60+phi^2/120)*l^2 -(11/210+11*phi/120+phi^2/24)*l (1/105+phi/60+phi^2/120)*l^2];

Mr = mpl*l/(1+phi)^2*(rg/l)^2*[6/5 (1/10-phi/2)*l -6/5 (1/10-phi/2)*l;
    (1/10-phi/2)*l (2/15+phi/6+phi^2/3)*l^2 -(1/10-phi/2)*l -(1/30+phi/6-phi^2/6)*l^2;
    -6/5 -(1/10-phi/2)*l 6/5 -(1/10-phi/2)*l;
    (1/10-phi/2)*l -(1/30+phi/6-phi^2/6)*l^2 -(1/10-phi/2)*l (2/15+phi/6+phi^2/3)*l^2];

M = Mt + Mr;

% Stiffness matrix，刚度阵
Kb = EI/(l^3*(1+phi)^2)*[12 6*l -12 6*l;
    6*l (4+2*phi+phi^2)*l^2 -6*l (2-2*phi-phi^2)*l^2;
    -12 -6*l 12 -6*l;
    6*l (2-2*phi-phi^2)*l^2 -6*l (4+2*phi+phi^2)*l^2];

Ks = kp*AG*phi^2/(4*l*(1+phi)^2)*[4 2*l -4 2*l;
    2*l l^2 -2*l l^2;
    -4 -2*l 4 -2*l;
    2*l l^2 -2*l l^2];

K = Kb + Ks;

Mtemp2 = M;
Ktemp2 = K;

% Build full mass and stiffness matrices
for cnt = 2:n
    % Concatenate left element matrices with required zeros
    right = zeros(cnt*2, 2);
    bottom = zeros(2, (cnt+1)*2);
    % Mass matrix
    temp = cat(2, M, right);
    Mtemp1 = cat(1, temp, bottom);
    % Stiffness matrix
    temp = cat(2, K, right);
    Ktemp1 = cat(1, temp, bottom);
    
    % Concatenate right element matrices with required zeros
    left = zeros(cnt*2, 2);
    top = zeros(2, (cnt+1)*2);
    % Mass matrix
    temp = cat(2, left, Mtemp2);
    Mtemp2 = cat(1, top, temp);
    % Stiffness matrix
    temp = cat(2, left, Ktemp2);
    Ktemp2 = cat(1, top, temp);
    
    % Add two matrices
    M = Mtemp1 + Mtemp2;
    K = Ktemp1 + Ktemp2;
end

% Calculate required direct and cross receptances for ends of beam 梁的原点和跨点频响
t1 = length(f);
R11 = zeros(2,2,t1);
R21 = zeros(2,2,t1);
R12 = zeros(2,2,t1);
R22 = zeros(2,2,t1);

for cnt = 1:t1
    w = f(cnt)*2*pi;        % frequency, rad/s
    D = inv(-M*w^2 + K);    % dynamic matrix
    R11(1,1,cnt) = D(1,1);
    R11(1,2,cnt) = -D(1,2);
    R11(2,1,cnt) = -D(2,1);
    R11(2,2,cnt) = D(2,2);
    
    R12(1,1,cnt) = D(1,2*(n+1)-1);
    R12(1,2,cnt) = -D(1,2*(n+1));
    R12(2,1,cnt) = -D(2,2*(n+1)-1);
    R12(2,2,cnt) = D(2,2*(n+1));
    
    R21(1,1,cnt) = D(2*(n+1)-1,1);
    R21(1,2,cnt) = -D(2*(n+1)-1,2);
    R21(2,1,cnt) = -D(2*(n+1),1);
    R21(2,2,cnt) = D(2*(n+1),2);
    
    R22(1,1,cnt) = D(2*(n+1)-1,2*(n+1)-1);
    R22(1,2,cnt) = -D(2*(n+1)-1,2*(n+1));
    R22(2,1,cnt) = -D(2*(n+1),2*(n+1)-1);
    R22(2,2,cnt) = D(2*(n+1),2*(n+1));
    
    clear D;
end

