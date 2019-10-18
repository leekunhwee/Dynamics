function MultiSemi(Dynamic,Simulation,CuttingCoef,ToolGeo,Cutting,Calc)
%% Parameter 
% Cutting Coefficient
Kt = CuttingCoef.Kt; % N/m^2
Kn = CuttingCoef.Kr;
% Tool geometry
Nt = ToolGeo.Nt;
D  = ToolGeo.D;
% Cutting condition
ae = Cutting.ae;
operation = Cutting.operation; 
% Dynamics data
fx     = Dynamic.f.x*2*pi;
zeta_x = Dynamic.zeta.x;
Px     = Dynamic.shape.x;
fy     = Dynamic.f.y*2*pi;
zeta_y = Dynamic.zeta.y;
Py     = Dynamic.shape.y;
% Build Dynamic Matrix
Cq = diag([zeta_x*2.*fx,zeta_y*2.*fy]);
Kq = diag([fx.*fx      ,      fy.*fy]);
UU=[Px,0,0;
    0,0,Py];
M = length(fx)+length(fy); % Number of modes in both directions

%% Cutting angle
% 1: up-milling, -1: down-milling
if operation == 1  % up-milling
    phist = 0;   % start angle % exit angle
    phiex = acos(1-2*(ae/D));
elseif operation == -1 % down-milling % start angle % exit angle
    phist = acos(2*(ae/D)-1); 
    phiex = pi;
end

%% computational parameters 
k = Calc.k;% number of discretization interval over one period
intk =Calc.intk;% number of numerical integration steps for Equation (37)
m = Calc.m;% since time delay = time period
wa = Calc.wa;% since time delay = time period
wb = Calc.wb;% since time delay = time period

D = zeros(M*m +2*M,M*m +2*M);% matrix D
d = ones(M*m +M, 1); 
d(1 : 2*M) = 0;
D = D + diag(d,-M); 
D((2*M+1):(2*M+M), 1:M)= eye(M); 

%% Simulation Parameters
step_speed = Simulation.step_speed;% steps of spindle speed
step_depth = Simulation.step_depth;% steps of depth of cut
depth_st   = Simulation.depth_st; % starting depth of cut (m) 
depth_fi   = Simulation.depth_fi;  % final depth of cut (m)
speed_st   = Simulation.speed_st;% starting spindle speed (rpm) Int.
speed_fi   = Simulation.speed_fi; % final spindle speed (rpm)
hxx = zeros(1,k); 
hxy = zeros(1,k); 
hyx = zeros(1,k);  
hyy = zeros(1,k); 
fi = zeros(1,intk); 
g = zeros(1,intk);
dtr = 2*pi/Nt/k; % Delta_Phi,if Phi_p = 2 Pai/N

%% numerical integration of specific cutting force coefficient 
for i = 1 : k
    for j = 1 : Nt   % loop for tooth j 
        for h = 1 : intk % loop for numerical integration of hi
            fi(h) = (i-1)*dtr +(j-1)*2*pi/Nt + h*dtr/intk; 
            if (fi(h)>= phist)&&(fi(h)<= phiex) 
                g(h) = 1;% tooth is in the cut
            else
                g(h) = 0;% tooth is out of cut
            end
        end
        hxx(i) = hxx(i)+sum(g.*(Kt.* cos(fi)+Kn.* sin(fi)).* sin(fi))/intk; 
        hxy(i) = hxy(i)+sum(g.*(Kt.* cos(fi)+Kn.* sin(fi)).* cos(fi))/intk; 
        hyx(i) = hyx(i)+sum(g.*(-Kt* sin(fi)+Kn.* cos(fi)).* sin(fi))/intk; 
        hyy(i) = hyy(i)+sum(g.*(-Kt* sin(fi)+Kn.* cos(fi)).* cos(fi))/intk;
    end
end
% start of computation 
for x = 1 : step_speed+1   % loop for spindle speeds
    speed = speed_st +(x-1)*(speed_fi-speed_st)/step_speed;% spindle speed
    tau = 60/speed/Nt;% time delay
    dt = tau/(m);% time step
    for y = 1 : step_depth+1 % loop for depth of cuts
        w = depth_st +(y-1)*(depth_fi-depth_st)/step_depth; % depth of cut 
        % construct transition matrix Fi
        Fi = eye(M*m + 2*M, M*m + 2*M); 
        for i = 1 : m 
            KK = w * UU' * [hxx(i),hxy(i);hyx(i),hyy(i)] * UU;
            L = zeros(2*M, 2*M);      % matrix Ri
            L(1:M, (M+1):2*M) = eye(M); 
            L((M+1):2*M, 1:M) = -Kq-KK;
            L((M+1):2*M, (M+1):2*M)=-Cq;
            R = zeros(2*M, 2*M);      % matrix Li
            R((M+1):2*M, 1:M) = KK;
            P = expm(L*dt);% matrix Pi
            Q = (expm(L*dt)-eye(2*M))/L * R; % matrix Ri 
            D(1 : 2*M, 1 : 2*M) = P; 
            D(1 : 2*M,(M*m + 1) : (M*m + M)) = wa*Q(1 : 2*M, 1 : M); 
            D(1 : 2*M,(M*m + M+1) : (M*m + 2*M)) = wb*Q(1 : 2*M, 1 : M); 
            Fi = D*Fi;  % transition matrix Phi
        end
        ss(x, y) = speed;   % matrix of spindle speeds
        dc(x, y) = w * 1000;   % matrix of depth of cuts
        ei(x, y) = max(abs(eig(Fi)));   % matrix of eigenvalues
    end
    step_speed+1-x
end
figure
contour(ss,dc,ei,[1, 1],'k')