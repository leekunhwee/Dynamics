function [R11, R21, R12, R22] = Beam_FRF(length, d_out, d_in, f1, f2, r, density_out, density_in, E_out, E_in, nu_out, nu_in, el)
%% function to get FRFs of a beam in Free-Free conditon with timoshenko beam element (Numerical method)
% This function can be used to calculate the Hollow Beam by considering the Externer Diameter and Inner Diameter
% While timo_free_free.m only can be used to calculate the Solid Cylindrical Beam
% 'l'(unit,mm) is the length of the beam.
% 'd'(unit,mm) is the diameter of the beam.
% [f1,f2](unit,HZ) is the range of calculated frequency, and 'r'(unit,HZ) is the resolution of the range, and 'num'(unit,1) is the number of natural frequencies you want to obtain.
% 'rho'(unit,kg/m^3) is the density of the tool, and 'E'(unit,pa) is the modulus of the tool, and 'nu_v'(unit,1) is Poissons ratio.
% [nf](unit,1) is the natural frequencies of the model.
% 'el' is the maxmun length of the Timoshenko element 

f = f1:r:f2;
%  Define beam parameters
eta = 0.001;                                          % structural damping
I_out = pi/64*(d_out^4-d_in^4);                       % second area moment of inertia for outer diameter, m^4
I_in = pi/64*(d_in^4);                                % second area moment of inertia for inner diameter, m^4
EIeq  = E_out*I_out + E_in*I_in;                      % composite structural stiffness, N-m^2
EI = EIeq*(1+1i*eta);                                 % damped structural stiffness, N-m^2
Aout = pi/4*(d_out^2 - d_in^2);                       % cross-sectional area of outer diameter, m^2
Ain = pi/4*(d_in^2);                                  % cross-sectional area of inner diameter, m^2

% Timoshenko terms
n = ceil(length/el);                                  % number of finite elements
Gout = E_out/(2*(1+nu_out));                          % shear modulus of outer diameter, N/m^2
Gin = E_in/(2*(1+nu_in));                             % shear modulus of inner diameter, N/m^2
A = Aout + Ain;                                       % composite cross-sectional area, m^2
AG = Gout*Aout + Gin*Ain;                             % produce of area and shear modulus for composite beam, N
nu = (nu_out*Aout + nu_in*Ain)/A;                     % composite Poissons ratio
I = I_out + I_in;                                     % composite second area moment of inertia, m^4
rg = (I/A)^0.5;                                       % radius of gyration, m

mpl = density_out*Aout + density_in*Ain;              % mass per unit length, kg/m
    if (d_in == 0) || ((d_in ~= 0) && (E_in ~= 0))    % solid cross-section
        kp = 6*(1+nu)^2/(7+12*nu+4*nu^2);             % shape factor
    else                                              % hollow cross-section
        num = 6*(1 + nu_out)^2*(d_in^2 + d_out^2)^2;
        den = 7*d_in^4 + 34*d_in^2*d_out^2 + 7*d_out^4 + nu_out*(12*d_in^4 + 48*d_in^2*d_out^2 + 12*d_out^4) + nu_out^2*(4*d_in^4 + 16*d_in^2*d_out^2 + 4*d_out^4);
        kp = num/den;
    end

% Determine free-free receptance of beam
[R11, R21, R12, R22] = timo_free_free(f, EI, length, AG, kp, rg, mpl, n);