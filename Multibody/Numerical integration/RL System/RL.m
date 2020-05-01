%% Integration of the differential equations of a RL circuit
% The differential equation is defined in CauchyForm.m

clc; clear;   % the best way to begin a Matlab script
close all;    % close all figures

% This allows to share variables with functions (rled in this case)
global R L U;

% The values of R, L and U are given here and then used in the function
% rled. It is a common way to send parameters to a function
R = 470;
L = 0.001; % the time constant is then 2.13E-6 s
U = 1;
tau = L/R;

% Definition of initial conditions
y0=0; % initial current=0

% Definition of integration time span (integration time interval)
tspan=[0 0.00002];% We give here the initial and final integration times
% It looks short but it is nearly 10 times the time constant, so sufficient
% to reach the steady-state value
% It is also possible to impose the instants where the solution is SAVED by
% tspan=[0:2e-7:2e-5];
% This does not impose the integration time step which is automatically
% adjusted to keep the error below the tolerances

% Launch of integration with ODE45
% rled is the function in which the equations to integrate are specified in
% the form dydt=f(t,y)
[t, y]=ode45(@CauchyForm,tspan,y0);
% After integration, t is a vector array with the time instants for which
% the solution has been computed (the time steps) and y is an array where
% each line contains the value of the state variables at the time indicated
% in the corresponding line of t
% y contains as many columns as the number of state variables (differential
% equations)
% The function could be specified like this as well
% [t y]=ode45('CauchyForm',tspan,y0);
fprintf('Number of time steps used: %3d\n', length(t))

% Theoretical solution (using the time array t issued from integration is
% a pragmatic solution)
t2=t;
% We could have used as well something like
% t2=0:2e-7:2e-5;
y2=U/R*(1-exp(-t2*R/L));

% Plot of the time history of the current and comparison with exact solution
figure,plot(t,y,'-o',t2,y2,'linewidth',2);
xlabel('Time [s]')
ylabel('Current (A)')
title('RL circuit')
legend('Numerical','exact','location','best');
grid on
set(gcf,'unit','centimeters','position',[28 5 13.53 9.03],'color','white');