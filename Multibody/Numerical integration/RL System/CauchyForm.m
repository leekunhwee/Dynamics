% File CauchyForm.m
% This file implements the differential equation governing the behaviour of a RL circuit
% L di/dt + R i =U
% with i the current, L the inductance, R the resistance and U the supply voltage 
% This function is called by program RL.m

% Definition of function
% Be careful the name of the function should be the one of the file
function dydt = CauchyForm(t,y)
global R L U; 
dydt = (U-R*y)/L; 
end
