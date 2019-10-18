function [f,w] = FrequencySequence(Simulation)
f_start = Simulation.f_start;
f_end   = Simulation.f_end;
df      = Simulation.df;
f = (f_start:df:f_end)';   % frequency, Hz
w = f*2*pi;                % frequency, rad