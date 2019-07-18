%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %---------------Copyright------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Jianhui Li                         %
% Time: 07/18/2019                           %
% University of British Columbia, BC, Canada %
% Affiliation:                               %
% Department of Mechanical Engineering       %
% Manufacturing Automation Laboratary        %
% E-mail: jianhui.li@alumni.ubc.ca           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R2b2b] = Experimental_IRCSA(H11, H2a1, s, R11, R21, R12, R22)
%% function of Inverse Receptance Coupling Substructure Analysis for Experimental Data
% Because the H11 or H12 are the only FRFs which can be obtained from experiment
% Here is a method (First Order Difference Method) to get the Full FRFs Matrix 
% from the direct and cross displacement/force FRFs.

% 'H11'(unit,m/N) is the tested FRF, which is a complex value; 's'(unit,m)
% is the distance from the tested point of cross response to the free end

% Parameter 's' is used to find the Full FRFs Matrix at the end of the
% Cantilever Beam based on the experimental direct and cross
% Displacement/Force FRFs

% Notice!!!: 's' can be different from the length of the part to be removed
% Here the s = L is a coincidence
% Refer to:Schmitz, T. L., & Smith, K. S. (2012). Mechanical Vibrations. 
% Springer. Boston, MA: Springer US. https://doi.org/10.1007/978-1-4614-0460-6
% Page 282

N11 = (H11 - H2a1)/s;
L11 = N11;
P11 = N11.*N11./H11;

N = length(H11);

% Assemble the Full FRFs Matrix
G11 = zeros(2,2,N);
G11(1,1,:) = H11;
G11(2,1,:) = N11;
G11(1,2,:) = L11;
G11(2,2,:) = P11;

%% Calculate the Full FRFs at the end of the cantilever beam which length L was removed
R2b2b = IRCSA(G11, R11, R21, R12, R22, N);