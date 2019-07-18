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

function R0 = IRCSA(G11, R11, R21, R12, R22, N)
%% Inverse Receptance Coupling Substructure Analysis
% Original end point FRFs G11 Removing a free-free beam with Direct and
% Cross FRFs including R11, R21, R12 and R22.
% Return the FRFs of the Removed end point R0
% N is the number of frequency data points

R0 = zeros(2,2,N);
for cnt = 1:N
    R0(:,:,cnt) = R21(:,:,cnt)/(R11(:,:,cnt) - G11(:,:,cnt))*R12(:,:,cnt) - R22(:,:,cnt);
end