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

function G11 = RCSA(R0, R11, R21, R12, R22, N)
%% Receptance Coupling Substructure Analysis
% Original end point FRFs R0 Coupling a free-free beam with Direct and
% Cross FRFs including R11, R21, R12 and R22.
% Return the FRFs of the Coupled end point G11
% N is the number of frequency data points

G11 = zeros(2,2,N);
for cnt = 1:N
	G11(:,:,cnt) = R11(:,:,cnt) - R12(:,:,cnt) * ((R0(:,:,cnt) + R22(:,:,cnt)) \ R21(:,:,cnt));
end