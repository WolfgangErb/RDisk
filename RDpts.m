function [rrRC, thRC] = RDpts(m)

% (C) Wolfgang Erb 01.09.2018
% Computes Rhodonea nodes with parameter m1, m2
%-------------------------------------------------------------------------
% INPUT   
% m=[m1,m2]    : parameters of the rhodonea curve
%
% OUTPUT  
% rrRC,thRC    : polar coordinates of the rhodonea nodes

% Determination of points
zrr = linspace(0,1,m(1)+1)*pi/2;
zth = linspace(0,1,4*m(2))*pi/m(2)*(4*m(2)-1)/2;

[RC2,RC1] = meshgrid(zth,zrr);
 
% Selection of used points
[M2,M1] = meshgrid(0:4*m(2)-1,0:m(1));
findM = find(mod(M1+M2+1,2));

rrRC = cos(RC1(findM));
thRC = RC2(findM);

return


