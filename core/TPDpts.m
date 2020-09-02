function [rrRC, thRC] = TPDpts(m)

% (C) Wolfgang Erb 01.09.2018
%     Version 0.4: 31.08.2020

% Computes tensor product grid of dimension (m1+1) x 4m2
%-------------------------------------------------------------------------
% INPUT   
% m=[m1,m2]    : parameters of the polar grid
%
% OUTPUT  
% rrRC,thRC    : polar coordinates of the grid

% Determination of points
zrr = linspace(0,1,m(1)+1)*pi/2;
zth = linspace(0,1,4*m(2))*pi/m(2)*(4*m(2)-1)/2;

[RC2,RC1] = meshgrid(zth,zrr);
 
rrRC = cos(RC1(:));
thRC = RC2(:);

return


