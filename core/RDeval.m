function Sf = RDeval(CR, m, r, theta)

% (C) Wolfgang Erb 01.09.2018
% Computes the interpolation function at polar coordinates (r,theta) 
%------------------------------------------------------------------
% INPUT    
% m = [m1,m2]  : parameters of the Rhodonea curve
% CR           : coefficient matrix of the interpolation polynomial
% (r,theta)    : evaluation points in polar coordinates
%
% OUTPUT  
% Sf           : evaluated function at point (r,theta)


% Computation of trigonometric polynomials evaluated at (r,theta)
Tr = T(2*m(1),r);
Ttheta = [cos([0:1:2*m(2)]'*theta);sin([2*m(2)-1:-1:1]'*theta)];

% Evaluation of interpolation polynomial via summation
Sf = sum((Tr'*CR).*Ttheta',2)';

return




