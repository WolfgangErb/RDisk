function G = TPDdatM(m,f)

% (C) Wolfgang Erb 01.09.2018
% Generates the data matrix from the function evaluations
%--------------------------------------------------------
% INPUT    
% m = [m1,m2]  : parameters of rhodonea curve
% f            : function values at TP-grid

% Output
% G            : (4 m1) x ( 4 m2) data matrix

% Generate extended grid G via symmetrization
% Go from (m1+1) x (4m2) to (4m1) x (4m2)

Gh = reshape(f,m(1)+1,4*m(2))/m(1)/m(2)/8;

Gh2 = Gh(end:-1:1,:);
Gh2 = circshift(Gh2',2*m(2))';
Gh2 = Gh2(2:end,:);

Gh3 = [Gh;Gh2];
Gh4 = Gh3(end-1:-1:2,:);

G = [Gh3;Gh4];

return