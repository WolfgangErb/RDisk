function G = RDdatM(m,f)

% (C) Wolfgang Erb 01.09.2018
% Generates the data matrix from the function evaluations
%--------------------------------------------------------
% INPUT    
% m = [m1,m2]  : parameters of rhodonea curve
% f            : function values at RC points

% Output
% G            : (4 m1) x ( 4 m2) data matrix

% Generation of Data Matrix
[M2,M1] = meshgrid(0:4*m(2)-1,0:m(1));
findM = find(mod(M1+M2+1,2));

Gh = zeros((m(1)+1)*4*m(2),1);
Gh(findM) = f/m(1)/m(2)/4;
Gh = reshape(Gh,size(M1));

Gh2 = Gh(end:-1:1,:);
Gh2 = circshift(Gh2',2*m(2))';
Gh2 = Gh2(2:end,:);

Gh3 = [Gh;Gh2];
Gh4 = Gh3(end-1:-1:2,:);

G = [Gh3;Gh4];

return