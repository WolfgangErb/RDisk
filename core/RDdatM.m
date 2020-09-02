function G = RDdatM(m,f)

% (C) Wolfgang Erb 01.09.2018
%     Version 0.4: 31.08.2020

% Generates the data matrix from the function evaluations
% at the rhodonea nodes
%--------------------------------------------------------
% INPUT    
% m = [m1,m2]  : parameters of rhodonea curve
% f            : function values at RC points

% Output
% G            : (4 m1) x ( 4 m2) data matrix

% Generation of Data Matrix
[M2,M1] = meshgrid(0:4*m(2)-1,0:m(1));
findM = find(mod(M1+M2+1,2));

% Generate part of data Matrix related to the index set I^{(m)}
% The function values for the center indices i_1 = m_1 are added already here
G_Im = zeros((m(1)+1)*4*m(2),1);
G_Im(findM) = f/m(1)/m(2)/8;
G_Im = reshape(G_Im,size(M1));

% Extend data Matrix to the index set K^{(m)} by using the flip operation
G_Imstar = G_Im(end:-1:1,:);
G_Imstar = circshift(G_Imstar',2*m(2))';
G_Imstar = G_Imstar(2:end,:);
G_Km = [G_Im;G_Imstar];

% Extend data Matrix symmetrically to the whole index set J^{(m)}  
G_Kmcross = G_Km(end-1:-1:2,:);
G = [G_Km;G_Kmcross];

%Note:
%G is the data matrix of function evalutions of g on the index set J^{(m)}
%the rows 1,...,2m(1)+2 of G correspond to the indices j_1 = 0,...,2m(1)+1,
%the rows 2*m(1)+2,...,4m(1) to the negative indices j_1 = -2m(1)+1,...-1.

return