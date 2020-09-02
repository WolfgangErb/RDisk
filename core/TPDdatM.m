function G = TPDdatM(m,f)

% (C) Wolfgang Erb 01.09.2018
%     Version 0.4: 31.08.2020

% Generates the data matrix from the function evaluation at tensor product
% polar grid
%--------------------------------------------------------
% INPUT    
% m = [m1,m2]  : parameters of rhodonea curve
% f            : function values at TP-grid

% Output
% G            : (4 m1) x ( 4 m2) data matrix

% Generate extended grid G via symmetrization
% Go from (m1+1) x (4m2) to (4m1) x (4m2)

G_TP = reshape(f,m(1)+1,4*m(2))/m(1)/m(2)/8;

% Extend data Matrix to the index set K^{(m)} by using the flip operation
G_TPstar = G_TP(end:-1:1,:);
G_TPstar = circshift(G_TPstar',2*m(2))';
G_TPstar = G_TPstar(2:end,:);
G_Km = [G_TP;G_TPstar];

% Extend data Matrix symmetrically to the whole index set J^{(m)} 
G_Kmcross = G_Km(end-1:-1:2,:);
G = [G_Km;G_Kmcross];

return