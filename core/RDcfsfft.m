function [CC,CR] = RDcfsfft(m,G,s,av)

% (C) Wolfgang Erb 01.09.2018
% Computes the coefficient matrix of the interpolating polynomial
% ----------------------------------------------------------------------
% INPUT    
% m = [m1,m2]  : parameters of rhodonea curve
% G            : (4 m1) x (4 m2) data matrix
% s            : Form of spectral index set: 'square' or 'triangle'
% av           : av = 0: no averaging, av = 1: averaging on boundary

% Output 
% CC           : (2 m1+1) x (4 m2) coefficient matrix (complex)
% CR           : (2 m1+1) x (4 m2) coefficient matrix (real)

% 2D Fourier transform of G
Gh = fft2(G);

% Generate mask for coefficients. This depends on the form of the spectral
% index set (rectangular or triangular) and whether we average the
% coefficients on the boundary or not. 

[M2,M1] = meshgrid(0:4*m(2)-1,0:4*m(1)-1);

if strcmp(s,'triangle') && (av==0)

Mask = 2*(double(M1*m(2)+M2*m(1)<2*m(1)*m(2)) + double(M1*m(2)+(4*m(2)-M2)*m(1)<2*m(1)*m(2)));
Meq  = 2*double(M1*m(2)+M2*m(1)==2*m(1)*m(2)).*double(M1 >= m(1))+ 2*double(M1*m(2)+(4*m(2)-M2)*m(1)==2*m(1)*m(2)).*double(M1 > m(1));
Mask = (Mask + Meq).*(1-mod(M1+M2,2));
Mask(1,1:end) = Mask(1,1:end)/2;
Mask(1:end,1) = Mask(1:end,1)/2;
Mask(2*m(1)+1,1) = Mask(2*m(1)+1,1)/2;
Mask(m(1)+1,m(2)+1) = Mask(m(1)+1,m(2)+1)/2;

elseif strcmp(s,'square') && (av==0)
    
Mask = 2*(double(M2 < m(2)) + double((4*m(2)-M2)<m(2)));
Meq  = 2*double(M2==m(2));
Mask = (Mask + Meq).*(1-mod(M1+M2,2));
Mask(1,1:end) = Mask(1,1:end)/2;
Mask(2*m(1)+1,1:end) = Mask(2*m(1)+1,1:end)/2;
Mask(1:end,1) = Mask(1:end,1)/2;
Mask(m(1)+1,m(2)+1) = Mask(m(1)+1,m(2)+1)/2;

elseif strcmp(s,'triangle') && (av==1)

Mask = 2*(double(M1*m(2)+M2*m(1)<2*m(1)*m(2)) + double(M1*m(2)+(4*m(2)-M2)*m(1)<2*m(1)*m(2)));
Meq  = double(M1*m(2)+M2*m(1)==2*m(1)*m(2))+ double(M1*m(2)+(4*m(2)-M2)*m(1)==2*m(1)*m(2));
Mask = (Mask + Meq).*(1-mod(M1+M2,2));
Mask(1,1:end) = Mask(1,1:end)/2;
Mask(1:end,1) = Mask(1:end,1)/2;
Mask(1,2*m(2)+1) = 0;

elseif strcmp(s,'square') && (av==1)
    
Mask = 2*(double(M2 < m(2)) + double((4*m(2)-M2)<m(2)));
Meq  = double(M2==m(2))+ double((4*m(2)-M2)==m(2));
Mask = (Mask + Meq).*(1-mod(M1+M2,2));
Mask(1,1:end) = Mask(1,1:end)/2;
Mask(2*m(1)+1,1:end) = Mask(2*m(1)+1,1:end)/2;
Mask(1:end,1) = Mask(1:end,1)/2;

end

% Coefficients for complex valued basis
CC = Gh.*Mask;
CC = CC(1:2*m(1)+1,:);

% Coefficients for real valued basis
CR = [real(CC(:,1:2*m(2)+1)),imag(CC(:,2*m(2)+2:4*m(2)))];

if strcmp(s,'square') && (av==0)
  CR(m(1)+2:2*m(1)+1,m(2)+1) = zeros(m(1),1);
  CR(m(1)+2:2*m(1)+1,3*m(2)+1) = -imag(CC(m(1)+2:2*m(1)+1,m(2)+1)); 
end
       
return