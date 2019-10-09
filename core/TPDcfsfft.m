function [CC,CR] = TPDcfsfft(m,G)

% (C) Wolfgang Erb 01.09.2018
% Computes the coefficient matrix of the interpolating polynomial
% For data given on a tensor-product grid
% ----------------------------------------------------------------------
% INPUT    
% m = [m1,m2]  : parameters of tensor-product grid
% G            : (4 m1) x (4 m2) data matrix

% Output 
% CC           : (2 m1+1) x (4 m2) coefficient matrix (complex)
% CR           : (2 m1+1) x (4 m2) coefficient matrix (real)


% Perform 2D Fourier transform of G
Gh = fft2(G);

% Generate mask for coefficients. 
Mask = 2*ones(2*m(1)+1,4*m(2));
Mask(1,1:end) = Mask(1,1:end)/2;
Mask(2*m(1)+1,1:end) = Mask(2*m(1)+1,1:end)/2;
Mask(1:end,1) = Mask(1:end,1)/2;
Mask(1:end,2*m(2)+1) = Mask(1:end,2*m(2)+1)/2;

% Coefficients for complex valued basis
CC = Gh(1:2*m(1)+1,:).*Mask;

% Coefficients for real valued basis
CR = [real(CC(:,1:2*m(2)+1)),imag(CC(:,2*m(2)+2:4*m(2)))];
    
return