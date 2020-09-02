function [CC,CR] = TPDcfsfft(m,G)

% (C) Wolfgang Erb 01.09.2018
%     Version 0.4: 31.08.2020

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
Gh = Gh(1:2*m(1)+1,:);

% Generate mask for complex valued coefficients. 
MaskC = ones(2*m(1)+1,4*m(2));
MaskC(1,1:end) = MaskC(1,1:end)/2;
MaskC(2*m(1)+1,1:end) = MaskC(2*m(1)+1,1:end)/2;

% Generate mask for real valued coefficients.
MaskR = 2*MaskC;
MaskR(1:end,1) = MaskR(1:end,1)/2;
MaskR(1:end,2*m(2)+1) = MaskR(1:end,2*m(2)+1)/2;

% Coefficients for complex valued basis
CC = Gh.*MaskC;

% Coefficients for real valued basis
CR = [real(Gh(:,1:2*m(2)+1)),-imag(Gh(:,2*m(2)+2:4*m(2)))];
CR = CR.*MaskR;
    
return