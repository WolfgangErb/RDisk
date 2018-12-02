% Plot Lagrange basis functions for rhodonea nodes 
% (C) Wolfgang Erb 01.09.2018

% Parameters 

clear all
close all

m = [30,31];              % Parameters of spherical Lissajous curve
N = 200;                  % Discretization for plot
K = 12;                   % Index for Lagrange basis function

% Parameters for spectral index set
sis = 'square';           % 'square' or 'triangle'
av = 0;                   % average on boundary or not

% Initialization for plot
x = linspace(-1,1,N);
[X,Y] = meshgrid(x,x);
[theta,r] = cart2pol(X,Y);
idx = r<=1;
Z = nan(size(X));

% Coordinates and weights of rhodonea nodes
tic; [rrRD, thRD] = RDpts(m);
        
% Extraction of Lagrange function values
f = zeros(length(thRD),1);
f(K) = 1;

% Computation of Coefficient Matrix
G = RDdatM(m,f); 
[~,CR] = RDcfsfft(m,G,sis,av);

% Values of the interpolation polynomial
Sf = RDeval(CR,m,r(idx)',theta(idx)');       
elapsedtime = toc;

% Plot the interpolated function
figure(1)
Z(idx) = Sf;
colormap(jet)
pcolor(x,x,Z), shading interp
set(gca,'XTick',[],'YTick',[])
axis square
title(['Lagrange basis function for rhodonea nodes'])