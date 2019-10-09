% RDisk: spectral interpolation on rhodonea curves
% Main example for spectral interpolation on rhodonea nodes
% (C) Wolfgang Erb 01.09.2018


clear all
close all

% Path
addpath(genpath('./core/'));

% Set parameters 

m       = [10,9];      % Frequency parameter of Rhodonea curve

nofun   = 9;            % Number of test function [1-10]
parfun  = [1,1];        % Additional parameter of test function

sis     = 'square';     % Form of spectral index set: 'triangle' or 'square'
av      = 0;            % av = 0: no averaging, av = 1 averaging on boundary

N       = 1000;          % Discretization for plot

% Polar coordinates of rhodonea nodes
tic; [rrRD, thRD] = RDpts(m);
        
% Extraction of function values
f = testfundisk(rrRD,thRD,nofun,parfun); 
      
% Computation of realvalued coefficient matrix
G = RDdatM(m,f); 
[~,CR] = RDcfsfft(m,G,sis,av);

% Calculate integral of function
clenshaw = CR(1:4:end,1);
weights = 4*pi./(4-(0:4:2*m(1)).^2);
quad = weights*clenshaw;

% Initialization for plot
x = linspace(-1,1,N);
[X,Y] = meshgrid(x,x);
[theta,r] = cart2pol(X,Y);
idx = r<=1;
Z = nan(size(X));

% Values of the interpolation polynomial for grid and on rhodonea nodes
Sf = RDeval(CR,m,r(idx)',theta(idx)');       
Sfint = RDeval(CR,m,rrRD',thRD'); 
elapsedtime = toc;

% Plot the interpolated function
figure
Z(idx) = Sf;
colormap(jet)
pcolor(x,x,Z), shading interp
set(gca,'XTick',[],'YTick',[])
colorbar
axis square
title('Spectral interpolation on rhodonea nodes','FontName','Avantgarde','FontSize',12)

figure
colormap(copper);
surfl(X, Y, Z);
shading interp
hold on;
plot3(rrRD.*cos(thRD),rrRD.*sin(thRD),f,'o','LineWidth',1,'markersize',5,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k');
title('Spectral interpolation on rhodonea nodes','FontName','Avantgarde','FontSize',12)
hold off
              
% Calculation of the maximal error between function and polynomial interpolation
maxError = max(abs(Sf-testfundisk(r(idx)',theta(idx)',nofun,parfun)));
maxErrorInt = max(abs(f'-Sfint));

fprintf('Main test example for interpolation on rhodonea nodes\n');
fprintf('-----------------------------------------------------------------\n');
fprintf('Elapsed time for interpolation: %10.5f \n',elapsedtime);
fprintf('Maximal error at interpolation points: %16.14f \n',maxErrorInt);
fprintf('Maximal approximation error: %16.14f \n\n',maxError);
fprintf('Integral over the disk: %16.14f \n\n',quad);
        


