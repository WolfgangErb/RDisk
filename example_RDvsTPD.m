% RDisk: spectral interpolation on rhodonea curves
% Comparison of rhodonea interpolation (RD scheme with a triangular spectral index set)
% on the rhodonea nodes (RDpts) 
% and a tensor product interpolation (TPD scheme with a rectangular spectral index set) 
% on a rectangular polar grid (TPDpts)
% (C) Wolfgang Erb 01.09.2018


clear all
close all

% Path
addpath(genpath('./core/'));

% Set parameters 

nofun   = 10;            % Number of test function [1-9]
parfun  = [1,1];        % Additional parameter of test function

av      = 0;            % av = 0: no averaging, av = 1 averaging on boundary
N       = 200;          % Discretization for plot

% Initialization for plot
x = linspace(-1,1,N);
[X,Y] = meshgrid(x,x);
[theta,r] = cart2pol(X,Y);
idx = r<=1;
Z = nan(size(X));


m = 2:2:40;               % Range of frequency parameters

maxError1 = zeros(length(m),1);
maxError2 = zeros(length(m),1);
np1 = zeros(length(m),1);
np2 = zeros(length(m),1);

for i = 1:length(m)

mm1 = [m(i),m(i)];
mm2 = [m(i),m(i)];
% Polar coordinates of rhodonea nodes
[rrRD1, thRD1] = RDpts(mm1);
[rrRD2, thRD2] = TPDpts(mm2);

np1(i) = length(rrRD1);
np2(i) = length(rrRD2);
        
% Extraction of function values
f1 = testfundisk(rrRD1,thRD1,nofun,parfun); 
f2 = testfundisk(rrRD2,thRD2,nofun,parfun); 
      
% Computation of realvalued coefficient matrix
G1 = RDdatM(mm1,f1);
G2 = TPDdatM(mm2,f2); 
[~,CR1] = RDcfsfft(mm1,G1,'triangle',av);
[~,CR2] = TPDcfsfft(mm2,G2);

% Values of the interpolation polynomial for grid and on rhodonea nodes
Sf1 = RDeval(CR1,mm1,r(idx)',theta(idx)'); 
Sf2 = RDeval(CR2,mm2,r(idx)',theta(idx)'); 
            
% Calculation of the maximal error between function and polynomial interpolation
maxError1(i) = max(abs(Sf1-testfundisk(r(idx)',theta(idx)',nofun,parfun)));
maxError2(i) = max(abs(Sf2-testfundisk(r(idx)',theta(idx)',nofun,parfun)));

end

% Plot the last interpolated function
figure
Z(idx) = Sf1;
colormap(jet)
pcolor(x,x,Z), shading interp
set(gca,'XTick',[],'YTick',[])
colorbar
axis square
title('Spectral interpolation (RD scheme)','FontName','AvantGarde','FontSize',12)

figure

h1 = plot(m,log10(maxError1),'o','LineWidth',1,'markersize',5,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k');
hold on 
h2 = plot(m,log10(maxError2),'o','LineWidth',1,'markersize',5,'MarkerFaceColor','r','MarkerEdgeColor','k');

hTitle  = title ('Comparison of RD scheme with tensor-product scheme');
hXLabel = xlabel('Degree $$m$$','Interpreter','latex' );
hYLabel = ylabel('$$\log_{\;\,10}\, (\mathrm{max error}\;\, ) $$','Interpreter','latex' );

hLegend = legend( ...
  [h1, h2], ...
  'Total degree (RD scheme)' ,  ...
  'Maximal degree (TPD scheme)'      , ...
  'location', 'NorthEast');

set( gca                    , ...
  'FontName'   , 'Helvetica', ...
  'FontSize'   ,   12       , ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.2 .2 .2], ...
  'YColor'      , [.2 .2 .2], ...
  'LineWidth'   , 1         );

set([hTitle, hXLabel, hYLabel], ...
    'FontName'   , 'AvantGarde' ,...
    'FontSize'   ,   12  );
set( hLegend             , ...
    'FontSize'   , 12    );
hold off  

%print('ComparisonDegreeFunction1','-dpng','-r600');

figure

h1 = plot(np1,log10(maxError1),'o','LineWidth',1,'markersize',5,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k');
hold on 
h2 = plot(np2,log10(maxError2),'o','LineWidth',1,'markersize',5,'MarkerFaceColor','r','MarkerEdgeColor','k');

hTitle  = title ('Comparison of RD scheme with tensor-product scheme');
hXLabel = xlabel('Number of nodes','Interpreter','latex' );
hYLabel = ylabel('$$\log_{\;\,10}\, (\mathrm{max error}\;\, ) $$','Interpreter','latex' );

hLegend = legend( ...
  [h1, h2], ...
  'Total degree (RD scheme)' ,  ...
  'Maximal degree (TPD scheme)' , ...
  'location', 'NorthEast');

set( gca                    , ...
  'FontName'   , 'Helvetica', ...
  'FontSize'   ,   12       , ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.2 .2 .2], ...
  'YColor'      , [.2 .2 .2], ...
  'LineWidth'   , 1         );

set([hTitle, hXLabel, hYLabel], ...
    'FontName'   , 'AvantGarde' ,...
    'FontSize'   ,   12  );
set( hLegend             , ...
    'FontSize'   , 12    );

hold off 

%print('ComparisonNodesFunction1','-dpng','-r600');


