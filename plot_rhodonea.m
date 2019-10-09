% RDisk: spectral interpolation on rhodonea curves
% Plot: This program plots a rose or rhodonea curve and its node points 
% (C) Wolfgang Erb 01.09.2018

close all

% Path
addpath(genpath('./core/'));

% Parameter of rhodonea variety

m = [4,4];

g = gcd(m(1),m(2));

sizeB = 100;

% Initialize
t = 0:0.001:2*pi;
x = zeros(2*g,length(t)); y = x;

% Calculate node points
[rrRC, thRC] = RDpts(m);

xRC = rrRC.*cos(thRC);
yRC = rrRC.*sin(thRC);

% Calculate rhodonea trajectories

for i = 1:2*g
  x(i,:) = cos(m(2)*t/g).*cos(m(1)*t/g+(i-1)*pi/m(2));
  y(i,:) = cos(m(2)*t/g).*sin(m(1)*t/g+(i-1)*pi/m(2));
end

% Start plot

figure

for i = 1:2*g
  plot(x(i,:),y(i,:), 'Color' ,[183,207,246]/255,'LineWidth',2);
  hold on
end


plot(xRC,yRC,'o','LineWidth',2,'MarkerSize',6,...
             'MarkerEdgeColor','k','MarkerFaceColor',[65,105,225]/255);

set(gca,'FontSize',16);
axis square;

xlabel('x'); ylabel('y');
title(['Rhodonea variety $\mathcal{R}^{(\underline{\mathbf{m}})}$ and nodes $\mathbf{RD}^{(\underline{\mathbf{m}})}$,$\underline{\mathbf{m}}=(',num2str(m(1)),',',num2str(m(2)),')$'],'interpreter','latex','fontsize',16)

hold off

