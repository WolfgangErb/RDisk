% RDisk: spectral interpolation on rhodonea curves
% Plot: 15 basis functions for the triangular and the 
%       rectangular spectral index sets 
% (C) Wolfgang Erb 01.09.2018

close all

% Path
addpath(genpath('./core/'));      

% This program 
      
x = -1:0.005:1;
N = 15;

[X,Y] = meshgrid(x,x);
[theta,r] = cart2pol(X,Y);
idx = r<=1;
z = nan(size(X));
gam1    = [ 0  0  0  0  0  1  1  1  1  2  2  2  3  3  4];
gam2    = [-4 -2  0  2  4 -3 -1  1  3 -2  0  2 -1  1  0];
gamsq1  = [ 0  0  0   1  1  1  2  2  2  3  3  3  4  4  4];
gamsq2  = [-2  0  2  -1  1  3 -2  0  2 -1  1  3 -2  0  2];

y = zeros(length(r(idx)),N);
ysq = zeros(length(r(idx)),N);
D1 = T(4,r(idx)')';
D2 = cos(theta(idx)*[0:4]);
D3 = sin(theta(idx)*[0:4]);
      
for k = 1:N
    if gam2(k) >= 0
       y(:,k) = D1(:,gam1(k)+1).*D2(:,gam2(k)+1);
    else
       y(:,k) = D1(:,gam1(k)+1).*D3(:,1-gam2(k));
    end
    if gamsq2(k) >= 0
       ysq(:,k) = D1(:,gamsq1(k)+1).*D2(:,gamsq2(k)+1);
    else
       ysq(:,k) = D1(:,gamsq1(k)+1).*D3(:,1-gamsq2(k));
    end       
       ysq(:,12) = D1(:,4).*D3(:,4);
end

figure('Units','normalized')
      for k = 1:N
          z(idx) = y(:,k);
          colormap(gray)
          subplot('position',[gam1(k)/5,0.4+gam2(k)/10,0.155,0.155])
          pcolor(x,x,z), shading interp
          set(gca,'XTick',[],'YTick',[])
          axis square
          title(['$\mathrm{X}_{\mathcal{R},(',num2str(gam1(k)),',',num2str(gam2(k)),')}$'],'interpreter','latex','FontSize',12);
      end
      
figure('Units','normalized')
      for k = 1:N
          z(idx) = ysq(:,k);
          colormap(gray)
          subplot('position',[gamsq1(k)/5,0.4+gamsq2(k)/10,0.155,0.155])
          pcolor(x,x,z), shading interp
          set(gca,'XTick',[],'YTick',[])
          axis square
          title(['$\mathrm{X}_{\mathcal{R},(',num2str(gamsq1(k)),',',num2str(gamsq2(k)),')}$'],'interpreter','latex','FontSize',14);
      end