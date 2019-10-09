% RDisk: spectral interpolation on rhodonea curves
% Example: evaluate black hole data 
% (C) Wolfgang Erb 01.07.2019

% Credit for the black hole data: EHT collaboration et.al.

clear all
close all

% Path
addpath(genpath('./data/'));
addpath(genpath('./core/'));

% Set parameters 

m       = [8,9];      % Frequency parameter of Rhodonea curve

sis     = 'square';     % Form of spectral index set: 'triangle' or 'square'
av      = 0;            % av = 0: no averaging, av = 1 averaging on boundary

% Original image of black hole
% You can download it from https://cdn.eso.org/images/publicationjpg/eso1907a.jpg  
rgbImage = imread('eso1907a.jpg');  
   
redChannel = double(rgbImage(:,:,1)); % Red channel
greenChannel = double(rgbImage(:,:,2)); % Green channel
blueChannel = double(rgbImage(:,:,3)); % Blue channel
grayChannel = double(rgb2gray(rgbImage)); % Gray channel

% Polar coordinates of rhodonea nodes
[rrRD, thRD] = RDpts(m);

center = [2007,1050];
radius = 900;            

xRD = round(center(1) + radius*rrRD.*cos(thRD)); 
yRD = round(center(2) + radius*rrRD.*sin(thRD));
        
% Extraction of function values

idRD = sub2ind(size(redChannel), yRD,xRD);

fred = redChannel(idRD);
fgreen = greenChannel(idRD);
fblue = blueChannel(idRD);
fgray = grayChannel(idRD);
    
% Computation of realvalued coefficient matrix
Gred = RDdatM(m,fred); [~,CRred] = RDcfsfft(m,Gred,sis,av);
Ggreen = RDdatM(m,fgreen); [~,CRgreen] = RDcfsfft(m,Ggreen,sis,av);
Gblue = RDdatM(m,fblue); [~,CRblue] = RDcfsfft(m,Gblue,sis,av);
Ggray = RDdatM(m,fgray); [~,CRgray] = RDcfsfft(m,Ggray,sis,av);

% Initialization for plot
x = linspace(center(1)-radius,center(1)+radius,radius+1);
y = linspace(center(2)-radius,center(2)+radius,radius+1);
[X,Y] = meshgrid(x,y);

Zred = 255*ones(size(X)); Ored = 255*ones(size(X));
Zgreen = 255*ones(size(X)); Ogreen = 255*ones(size(X));
Zblue = 255*ones(size(X)); Oblue = 255*ones(size(X));
Zgray = 255*ones(size(X)); Ogray = 255*ones(size(X));

% Polar coordinates and index sets

[theta,r] = cart2pol((X-center(1))/radius,(Y-center(2))/radius);
idx = r<=1;
idgrid = sub2ind(size(redChannel), Y(idx),X(idx));

% Original image in disk

Ored(idx)=redChannel(idgrid);
Ogreen(idx)=greenChannel(idgrid);
Oblue(idx)=blueChannel(idgrid);
Ogray(idx)=grayChannel(idgrid);

O = cat(3, uint8(Ored), uint8(Ogreen), uint8(Oblue));

% Values of the interpolant at grid
Sfred = RDeval(CRred,m,r(idx)',theta(idx)');
Sfgreen = RDeval(CRgreen,m,r(idx)',theta(idx)');   
Sfblue = RDeval(CRblue,m,r(idx)',theta(idx)');
Sfgray = RDeval(CRgray,m,r(idx)',theta(idx)');

Itfred = RDeval(CRred,m,rrRD',thRD');
Itfgreen = RDeval(CRgreen,m,rrRD',thRD');   
Itfblue = RDeval(CRblue,m,rrRD',thRD'); 
Itfgray = RDeval(CRgray,m,rrRD',thRD');

Zred(idx) = Sfred;
Zgreen(idx) = Sfgreen;
Zblue(idx) = Sfblue;
Zgray(idx) = Sfgray;

Z = cat(3, uint8(Zred), uint8(Zgreen), uint8(Zblue));

maxnorm = max([max(max(abs(Ored))),max(max(abs(Ogreen))),max(max(abs(Oblue)))]);
maxnormred = max(max(abs(Ored)));
maxnormgray = max(max(abs(Ogray)));

% Plot the original function
figure(1)
imagesc(O)
set(gca,'XTick',[],'YTick',[])
axis square
title('Original black hole data from EHT collaboration','FontName','Avantgarde','FontSize',12);
% TT = get(gca,'TightInset');
% set(gca,'Position', [TT(1) TT(2), 1 - TT(1)-TT(3),1-TT(2)-TT(4)]);


% Plot the interpolation function
figure(2)
imagesc(Z)
set(gca,'XTick',[],'YTick',[])
axis square
title('Interpolated data on rhodonea nodes','FontName','Avantgarde','FontSize',12);
% TT = get(gca,'TightInset');
% set(gca,'Position', [TT(1) TT(2), 1 - TT(1)-TT(3),1-TT(2)-TT(4)]);
% print('data/Blackhole0809','-dpng','-r600');

% Plot the interpolation error
figure(3)
imagesc(abs(Ored-Zred)/maxnormred)
set(gca,'XTick',[],'YTick',[])
colormap(hot)
colorbar
axis square
title('Interpolation error (red channel)','FontName','Avantgarde','FontSize',12);
% TT = get(gca,'TightInset');
% set(gca,'Position', [TT(1) TT(2), 1 - TT(1)-TT(3),1-TT(2)-TT(4)]);
% print('data/Blackhole0809error','-dpng','-r600');

% Print interpolation errors
maxEred = max(max(abs(Zred-Ored)));
maxEgreen = max(max(abs(Zgreen-Ogreen)));
maxEblue = max(max(abs(Zblue-Oblue)));
maxEgray = max(max(abs(Zgray-Ogray)))/maxnormgray;
maxError = max([maxEred,maxEgreen,maxEblue])/maxnorm;

maxEredInt = max(abs(fred'-Itfred));
maxEgreenInt = max(abs(fgreen'-Itfgreen));
maxEblueInt = max(abs(fblue'-Itfblue));
maxEgrayInt = max(abs(fgray'-Itfgray))/maxnormgray;
maxErrorInt = max([maxEredInt,maxEgreenInt,maxEblueInt])/maxnorm;
             
fprintf('Main test example for interpolation on rhodonea nodes\n');
fprintf('-----------------------------------------------------------------\n');
fprintf('Maximal error at interpolation points: %16.14f \n',maxErrorInt);
fprintf('Maximal approximation error: %16.14f \n\n',maxError);
        


