clear;close all;clc;

domainrbsb = '../inputdata/domain_lnd_columbia_half_hexagon_RBSB_c210406.nc';
mosartrbsb = '../inputdata/MOSART_columbia_half_hexagon_RBSB_c210406.nc';
domaintbsb = '../inputdata/domain_lnd_columbia_half_hexagon_TBSB_c210406.nc';
mosarttbsb = '../inputdata/MOSART_columbia_half_hexagon_TBSB_c210406.nc';
domainsqu = '../inputdata/domain_lnd_columbia_half_square_c201016.nc';
mosartsqu = '../inputdata/MOSART_columbia_half_square_c201016.nc';
% Read data from Hexagon mesh
xv_rbsb = ncread(domainrbsb,'xv');
yv_rbsb = ncread(domainrbsb,'yv');
rslp_rbsb = ncread(mosartrbsb,'rslp');
rlen_rbsb = ncread(mosartrbsb,'rlen');
rdep_rbsb = ncread(mosartrbsb,'rdep');
rwid_rbsb = ncread(mosartrbsb,'rwid');

xv_tbsb = ncread(domaintbsb,'xv');
yv_tbsb = ncread(domaintbsb,'yv');
rslp_tbsb = ncread(mosarttbsb,'rslp');
rlen_tbsb = ncread(mosarttbsb,'rlen');
rdep_tbsb = ncread(mosarttbsb,'rdep');
rwid_tbsb = ncread(mosarttbsb,'rwid');

% Read data from tradition mesh
xv = ncread(domainsqu,'xv');
yv = ncread(domainsqu,'yv');
rslp = ncread(mosartsqu,'rslp');
rlen = ncread(mosartsqu,'rlen');
rdep = ncread(mosartsqu,'rdep');
rwid = ncread(mosartsqu,'rwid');

figure; set(gcf,'Position',[10 10 1000 1000]);

subplot(4,3,1);
patch(xv,yv,rslp); axis equal; colorbar; caxis([0 0.03]);
title('Square','FontSize',15,'FontWeight','bold');
ylabel('river slope [m/m]','FontSize',15,'FontWeight','bold');
subplot(4,3,2);
patch(xv_rbsb,yv_rbsb,rslp_rbsb); axis equal;colorbar;caxis([0 0.03]);
title('RBSB','FontSize',15,'FontWeight','bold');
subplot(4,3,3);
patch(xv_tbsb,yv_tbsb,rslp_tbsb); axis equal;colorbar;caxis([0 0.03]);
title('TBSB','FontSize',15,'FontWeight','bold');

subplot(4,3,4);
patch(xv,yv,rlen); axis equal; colorbar; caxis([0 2e5]);
ylabel('river length [m]','FontSize',15,'FontWeight','bold');
subplot(4,3,5);
patch(xv_rbsb,yv_rbsb,rlen_rbsb); axis equal;colorbar;caxis([0 2e5]);
subplot(4,3,6);
patch(xv_tbsb,yv_tbsb,rlen_tbsb); axis equal;colorbar;caxis([0 2e5]);

subplot(4,3,7);
patch(xv,yv,rwid); axis equal; colorbar; %caxis([0 2e5]);
ylabel('river width [m]','FontSize',15,'FontWeight','bold');
subplot(4,3,8);
patch(xv_rbsb,yv_rbsb,rwid_rbsb); axis equal;colorbar;%caxis([0 2e5]);
subplot(4,3,9);
patch(xv_tbsb,yv_tbsb,rwid_tbsb); axis equal;colorbar;%caxis([0 2e5]);

subplot(4,3,10);
patch(xv,yv,rdep); axis equal; colorbar; %caxis([0 2e5]);
ylabel('river depth [m]','FontSize',15,'FontWeight','bold');
subplot(4,3,11);
patch(xv_rbsb,yv_rbsb,rdep_rbsb); axis equal;colorbar;%caxis([0 2e5]);
subplot(4,3,12);
patch(xv_tbsb,yv_tbsb,rdep_tbsb); axis equal;colorbar;%caxis([0 2e5]);

figure;
nbins1 = 20;
nbins2 = 40;
subplot(1,2,1);
histogram(rlen_rbsb,nbins1); hold on;
histogram(rlen_tbsb,nbins1);
histogram(rlen,nbins2);
title('length','FontSize',15,'FontWeight','bold');
subplot(1,2,2);
histogram(rslp_rbsb,nbins1); hold on;
histogram(rslp_tbsb,nbins1);
histogram(rslp,nbins2);
title('slope','FontSize',15,'FontWeight','bold');


% verify shapefile and netcdf file for slope
rslpb_rbsb_sp = [shaperead('../RBSB/slope_between_polygon.shp').slopb]';
rslpw_rbsb_sp = [shaperead('../RBSB/slope_within_polygon.shp').slopw]';
rslpb_rbsb_nc = ncread('../RBSB/hexwatershed.nc','SlopeB');
rslpw_rbsb_nc = ncread('../RBSB/hexwatershed.nc','SlopeW');
db_rbsb = rslpb_rbsb_sp - rslpb_rbsb_nc;
dw_rbsb = rslpw_rbsb_sp - rslpw_rbsb_nc;
rslpb_tbsb_sp = [shaperead('../TBSB/slope_between_polygon.shp').slopb]';
rslpw_tbsb_sp = [shaperead('../TBSB/slope_within_polygon.shp').slopw]';
rslpb_tbsb_nc = ncread('../TBSB/hexwatershed.nc','SlopeB');
rslpw_tbsb_nc = ncread('../TBSB/hexwatershed.nc','SlopeW');
db_tbsb = rslpb_tbsb_sp - rslpb-tbsb_nc;
dw_tbsb = rslpw_tbsb_sp - rslpw_tbsb_nc;

