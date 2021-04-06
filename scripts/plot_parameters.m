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