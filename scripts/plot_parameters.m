clear;close all;clc;

% Read data from Hexagon mesh
xv_hex = ncread('inputdata/domain_lnd_columbia_half_hexagon_case20210322014_c210324.nc','xv');
yv_hex = ncread('inputdata/domain_lnd_columbia_half_hexagon_case20210322014_c210324.nc','yv');
rslp_hex = ncread('inputdata/MOSART_columbia_half_hexagon_case20210322014_c210324.nc','rslp');
rdep_hex = ncread('inputdata/MOSART_columbia_half_hexagon_case20210322014_c210324.nc','rdep');
rwid_hex = ncread('inputdata/MOSART_columbia_half_hexagon_case20210322014_c210324.nc','rwid');

% Read data from tradition mesh
xv = ncread('../Columbia_River_Basin/domain_lnd_columbia_river_basin_half_c201016.nc','xv');
yv = ncread('../Columbia_River_Basin/domain_lnd_columbia_river_basin_half_c201016.nc','yv');
rslp = ncread('../Columbia_River_Basin/MOSART_columbia_river_basin_half_c201016.nc','rslp');
rdep = ncread('../Columbia_River_Basin/MOSART_columbia_river_basin_half_c201016.nc','rdep');
rwid = ncread('../Columbia_River_Basin/MOSART_columbia_river_basin_half_c201016.nc','rwid');

figure; set(gcf,'Position',[10 10 1400 600]);

subplot(1,2,1);
patch(xv,yv,rslp); axis equal; colorbar; caxis([0 0.03]);
subplot(1,2,2);
patch(xv_hex,yv_hex,rslp_hex); axis equal;colorbar;caxis([0 0.03]);