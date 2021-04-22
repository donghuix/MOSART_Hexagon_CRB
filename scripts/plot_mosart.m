clear;close all;clc;

files1 = dir('../output/*Columbia_half_RBSB*.h0*.nc');
files2 = dir('../output/*Columbia_half_TBSB*.h0*.nc');
files3 = dir('../output/*Columbia_half_DRT*.h0*.nc');

fname1 = '../inputdata/MOSART_columbia_half_hexagon_RBSB_c210406.nc';
fname2 = '../inputdata/MOSART_columbia_half_hexagon_TBSB_c210406.nc';
fname3 = '../inputdata/MOSART_columbia_half_square_c201016.nc';

area2 = ncread(fname2,'area');
[mosart1,iout1] = cat_mosart(files1,{'RIVER_DISCHARGE_OVER_LAND_LIQ','RIVER_DISCHARGE_TO_OCEAN_LIQ','QSUR_LIQ','QSUB_LIQ'});
[mosart2,iout2] = cat_mosart(files2,{'RIVER_DISCHARGE_OVER_LAND_LIQ','RIVER_DISCHARGE_TO_OCEAN_LIQ','QSUR_LIQ','QSUB_LIQ'});
[mosart3,iout3] = cat_mosart(files3,{'RIVER_DISCHARGE_OVER_LAND_LIQ','RIVER_DISCHARGE_TO_OCEAN_LIQ','QSUR_LIQ','QSUB_LIQ'});

q1(:,:,1) = mosart1.RIVER_DISCHARGE_OVER_LAND_LIQ;
q1(:,:,2) = mosart1.RIVER_DISCHARGE_TO_OCEAN_LIQ;
q1 = nanmean(q1,3);

q2(:,:,1) = mosart2.RIVER_DISCHARGE_OVER_LAND_LIQ;
q2(:,:,2) = mosart2.RIVER_DISCHARGE_TO_OCEAN_LIQ;
q2 = nanmean(q2,3);

q3(:,:,1) = mosart3.RIVER_DISCHARGE_OVER_LAND_LIQ;
q3(:,:,2) = mosart3.RIVER_DISCHARGE_TO_OCEAN_LIQ;
q3 = nanmean(q3,3);

    
addpath('/Users/xudo627/projects/topotoolbox/colormaps/');


corr = [45.63333333	-121.9544444; ...
        45.65555556	-121.5472222; ...
        45.71472222	-120.6936111; ...
        45.75666667	-121.2088889; ...
        45.9347  	-119.2958;    ...
        44.7289  	-121.2572;    ...	
        44.6039	    -121.2778;    ...
        45.6075	    -121.1722;    ... 
        45.7522	    -121.5258;    ...
        45.24166667	-121.0938889; ...
        46.0278	    -118.7286];
ID = {'BON','HOD','JDA','KLC','MCN', ...
      'PEL','ROU','TDA','WHS','WHT','WWA'}; 

T = readtable('/Users/xudo627/projects/Columbia_River_Basin/BPA_NRNI_flow/NRNI_Flows_1929-2008_Corrected_04-2017.csv','HeaderLines',1);
time = datenum(T.Var2(6:end));
[Yr,Mo,Da] = datevec(time);
ind_obs = find(Yr >=1979 & Yr <=2008); 
Yr = Yr(ind_obs);
Mo = Mo(ind_obs);
Da = Da(ind_obs);

k = 1;
for i = 1979 : 2008
    for j = 1 : 12
        t_sim(k,1) = datenum(i,j,15,1,1,1);
        k = k + 1;
    end
end
[yr_sim,mo_sim] = datevec(t_sim);

isub = 1;
figure(1); set(gcf,'Position',[10 10 1400 400]);
figure(2); set(gcf,'Position',[10 10 1400 400]);
for iGauge = 8%1 : length(ID)
    if isempty(find(strcmp(T.Properties.VariableNames,ID{iGauge})))
        ind(iGauge) = 0;
    else
        tmp = T.(ID{iGauge});
        tmp = tmp(6:end);
        data(1).(ID{iGauge}) = cellfun(@str2num,tmp(ind_obs));
        ind(iGauge) = 1;
        %[ioutlet, icontributing] = find_mosart_cell(fname,corr(iGauge,2),corr(iGauge,1));
    end
    
    if ind(iGauge) == 1
        k = 1;
        for i = 1979 : 2008
            for j = 1 : 12
                ii = find(Yr == i & Mo == j);
                data_mon(k,1) = nanmean(data.(ID{iGauge})(ii))*0.0283168;
                t_mon(k,1) = datenum(i,j,15,1,1,1);
                k = k + 1;
            end
        end
        
        [ioutlet1, icontributing1] = find_mosart_cell(fname1,corr(8,2),corr(8,1));
        [ioutlet2, icontributing2] = find_mosart_cell(fname2,corr(8,2),corr(8,1));
        [ioutlet3, icontributing3] = find_mosart_cell(fname3,corr(8,2),corr(8,1));

        %subplot(6,2,isub);
        figure(1);
        plot(t_mon,data_mon,'k-','LineWidth',2); hold on; grid on;
        plot(t_mon,q1(ioutlet1,:),'b--','LineWidth',1.5);
        plot(t_mon,q2(ioutlet2,:),'r--','LineWidth',1.5);
        plot(t_mon,q3(ioutlet3,:),'g--','LineWidth',1.5);
        xlim([t_mon(1),t_mon(end)]);
%         [~,~,R2] = LSE(data_mon,qall);
%         title([ID{iGauge} ', R2 = ' num2str(R2)],'FontSize',15,'FontWeight','bold');
        xticks([t_mon(1:24:end)]);
        datetick('x','mmmyy','keepticks');
        isub = isub + 1;

        leg = legend('NRNI Flows','RBSB','TBSB','DRT');
        leg.FontSize = 15;
        leg.FontWeight = 'bold';

        set(gca,'FontSize',15);
        
        figure(2);
        plot(1:12,nanmean(reshape(data_mon,[12,30]),2),'k-','LineWidth',2); hold on; grid on;
        plot(1:12,nanmean(reshape(q1(ioutlet1,:)',[12,30]),2),'b--','LineWidth',1.5);
        plot(1:12,nanmean(reshape(q2(ioutlet2,:)',[12,30]),2),'r--','LineWidth',1.5);
        plot(1:12,nanmean(reshape(q3(ioutlet3,:)',[12,30]),2),'g--','LineWidth',1.5);
        xlim([1 12]);
        xticks([1 : 12]);
        xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
        
        leg = legend('NRNI Flows','RBSB','TBSB','DRT');
        leg.FontSize = 15;
        leg.FontWeight = 'bold';

        set(gca,'FontSize',15);
        
    end

end

xc1 = ncread('../inputdata/domain_lnd_columbia_half_hexagon_RBSB_c210406.nc','xc');
yc1 = ncread('../inputdata/domain_lnd_columbia_half_hexagon_RBSB_c210406.nc','yc');
[gsim,gsim2] = search_GSIM_gauge(xc1,yc1,xc1(ioutlet1),yc1(ioutlet1),1);
[lon_gsim,lat_gsim,S,yr_gsim,mo_gsim,da_gsim,mu_gsim,sd_gsim,cv_gsim] = ...
             get_GSIM_discharge(gsim,2);
ind = find(yr_gsim >= 1979 & yr_gsim <= 2008);

[ioutlet1, icontributing1] = find_mosart_cell(fname1,lon_gsim,lat_gsim);
[ioutlet2, icontributing2] = find_mosart_cell(fname2,lon_gsim,lat_gsim);
[ioutlet3, icontributing3] = find_mosart_cell(fname3,lon_gsim,lat_gsim);
        
figure(3); set(gcf,'Position',[10 10 1400 400]);
plot(t_mon,mu_gsim(ind),'k-','LineWidth',2); hold on; grid on;
plot(t_mon,q1(ioutlet1,:),'b--','LineWidth',1.5);
plot(t_mon,q2(ioutlet2,:),'r--','LineWidth',1.5);
plot(t_mon,q3(ioutlet3,:),'g--','LineWidth',1.5);
xlim([t_mon(1),t_mon(end)]);
%         [~,~,R2] = LSE(data_mon,qall);
%         title([ID{iGauge} ', R2 = ' num2str(R2)],'FontSize',15,'FontWeight','bold');
xticks([t_mon(1:24:end)]);
datetick('x','mmmyy','keepticks');
isub = isub + 1;

leg = legend('USGS','RBSB','TBSB','DRT');
leg.FontSize = 15;
leg.FontWeight = 'bold';

set(gca,'FontSize',15);

figure(4)
plot(1:12,nanmean(reshape(mu_gsim(ind),[12,30]),2),'k-','LineWidth',2); hold on; grid on;
plot(1:12,nanmean(reshape(q1(ioutlet1,:)',[12,30]),2),'b--','LineWidth',1.5);
plot(1:12,nanmean(reshape(q2(ioutlet2,:)',[12,30]),2),'r--','LineWidth',1.5);
plot(1:12,nanmean(reshape(q3(ioutlet3,:)',[12,30]),2),'g--','LineWidth',1.5);
xlim([1 12]);
xticks([1 : 12]);
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});

leg = legend('USGS','RBSB','TBSB','DRT');
leg.FontSize = 15;
leg.FontWeight = 'bold';

set(gca,'FontSize',15);
        
% xv  = ncread('inputdata/domain_lnd_columbia_half_hexagon_case20210322014_c210324.nc','xv');
% yv  = ncread('inputdata/domain_lnd_columbia_half_hexagon_case20210322014_c210324.nc','yv');
% xv2  = ncread('/Users/xudo627/projects/Columbia_River_Basin/domain_lnd_columbia_river_basin_half_c201016.nc','xv');
% yv2  = ncread('/Users/xudo627/projects/Columbia_River_Basin/domain_lnd_columbia_river_basin_half_c201016.nc','yv');
% % xv2 = ncread('../UQ_test/inputdata/domain_lnd_columbia_half_c200831.nc','xv');
% % yv2 = ncread('../UQ_test/inputdata/domain_lnd_columbia_half_c200831.nc','yv');
% figure; set(gcf,'Position',[10 10 1400 600]);
% subplot(1,2,1);
% patch(xv,yv,qhex); axis equal; colormap(flowcolor(256));
% caxis([0 3000]);
% % subplot(1,2,2);
% % patch(xv2,yv2,qsqu); axis equal; colormap(flowcolor(256));
% % caxis([0 300]);
% cl = colorbar('east');
% cl.Position = cl.Position + [0.07 0 0 0];
% 
% subplot(1,2,2);
% patch(xv2,yv2,qsqu); axis equal; colormap(flowcolor(256));
% caxis([0 3000]);
% % subplot(1,2,2);
% % patch(xv2,yv2,qsqu); axis equal; colormap(flowcolor(256));
% % caxis([0 300]);
% cl = colorbar('east');
% cl.Position = cl.Position + [0.07 0 0 0];
% 
% % 
% % figure;
% % patch(xv,yv,qhex); axis equal; colormap(flowcolor(256)); hold on;
% % %patch(xv2,yv2,qsqu); axis equal; colormap(flowcolor(256));
% % plot(xc2,yc2,'r.');
% % 
% % figure; set(gcf,'Position',[10 10 1400 600]);
% % subplot(1,2,1);
% % patch(xv,yv,nanmean(rtot,2)); axis equal;
% % caxis([0 100]);
% % subplot(1,2,2);
% % patch(xv2,yv2,nanmean(rtot2,2)); axis equal;
% % caxis([0 100]);
% % cl = colorbar('east');
% % cl.Position = cl.Position + [0.07 0 0 0];
% 
% q1 = qsub+qsur;
% q2 = mosart2.QSUB_LIQ + mosart2.QSUR_LIQ;
% 
% q1 = nanmean(q1); q1 = q1(:);
% q2 = nanmean(q2); q2 = q2(:);
% 
% figure;
% plot(q1(25:end),'r--','LineWidth',2); hold on;
% plot(q2(1:end-24),'b-','LineWidth',2);
% 
% q1 = qsub+qsur;
% q2 = mosart2.QSUB_LIQ + mosart2.QSUR_LIQ;
% 
% q1 = nanmean(q1(:,25:end),2);
% q2 = nanmean(q2(:,:,1 : end-24),3);
% 
% figure;set(gcf,'Position',[10 10 1400 600]);
% subplot(1,2,1);
% patch(xv,yv,q1); axis equal; colormap(flowcolor(256)); colorbar; caxis([0 120])
% subplot(1,2,2);
% patch(xv2,yv2,q2); axis equal; colormap(flowcolor(256)); colorbar; caxis([0 120])

