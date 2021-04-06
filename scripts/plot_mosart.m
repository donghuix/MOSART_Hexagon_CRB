clear;close all;clc;

files = dir('./output/*Columbia_half_Hexagon*.h0*.nc');
files2 = dir('/Users/xudo627/projects/Columbia_River_Basin/output/half/latlon/*.nc');
fname2 = '/Users/xudo627/projects/Columbia_River_Basin/MOSART_columbia_river_basin_half_c201016.nc';
area2 = ncread(fname2,'area');
[mosart2,iout2] = cat_mosart(files2,{'RIVER_DISCHARGE_OVER_LAND_LIQ','RIVER_DISCHARGE_TO_OCEAN_LIQ','QSUR_LIQ','QSUB_LIQ'});
qsqu(:,:,1) = mosart2.RIVER_DISCHARGE_OVER_LAND_LIQ;
qsqu(:,:,2) = mosart2.RIVER_DISCHARGE_TO_OCEAN_LIQ;
qsqu = nanmean(qsqu,3);
qsqu = mean(qsqu,2);
    
addpath('/Users/xudo627/projects/topotoolbox/colormaps/');

fname = 'inputdata/MOSART_columbia_half_hexagon_case20210322014_c210324.nc';
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
  
[ioutlet2, icontributing2] = find_mosart_cell(fname2,corr(8,2),corr(8,1));

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
figure; set(gcf,'Position',[10 10 1400 400]);
for iGauge = 8%1 : length(ID)
    if isempty(find(strcmp(T.Properties.VariableNames,ID{iGauge})))
        ind(iGauge) = 0;
    else
        tmp = T.(ID{iGauge});
        tmp = tmp(6:end);
        data(1).(ID{iGauge}) = cellfun(@str2num,tmp(ind_obs));
        ind(iGauge) = 1;
        [ioutlet, icontributing] = find_mosart_cell(fname,corr(iGauge,2),corr(iGauge,1));
    end
    
    if ind(iGauge) == 1
    for i = 1 : length(files)
        filename = fullfile(files(i).folder,files(i).name);
        if ( i == 1 ) 
            areatotal = ncread(filename,'areatotal2');
            area = ncread(filename,'area');
            xc   = ncread(filename,'lon');
            yc   = ncread(filename,'lat');
%             fprintf(['Hexgon mesh contributing area is ' num2str(areatotal(ioutlet)/1e6) ' km^{2}\n']);
%             fprintf(['xc = ' num2str(xc(ioutlet)) ', yc = ' num2str(yc(ioutlet)) '\n']);
        end
        rdo(:,i) = ncread(filename,'RIVER_DISCHARGE_TO_OCEAN_LIQ');
        qr(i) = rdo(ioutlet,i);%*35.3146667;
        rdl(:,i) = ncread(filename,'RIVER_DISCHARGE_OVER_LAND_LIQ');
        qsur(:,i) = ncread(filename,'QSUR_LIQ');
        qsub(:,i) = ncread(filename,'QSUB_LIQ');
        
        ql(i) = rdl(ioutlet,i);%*35.3146667; 
        qall(i,1) = rdl(ioutlet,i);
        qall(i,2) = rdo(ioutlet,i);
    end
    qhex(:,:,1) = rdo;
    qhex(:,:,2) = rdl;
    qhex = nanmean(qhex,3);
    qhex = mean(qhex,2);
    qall = nanmean(qall,2);
    
    k = 1;
    for i = 1979 : 2008
        for j = 1 : 12
            ii = find(Yr == i & Mo == j);
            data_mon(k,1) = nanmean(data.(ID{iGauge})(ii))*0.0283168;
            t_mon(k,1) = datenum(i,j,15,1,1,1);
            k = k + 1;
        end
    end
    
    %subplot(6,2,isub);
    plot(t_mon,data_mon,'k-','LineWidth',2); hold on; grid on;
    plot(t_mon,qall,'b--','LineWidth',1.5);
    plot(t_mon(25:end),mosart2.RIVER_DISCHARGE_OVER_LAND_LIQ(ioutlet2,1:end-24),'r--','LineWidth',2);
    xlim([t_mon(1),t_mon(end)]);
    [~,~,R2] = LSE(data_mon,qall);
    title([ID{iGauge} ', R2 = ' num2str(R2)],'FontSize',15,'FontWeight','bold');
    xticks([t_mon(1:24:end)]);
    datetick('x','mmmyy','keepticks');
    isub = isub + 1;
    
    leg = legend('NRNI Flows','MOSART-Hexagon','MOSART-Latlon');
    leg.FontSize = 15;
    leg.FontWeight = 'bold';
    
    set(gca,'FontSize',15);
    
    end

end
%     [gsim,gsim2] = search_GSIM_gauge(xc,yc,xc(ioutlet),yc(ioutlet),1);
%     [lon_gsim,lat_gsim,S,yr_gsim,mo_gsim,da_gsim,mu_gsim,sd_gsim,cv_gsim] = ...
%             get_GSIM_discharge(gsim2{2},2);
%     if min(yr_gsim) <= min(yr_sim)
%         i1_gsim = min(find(yr_gsim == min(yr_sim)));
%         i1_sim  = 1;
%     else
%         i1_gsim = 1;
%         i1_sim  = min(find(yr_sim == min(yr_gsim)));
%     end
%     if max(yr_gsim) >= max(yr_sim)
%         i2_gsim = max(find(yr_gsim == max(yr_sim)));
%         i2_sim  = length(yr_sim);
%     else
%        i2_gsim = length(yr_gsim);
%        i2_sim  = max(find(yr_sim == max(yr_gsim)));
%     end
xv  = ncread('inputdata/domain_lnd_columbia_half_hexagon_case20210322014_c210324.nc','xv');
yv  = ncread('inputdata/domain_lnd_columbia_half_hexagon_case20210322014_c210324.nc','yv');
xv2  = ncread('/Users/xudo627/projects/Columbia_River_Basin/domain_lnd_columbia_river_basin_half_c201016.nc','xv');
yv2  = ncread('/Users/xudo627/projects/Columbia_River_Basin/domain_lnd_columbia_river_basin_half_c201016.nc','yv');
% xv2 = ncread('../UQ_test/inputdata/domain_lnd_columbia_half_c200831.nc','xv');
% yv2 = ncread('../UQ_test/inputdata/domain_lnd_columbia_half_c200831.nc','yv');
figure; set(gcf,'Position',[10 10 1400 600]);
subplot(1,2,1);
patch(xv,yv,qhex); axis equal; colormap(flowcolor(256));
caxis([0 3000]);
% subplot(1,2,2);
% patch(xv2,yv2,qsqu); axis equal; colormap(flowcolor(256));
% caxis([0 300]);
cl = colorbar('east');
cl.Position = cl.Position + [0.07 0 0 0];

subplot(1,2,2);
patch(xv2,yv2,qsqu); axis equal; colormap(flowcolor(256));
caxis([0 3000]);
% subplot(1,2,2);
% patch(xv2,yv2,qsqu); axis equal; colormap(flowcolor(256));
% caxis([0 300]);
cl = colorbar('east');
cl.Position = cl.Position + [0.07 0 0 0];

% 
% figure;
% patch(xv,yv,qhex); axis equal; colormap(flowcolor(256)); hold on;
% %patch(xv2,yv2,qsqu); axis equal; colormap(flowcolor(256));
% plot(xc2,yc2,'r.');
% 
% figure; set(gcf,'Position',[10 10 1400 600]);
% subplot(1,2,1);
% patch(xv,yv,nanmean(rtot,2)); axis equal;
% caxis([0 100]);
% subplot(1,2,2);
% patch(xv2,yv2,nanmean(rtot2,2)); axis equal;
% caxis([0 100]);
% cl = colorbar('east');
% cl.Position = cl.Position + [0.07 0 0 0];

q1 = qsub+qsur;
q2 = mosart2.QSUB_LIQ + mosart2.QSUR_LIQ;

q1 = nanmean(q1); q1 = q1(:);
q2 = nanmean(q2); q2 = q2(:);

figure;
plot(q1(25:end),'r--','LineWidth',2); hold on;
plot(q2(1:end-24),'b-','LineWidth',2);

q1 = qsub+qsur;
q2 = mosart2.QSUB_LIQ + mosart2.QSUR_LIQ;

q1 = nanmean(q1(:,25:end),2);
q2 = nanmean(q2(:,:,1 : end-24),3);

figure;set(gcf,'Position',[10 10 1400 600]);
subplot(1,2,1);
patch(xv,yv,q1); axis equal; colormap(flowcolor(256)); colorbar; caxis([0 120])
subplot(1,2,2);
patch(xv2,yv2,q2); axis equal; colormap(flowcolor(256)); colorbar; caxis([0 120])

