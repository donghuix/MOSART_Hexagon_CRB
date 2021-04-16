clear;close all;clc

addpath('/Users/xudo627/OneDrive - PNNL/donghui/CODE/Setup-E3SM-Mac/matlab-scripts-to-process-inputs/');
curdir = pwd;
global re; %earth radius [m]
re = 6.37122e6;

check_flow_network = 0;

scenario = 'RBSB';
%mosart_gridded_surfdata_filename = '/Users/xudo627/projects/cesm-inputdata/MOSART_Global_half_20200221.nc';
mosart_gridded_surfdata_filename = '../inputdata/MOSART_columbia_half_square_c201016.nc';
mosart_usrdat_name = ['columbia_half_hexagon_' scenario];
hexfile = ['../' scenario '/hexwatershed.nc'];
x    = ncread(hexfile,'X');
y    = ncread(hexfile,'Y');
ID   = ncread(hexfile,'GlobalID');
dnID = ncread(hexfile,'GlobalID_down');
slopew  = ncread(hexfile,'SlopeW');   % [%]
rlength = ncread(hexfile,'Length');   % [m]
rdep    = ncread(hexfile,'Depth'); rdep = rdep ./ 10^2.4;
rwid    = ncread(hexfile,'Width'); rwid = rwid ./ 10^3.6;

slopew = slopew ./ 100; % convert from [%] to [m/m]

ID2  = NaN(length(ID),1);
area = NaN(length(ID),1);
fdir = ones(length(ID),1);
areatotal = ones(length(ID),1);

for i = 1 : length(ID)
    j = find(ID == dnID(i));
    if isempty(j)
        dnID(i) = -9999;
        fdir(i) = 0; 
    else
        dnID(i) = j;
    end
    ID2(i) = i;
end
ID = ID2;

coord = [x y];
save(fullfile(curdir,'hexagon_xy.dat'),'coord','-ascii');
cd('/Users/xudo627/OneDrive - PNNL/donghui/mylib/p/');
args = ['epsg:5070 epsg:4326' ' ' ...
        fullfile(curdir,'hexagon_xy.dat') ' ' ...
        fullfile(curdir,'hexagon_lonlat.dat')];
cmd = ['/Users/xudo627/anaconda3/bin/python xy2lonlat.py ' args];
[status,cmdout] = system(cmd,'-echo');
cd(curdir);
data = load(fullfile(curdir,'hexagon_lonlat.dat'));
delete(fullfile(curdir,'hexagon_xy.dat'));
delete(fullfile(curdir,'hexagon_lonlat.dat'));
lon = data(:,1);
lat = data(:,2);

xv = NaN(length(x),6);
yv = NaN(length(y),6);
lonv = NaN(length(lon),6);
latv = NaN(length(lat),6);
%r  = 6.206518287578132e+03;
% estimate radius from the mesh nodes
dist = pdist2([x y],[x y]);
dist(dist == 0) = NaN;
dist = min(dist,[],2);
r = mean(dist);
assert(sum(abs(dist - r)) < 1e-3);
r = r/sqrt(3);

for i = 1 : length(x)
    xv(i,1) = x(i) - r;
    xv(i,2) = x(i) - r/2;
    xv(i,3) = x(i) + r/2;
    xv(i,4) = x(i) + r;
    xv(i,5) = x(i) + r/2;
    xv(i,6) = x(i) - r/2;
    
    yv(i,1) = y(i);
    yv(i,2) = y(i) + sqrt(3)/2*r;
    yv(i,3) = y(i) + sqrt(3)/2*r;
    yv(i,4) = y(i);
    yv(i,5) = y(i) - sqrt(3)/2*r;
    yv(i,6) = y(i) - sqrt(3)/2*r;
end

for i = 1 : 6
    coord = [xv(:,i) yv(:,i)];
    save(fullfile(curdir,'hexagon_xy.dat'),'coord','-ascii');
    cd('/Users/xudo627/OneDrive - PNNL/donghui/mylib/p/');
    args = ['epsg:5070 epsg:4326' ' ' ...
            fullfile(curdir,'hexagon_xy.dat') ' ' ...
            fullfile(curdir,'hexagon_lonlat.dat')];
    cmd = ['/Users/xudo627/anaconda3/bin/python xy2lonlat.py ' args];
    [status,cmdout] = system(cmd,'-echo');
    cd(curdir);
    data = load(fullfile(curdir,'hexagon_lonlat.dat'));
    delete(fullfile(curdir,'hexagon_xy.dat'));
    delete(fullfile(curdir,'hexagon_lonlat.dat'));
    lonv(:,i) = data(:,1);
    latv(:,i) = data(:,2);
end

for i = 1 : length(x)
    area(i,1) =  areaint(latv(i,:),lonv(i,:))*4*pi; 
end

aream2 = area.*(re^2);

for i = 1 : length(areatotal)
    j = i; k = 1;
    while dnID(j) ~= -9999 && k < length(areatotal)
        j = find(ID == dnID(j));
        areatotal(j) = areatotal(j) + 1;
    end
end
areatotal = areatotal.*mean(aream2);
ioutlet = find(areatotal == max(areatotal));

out_netcdf_dir = '../inputdata';
if ~exist(out_netcdf_dir,'dir')
    mkdir(out_netcdf_dir);
end

% Assign copied variables
user_defined_vars = struct([]);
user_defined_vars(1).latixy = lat;
user_defined_vars(1).longxy = lon;
user_defined_vars(1).ID     = ID;
user_defined_vars(1).dnID   = dnID;
user_defined_vars(1).fdir   = fdir;
user_defined_vars(1).area   = aream2;
user_defined_vars(1).rslp   = slopew;
user_defined_vars(1).tslp   = slopew;
user_defined_vars(1).rlen   = rlength;
user_defined_vars(1).areaTotal2 = areatotal;
user_defined_vars(1).rdep   = rdep;
user_defined_vars(1).rwid   = rwid;
user_defined_vars(1).rwid0  = 5.*rwid;

fname_out = CreateMOSARTUgridInputForE3SM3(     ...
            mosart_gridded_surfdata_filename,   ...
            out_netcdf_dir, mosart_usrdat_name, ...
            user_defined_vars);

frac = ncread(fname_out,'frac');
mask = zeros(length(frac),1);
mask(frac > 0) = 1;

fname_out2 =  sprintf('%s/domain_lnd_%s_%s.nc',out_netcdf_dir, ...
                      mosart_usrdat_name,datestr(now, 'cyymmdd'));
                  
area2 = generate_lnd_domain(lon,lat,lonv',latv',frac,mask,fname_out2);



% figure;
% plot(lon,lat,'k.'); hold on
% plot(lon(ioutlet),lat(ioutlet),'ro');
% plot(lonv,latv,'r.');

if check_flow_network
    figure;
    for i = 1 : length(dnID)
        if dnID(i) == -9999
            continue;
        else
            j = find(ID == dnID(i));
            plot([x(i) x(j)], [y(i) y(j)], 'b-','LineWidth',1); hold on;
        end
    end
end