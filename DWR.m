%-----------------------------
% Name:  DWR.m 
%
% Modified: 11.21.2011 dssserke & chrisj
%-----------------------------

%-----------------------------
%   date parameters 
%-----------------------------
%values for 20110518 @1951 UTC
%date                = '20110518';
%ceilometer_height   = 0.27432;
%start_rhi_gate      = 1;
%end_rhi_gate        = 136;
%snr_dbz_value       = -20;

% values for 20101215 @2252
date                = '20101215';
ceilometer_height   = 1.28016;
start_rhi_gate      = 1;
end_rhi_gate        = 96;
snr_dbz_value       = -20;

% values for 20101230 @2052
%date                = '20101230';
%ceilometer_height   = 0.335;
%start_rhi_gate      = 1;
%end_rhi_gate        = 90;
%snr_dbz_value       = -20;

% values for 20101231 @0022
%date                = '20101231';
%ceilometer_height   = 0.304;
%start_rhi_gate      = 1;
%end_rhi_gate        = 90;
%snr_dbz_value       = -20;

%-----------------------------
%  data dirs 
%-----------------------------
k_band_dir = '/d1/serke/projects/NIRSS_NASA/data/K_Band_txt/';
%k_band_dir = '/d1/chrisj/NIRSS_Platteville/2012_AMS/NIRSS_radar/';
%s_band_dir = '/d1/chrisj/NIRSS_Platteville/2012_AMS/CHILL_radar/NETCDF/20101231/';
s_band_dir = '/d1/serke/projects/IHL_MIT/data/CHILL@Platte/20101215/';
%s_band_dir = '/d1/serke/projects/IHL_MIT/data/CHILL@Platte/20101231/';
%s_band_dir = '/d1/serke/projects/IHL_MIT/data/CHILL@Platte/20110518/';
%s_band_dir = '/d1/serke/projects/IHL_MIT/data/CHILL@Platte/20101230/';
%s_band_dir = '/d1/chrisj/NIRSS_Platteville/2012_AMS/CHILL_radar/NETCDF/2010_1215/';
%s_band_dir = '/d1/chrisj/NIRSS_Platteville/2012_AMS/CHILL_radar/NETCDF/2010_1230/';
%s_band_dir = '/d1/chrisj/NIRSS_Platteville/2012_AMS/CHILL_radar/NETCDF/2011_0518/';

%-----------------------------
%  read in S band data 
%-----------------------------
%ncload([s_band_dir 'cfrad.20101231_002209_000_v73_RHI.nc']);
ncload([s_band_dir 'cfrad.20101215_225219_000_v1_RHI.nc']);
%ncload([s_band_dir 'cfrad.20101230_205210_000_v2_RHI.nc']);
%ncload([s_band_dir 'cfrad.20110518_195135_000_v5_RHI.nc']);

%-----------------------------
%  remove S-band values near noise threshold  
%-----------------------------
% CHANGED to DZ instead of DBZ in line below 
%if date == '20110518' 
%  DBZ_s        = DZ(start_rhi_gate:end_rhi_gate,200);
%else
%  DBZ_s        = DBZ(start_rhi_gate:end_rhi_gate,200);
%end
DBZ_s = DBZ(start_rhi_gate:end_rhi_gate,200);

ind          = find(DBZ_s/1000 <= snr_dbz_value);
DBZ_s(ind)   = nan;

%-----------------------------
% remove S-band values below ceil hgt 
%-----------------------------
alt_s		    = tan(elevation(start_rhi_gate:end_rhi_gate)*pi/180)*29.5;
%ind_s               = find(alt_s <= ceilometer_height);
%DBZ_s(ind_s)        = nan;

%-----------------------------
%  read in Ka band data 
%-----------------------------
%kaband = load([k_band_dir '20101231_0022.txt']);
kaband = load([k_band_dir '20101215_2252.txt']);
%kaband = load([k_band_dir '20101230_2052.txt']);
%kaband = load([k_band_dir '20110518_1951.txt']);

%-----------------------------
%  mess with Ka band data 
%-----------------------------
alt_k = kaband(:,2)/1000 - 0.4;
DBZ_k = kaband(:,3);

figure
plot(DBZ_s/1000,alt_s,'b*-');
hold on;
plot(DBZ_k,alt_k,'r-');
grid on;

%-----------------------------
% remove ka-band values below ceil hgt 
%-----------------------------
%ind_k               = find(alt_k <= ceilometer_height);
%DBZ_k(ind_k)        = nan;

%-----------------------------
% remove ka-band values above RAYLEIGH limit 
%-----------------------------
DBZ_k_orig          = DBZ_k;
ind_k               = find(DBZ_k >= 25);
DBZ_k(ind_k)        = nan;

%-----------------------------
% find lowest non-nan S value,
% and find diff in first s and K at that height
% then adjust k profile to that difference 
%-----------------------------
DBZ_s_smooth = smooth(DBZ_s,5);
DBZ_k_smooth = smooth(DBZ_k,5);
for kk=1:length(DBZ_k_smooth)
  if DBZ_k_smooth(kk) == 0
    DBZ_k_smooth(kk) = NaN;
  end
end

done = 0;
for j=1:length(DBZ_k_smooth)
  if ~isnan(DBZ_k_smooth(j)) & done == 0
    first_DBZ_k = DBZ_k_smooth(j);
    first_j     = j;
    done        = 1;
  end
end

if date == '20110518'
  first_j = 216;
  first_DBZ_k = DBZ_k_smooth(first_j);
  DBZ_k_smooth(1:100) = NaN;
end

done = 0;
for i=1:length(DBZ_s_smooth)
  if ~isnan(DBZ_s_smooth(i)) & done == 0
    first_DBZ_s = DBZ_s_smooth(i);
    first_i     = i;
    done        = 1;
  end
end

min_alt_diff = min(abs(alt_k(first_j) - alt_s));
ind_s_good   = find(abs(alt_k(first_j) - alt_s) == min_alt_diff);
diff_ks      = first_DBZ_k - DBZ_s_smooth(ind_s_good)/1000;
DBZ_k_match  = DBZ_k_smooth - diff_ks;
DBZ_k_match(1) = NaN;

%-----------------------------
% 
%-----------------------------
ind_hi_s_refl = find(~isnan(DBZ_s_smooth));
ind_hi_s_refl = ind_hi_s_refl(end-1);

for ii = 1:ind_hi_s_refl - first_i
  alt_s_find(ii)  = alt_s(first_i+ii-1);
  dif             = abs(alt_k - alt_s_find(ii));
  ind_alt_k(ii) = find(dif == min(dif));
  dwr_dif(ii)   = DBZ_s_smooth(first_i+ii-1)/1000 - DBZ_k_match(ind_alt_k(ii)); 
  %DBZ_s(first_i)/1000
  %DBZ_k_smooth(ind_alt_k(ii)
end

%-----------------------------
% plotting
%-----------------------------
figure;
plot(DBZ_s_smooth/1000,alt_s,'b*-');
hold on;
plot(DBZ_k_match,alt_k,'r*-');
grid on;
xlabel('Reflectivity [dBZ]');
ylabel('Height [AGL km]');
legend('CHILL S-band','METEK Ka-band');

figure;
if date == '20101215'
  dwr_dif(1) = NaN;
  plot(dwr_dif(1:end-3),alt_s_find(1:end-3),'r*-');
else
  plot(dwr_dif,alt_s_find,'r*-');
end
grid on;

