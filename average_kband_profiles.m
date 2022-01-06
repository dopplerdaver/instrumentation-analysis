%-----------------------------
% Name:  average_kband_profiles.m 
%
% Modified: 12.28.2011 dssserke
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

%-----------------------------
%  data dirs 
%-----------------------------
k_band_dir = '/d1/serke/projects/NIRSS_NASA/data/K_Band_txt/AMS/';
mat_dir    = '/d1/serke/projects/NIRSS_NASA/code/mat/';
%k_band_dir = '/d1/fripp_d3/nirss/nasa_glen/CLE_2008_2009/2008/vertical_prof/';
 
cd_string_datafiles = ['cd ' k_band_dir];
cd_string_matfiles  = ['cd ' k_band_dir];
eval(cd_string_datafiles);
unix('ls 200* > kband_files.txt');
load('kband_files.txt');

%keyboard;

%-----------------------------
%  read in Ka band data 
%-----------------------------
%for ii=1:2
for ii=1:length(kband_files)

  eval(cd_string_datafiles);
  kband_file  = [num2str(kband_files(ii)) '.txt'];  
  kaband_data = load(kband_file);

  %-----------------------------
  %  mess with Ka band data 
  %-----------------------------
  alt_k = kaband_data(:,2)/1000 - 0.4;
  DBZ_k = kaband_data(:,3);
  
  %-----------------------------
  % remove ka-band values above RAYLEIGH limit 
  %-----------------------------
  DBZ_k_orig          = DBZ_k;
  ind_k               = find(DBZ_k >= 25);
  DBZ_k(ind_k)        = nan;

  DBZ_k_smooth        = smooth(DBZ_k,5);

  %-----------------------------
  %  initial plotting 
  %-----------------------------
  figure
  plot(DBZ_k_smooth,alt_k,'r-');
  grid on;
  hold on;
  plot(DBZ_k_orig,alt_k,'b-');

  %-----------------------------
  %  
  %-----------------------------
  ind_inf                 = find(isinf(DBZ_k_smooth)); 
  DBZ_k_smooth(ind_inf)   = NaN;

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

  done = 0;
  for j=length(DBZ_k_smooth)-10:-1:1
    if ~isnan(DBZ_k_smooth(j)) & done == 0
      last_DBZ_k = DBZ_k_smooth(j);
      last_j     = j;
      done       = 1;
    end
  end
  
  %-----------------------------
  % find mean slope of smoothed Ka-band REFL  
  %-----------------------------
  disp('------------------------------');
  disp(ii);
  k_smooth_refl_last(ii)  = DBZ_k_smooth(last_j)
  k_smooth_refl_first(ii) = max(DBZ_k_smooth)
  k_smooth_alt_last(ii)   = alt_k(last_j)
  ind_max_dbz             = find(DBZ_k_smooth == max(DBZ_k_smooth));
  if size(ind_max_dbz,1) > 1
    ind_max_dbz = max(ind_max_dbz);
  end
  k_smooth_alt_first(ii)  = alt_k(ind_max_dbz)
  alt_diff(ii)            = k_smooth_alt_last(ii) - k_smooth_alt_first(ii)
  refl_diff(ii)           = k_smooth_refl_last(ii) - k_smooth_refl_first(ii)
  slope(ii)               = refl_diff(ii)/alt_diff(ii)

  %-----------------------------
  % cd back to mat code dir  
  %-----------------------------
  eval(cd_string_datafiles);

  %keyboard;

end

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
plot(DBZ_k_match,alt_k,'r*-');
grid on;
xlabel('Reflectivity [dBZ]');
ylabel('Height [AGL km]');
legend('CHILL S-band','METEK Ka-band');

