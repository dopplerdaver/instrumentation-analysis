%----------------------------------
% Name:      loadplot_slwsonde.m
%
% Purpose:   load slwsonde .csv file
%            load corresponding radiometer file
%            sort fields
%            compute slwsonde LWC 
%            plot required fields
%
% Notes:     for csvread to work properly,
%            need to remove text header
%            and change blank values to NaN 
%            and remove 'Infinities'
%            and remove ':' in time string
%            and remove '-' in date string
%
% Usage:     'loadplot_slwsonde', as long as all internal file  and dir paths are set
%
% Input SLW-sonde files for SNOWIE:
%            L1data = all pre-release data lines removed 
%            L2data = all matlab non-readable characters removed, as described in 'Notes' above
%            L3data = all resonant freqs (values < 10 Hz) removed
%
% Functions 
%      used: csvread, smooth, find
%
% Created:   3.1.2012 - dserke
% Modified:  2.1.2017 - dserke - for SNOWIE sonde flights
%
%----------------------------------

addpath /d1/serke/projects/SLWprobe_NASA/code

%----------------------------------
% set parameters 
%----------------------------------
% hardcoded constants
b0                       = 44.366 % from Hill '94 JOAT
wire_diam                = 0.061  % from Bognar '11
lwc_smooth_points        = 5;
freq_smooth_points       = 11;

%% 2014 OWLeS SLWsonde flights
%project_data_dir         = '/d1/serke/projects/SLWprobe_NASA/data/OWLeS_2014/';
%date_str                 = '20140109/ow002_20140109';
%flight_portion           = 'asc';  
%date_flight_num          = 2;     % can be 1 through 4 
%drop_diam                = 10; 
%%drop_diam                = 10; % used for flt 3 
%%drop_diam                = 7.5; % used for flt 4 
%radiom_exist             = 0;

%% UNCOMMENT LINES BELOW FOR 2014 FLIGHTS AT NASA GRC 
%project_data_dir         = '/d1/serke/projects/SLWprobe_NASA/data/NASA_early2014/';
%date_str                 = '20140301';
%flight_portion           = 'asc';  
%date_flight_num          = 1;     % can be 1 through 3 
%drop_diam                = 50; 
%radiom_exist             = 0;

%% UNCOMMENT LINES BELOW FOR 2012 FLIGHTS AT RADIOMETRICS
%project_data_dir         = '/d1/serke/projects/SLWprobe_NASA/data/SLW_@Radiometrics_2012/';
%date_str                 = '20120307';
%date_flight_num          = 1;     % can be 1 or 2 
%flight_portion           = 'asc'; % can be 'asc' or 'dsc' 
%drop_diam                = 54
%radiom_exist             = 1;

%% UNCOMMENT LINES BELOW FOR 2017 FLIGHTS FOR SNOWIE @ BOISE, ID, 2017.01.18 
%project_data_dir         = '/d1/serke/projects/case_studies/SNOWIE/sonde_files_for_EOL/renamed_for_EOL_files/L3data/';
%date_str                 = '20170118';
%flight_portion           = 'asc';  
%date_flight_num          = 1;     % can be 1 through 3 
%drop_diam                = 50; 
%radiom_exist             = 0;

%% UNCOMMENT LINES BELOW FOR 2017 FLIGHTS FOR SNOWIE @ BOISE, ID, 2017.01.19 
%project_data_dir         = '/d1/serke/projects/case_studies/SNOWIE/sonde_files_for_EOL/renamed_for_EOL_files/L3data/';
%date_str                 = '20170119';
%flight_portion           = 'asc';  
%date_flight_num          = 1;     % can be 1 through 3 
%drop_diam                = 50; 
%radiom_exist             = 0;

% UNCOMMENT LINES BELOW FOR 2017 FLIGHTS FOR SNOWIE @ BOISE, ID, 2017.01.22 
project_data_dir         = '/d2/d1/serke/projects/case_studies/SNOWIE/sonde_files_for_EOL/renamed_for_EOL_files/L3data/';
%date_str                 = '20170118';
%date_str                 = '20170119';  % 3 sondes this date (2 IOP5, 1 IOP6)
date_str                 = '20170122';
%date_str                 = '20170307';
flight_portion           = 'asc';  
date_flight_num          = 3;     % can be 1 through 3 
drop_diam                = 50; 
radiom_exist             = 0;

%----------------------------------
% define paths 
%----------------------------------
if str2num(date_str(1:8)) == 20140109

  if date_flight_num == 1
    prof_start_line    = 1;
    slwsonde_data_file = 'ow001_20140109.noheader.csv'

    if flight_portion == 'asc'
      start                     = 1225;
      stop                      = 1487;
      radiom_data_file          = 'NA';
      title_str                 = 'OWLeS ascent #1, ????Z, January 9, 2014';
    end
  elseif date_flight_num == 2
    prof_start_line    = 1;
    slwsonde_data_file = 'ow002_20140109.noheader.csv'

    if flight_portion == 'asc'
      start                     = 408;
      stop                      = 747;
      radiom_data_file          = 'NA';
      title_str                 = 'OWLeS ascent #2, ????Z, January 9, 2014';
    end
  elseif date_flight_num == 3
    prof_start_line    = 1;
    slwsonde_data_file = 'ow003_20140109.noheader.csv'

    if flight_portion == 'asc'
      start                     = 81;
      stop                      = 646;
      radiom_data_file          = 'NA';
      title_str                 = 'OWLeS ascent #3, ????Z, January 9, 2014';
    end
  elseif date_flight_num == 4
    prof_start_line    = 1;
    slwsonde_data_file = 'ow004_20140109.noheader.csv'

    if flight_portion == 'asc'
      start                     = 189;
      stop                      = 1066;
      radiom_data_file          = 'NA';
      title_str                 = 'OWLeS ascent #4, ????Z, January 9, 2014';
    end
  end     % end of if date_flight_num=1

elseif str2num(date_str) == 20140301
  if date_flight_num == 1
    prof_start_line    = 1;
    slwsonde_data_file = 'ma001_20140301.nohdr.csv'
    if flight_portion == 'asc'
      start                     = 191;
      stop                      = 1066;
      radiom_data_file          = 'NA';
      title_str                 = 'NASA GRC ascent #1, 20:14Z, March 1, 2014';
    end
  end

elseif str2num(date_str) == 20120307

  if date_flight_num == 1
    prof_start_line    = 1;
    slwsonde_data_file = 'sw002_20120307.noheader.csv';

    if flight_portion == 'asc'
      start                     = 420;
      stop                      = 1046;
      radiom_data_file          = '20120307_1408';
      radiom_alt_data_file_path = ([project_data_dir, date_str, '/radiometers/3107/', radiom_data_file, '_alt.txt']);
      %radiom_LWC_data_file_path = ([project_data_dir, date_str, '/radiometers/3107/', radiom_data_file, '_LWCprof.txt']);
      radiom_LWC_data_file_path = ([project_data_dir, date_str, '/radiometers/3107/', radiom_data_file, '_adjustLWCprof.txt']);
      title_str                 = 'ascent #1, 14:07Z, March 7, 2012';
    elseif flight_portion == 'dsc'
      start                     = 1046;
      stop                      = 1415;
      radiom_data_file          = '20120307_1418';
      radiom_alt_data_file_path = ([project_data_dir, date_str, '/radiometers/3058/', radiom_data_file, '_alt.txt']);
      radiom_LWC_data_file_path = ([project_data_dir, date_str, '/radiometers/3058/', radiom_data_file, '_LWCprof.txt']);
      title_str                 = 'descent #1, 14:18Z, March 7, 2012';
    end

  elseif date_flight_num == 2
    prof_start_line    = 1;
    slwsonde_data_file = 'sw003_20120307.noheader.csv';

    if flight_portion == 'asc'
      start                     = 90;
      stop                      = 590;
      radiom_data_file          = '20120307_1422';
      radiom_alt_data_file_path = ([project_data_dir, date_str, '/radiometers/3107/', radiom_data_file, '_alt.txt']);
      %radiom_LWC_data_file_path = ([project_data_dir, date_str, '/radiometers/3107/', radiom_data_file, '_LWCprof.txt']);
      radiom_LWC_data_file_path = ([project_data_dir, date_str, '/radiometers/3107/', radiom_data_file, '_adjustLWCprof.txt']);
      title_str                 = 'ascent #2, 15:57Z, March 7, 2012';
    elseif flight_portion == 'dsc'
      start                     = 590;
      stop                      = 1175;
      radiom_data_file          = '20120307_1605';
      radiom_alt_data_file_path = ([project_data_dir, date_str, '/radiometers/3058/', radiom_data_file, '_alt.txt']);
      radiom_LWC_data_file_path = ([project_data_dir, date_str, '/radiometers/3058/', radiom_data_file, '_LWCprof.txt']);
      title_str                 = 'descent #2, 16:05Z, March 7, 2012';
    end
  end

elseif str2num(date_str) == 20120228
  slwsonde_data_file = 'sw001_20120228_2.noheader.csv';
  prof_start_line    = 8500;

elseif str2num(date_str) == 20170118

  if date_flight_num == 1
    prof_start_line    = 1;
    slwsonde_data_file = 'upperair.NCAR_SLW_sonde.201701182313_Horseshoe_L3data.csv';

    if flight_portion == 'asc'
      start                     = 1;
      stop                      = 68;
      %radiom_data_file          = '20120307_1408';
      %radiom_alt_data_file_path = ([project_data_dir, date_str, '/radiometers/3107/', radiom_data_file, '_alt.txt']);
      title_str                 = 'ascent #1, 23:13Z, January 18, 2017';
    end
  end

elseif str2num(date_str) == 20170119

  if date_flight_num == 1
    prof_start_line    = 1;
    slwsonde_data_file = 'upperair.NCAR_SLW_sonde.201701191528_Horseshoe_L3data.csv';

    if flight_portion == 'asc'
      start                     = 1;
      stop                      = 51;
      %radiom_data_file          = '20120307_1408';
      %radiom_alt_data_file_path = ([project_data_dir, date_str, '/radiometers/3107/', radiom_data_file, '_alt.txt']);
      title_str                 = 'ascent #1, 15:28Z, January 19, 2017';
    end

  elseif date_flight_num == 2
    prof_start_line    = 1;
    slwsonde_data_file = 'upperair.NCAR_SLW_sonde.201701191706_Horseshoe_L3data.csv';

    if flight_portion == 'asc'
      start                     = 1;
      stop                      = 81;
      %radiom_data_file          = '20120307_1408';
      %radiom_alt_data_file_path = ([project_data_dir, date_str, '/radiometers/3107/', radiom_data_file, '_alt.txt']);
      title_str                 = 'ascent #1, 17:06Z, January 19, 2017';
    end

  elseif date_flight_num == 3
    prof_start_line    = 1;
    slwsonde_data_file = 'upperair.NCAR_SLW_sonde.201701192345_Horseshoe_L3data.csv';

    if flight_portion == 'asc'
      start                     = 1;
      stop                      = 114;
      %radiom_data_file          = '20120307_1408';
      %radiom_alt_data_file_path = ([project_data_dir, date_str, '/radiometers/3107/', radiom_data_file, '_alt.txt']);
      title_str                 = 'ascent #1, 23:45Z, January 19, 2017';
    end
  end

elseif str2num(date_str) == 20170122

  if date_flight_num == 1
    prof_start_line    = 1;
    slwsonde_data_file = 'upperair.NCAR_SLW_sonde.201701222118_Horseshoe_L3data.csv';

    if flight_portion == 'asc'
      start                     = 1;
      stop                      = 44;
      %radiom_data_file          = '20120307_1408';
      %radiom_alt_data_file_path = ([project_data_dir, date_str, '/radiometers/3107/', radiom_data_file, '_alt.txt']);
      title_str                 = 'ascent #1, 21:18Z, January 22, 2017';
    end

  elseif date_flight_num == 2
    prof_start_line    = 1;
    slwsonde_data_file = 'upperair.NCAR_SLW_sonde.201701222203_Horseshoe_L3data.csv';

    if flight_portion == 'asc'
      start                     = 1;
      stop                      = 115;
      %radiom_data_file          = '20120307_1422';
      %radiom_alt_data_file_path = ([project_data_dir, date_str, '/radiometers/3107/', radiom_data_file, '_alt.txt']);
      title_str                 = 'ascent #2, 22:03Z, January 22, 2017';
    end
  elseif date_flight_num == 3
    prof_start_line    = 1;
    slwsonde_data_file = 'upperair.NCAR_SLW_sonde.201701222257_Horseshoe_L3data.csv';

    if flight_portion == 'asc'
      start                     = 1;
      stop                      = 121;
      %radiom_data_file          = '20120307_1422';
      %radiom_alt_data_file_path = ([project_data_dir, date_str, '/radiometers/3107/', radiom_data_file, '_alt.txt']);
      title_str                 = 'ascent #3, 22:57Z, January 22, 2017';
    end
  end
  
elseif str2num(date_str) == 20170307
  disp('hello')
  if date_flight_num == 1
    prof_start_line    = 1;
    slwsonde_data_file = 'upperair.NCAR_SLW_sonde.201703071532_Horseshoe_L3data.csv';

    if flight_portion == 'asc'
      start                     = 1;
      stop                      = 208;
      %radiom_data_file          = '20120307_1408';
      %radiom_alt_data_file_path = ([project_data_dir, date_str, '/radiometers/3107/', radiom_data_file, '_alt.txt']);
      title_str                 = 'ascent #1, 18:32Z, March 7, 2017';
    end
  
  end
  
else
  disp('flight date does not exist');
end

%----------------------------------
% open desired data files 
%----------------------------------
%slwsonde_data_file_path = ([project_data_dir, date_str, '/', slwsonde_data_file]);
%slwsonde_data           = csvread(slwsonde_data_file_path);

slwsonde_data_file_path = ([project_data_dir, slwsonde_data_file]);

fid = fopen(slwsonde_data_file_path);
slwsonde_data_cell      = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','delimiter',',');
slwsonde_data           = cell2mat(slwsonde_data_cell);
fclose(fid);

if radiom_exist == 1
  radiom_alt_data         = csvread(radiom_alt_data_file_path);
  radiom_LWC_data         = csvread(radiom_LWC_data_file_path);
end

%----------------------------------
% sort fields from read file 
%----------------------------------
if str2num(date_str) == 20140301
  secs_pmid    = slwsonde_data(prof_start_line:end,5+3);
  alt          = slwsonde_data(prof_start_line:end,7+3);
  pres         = slwsonde_data(prof_start_line:end,8+3);
  temp_cor_c   = slwsonde_data(prof_start_line:end,9+3);
  temp_cor_f   = (temp_cor_c*9/5) + 32;
  humidity     = slwsonde_data(prof_start_line:end,11+3);
  temp_frost_c = slwsonde_data(prof_start_line:end,12+3);
  theta        = slwsonde_data(prof_start_line:end,15+3);
  ascent_rate  = slwsonde_data(prof_start_line:end,18+3);
  freq         = slwsonde_data(prof_start_line:end,42+3);
elseif str2num(date_str) == 20170118 || str2num(date_str) == 20170119 || str2num(date_str) == 20170122 || str2num(date_str) == 20170307
  secs_pmid    = slwsonde_data(prof_start_line:end,4);
  alt          = slwsonde_data(prof_start_line:end,6);
  pres         = slwsonde_data(prof_start_line:end,7);
  temp_cor_c   = slwsonde_data(prof_start_line:end,9);
  temp_cor_f   = (temp_cor_c*9/5) + 32;
  humidity     = slwsonde_data(prof_start_line:end,10);
  temp_frost_c = slwsonde_data(prof_start_line:end,11);
  theta        = slwsonde_data(prof_start_line:end,14);
  ascent_rate  = slwsonde_data(prof_start_line:end,17);
  freq         = slwsonde_data(prof_start_line:end,39);
else
  secs_pmid    = slwsonde_data(prof_start_line:end,5);
  alt          = slwsonde_data(prof_start_line:end,7);
  pres         = slwsonde_data(prof_start_line:end,8);
  temp_cor_c   = slwsonde_data(prof_start_line:end,9);
  temp_cor_f   = (temp_cor_c*9/5) + 32;
  humidity     = slwsonde_data(prof_start_line:end,11);
  temp_frost_c = slwsonde_data(prof_start_line:end,12);
  theta        = slwsonde_data(prof_start_line:end,15);
  ascent_rate  = slwsonde_data(prof_start_line:end,18);
  freq         = slwsonde_data(prof_start_line:end,42);
end

%---------------------------------
% crop relevant fields to desired flight portion 
%---------------------------------
secs_new        = secs_pmid;
alt_new         = alt;
freq_new        = freq;
ascent_rate_new = ascent_rate;
freq0           = freq_new(1);
%secs_new        = secs_pmid(start:stop);
%alt_new         = alt(start:stop);
%freq_new        = freq(start:stop);
%ascent_rate_new = ascent_rate(start:stop);
%freq0           = freq_new(1);

%---------------------------------
% find time diff between each valid (freq>10Hz)  freq value 
%---------------------------------
count             = 0;
secs_since_new(1) = 0;
%index_good_freq   = find(freq_new > 10);

for ii=2:length(secs_new) 
  freq_diff(ii) = freq_new(ii) - freq_new(ii-1);

  % same as the last value
  if freq_diff(ii) == 0
    count              = count + 1;

  % different from the last value
  else
    count = 0;

  end 

  secs_since_new(ii) = count;

end

ind      = find(secs_since_new == 0);

%keyboard;

freq_new        = freq_new(ind);
alt_new         = alt_new(ind);
secs_new        = secs_new(ind);
ascent_rate_new = ascent_rate_new(ind);

%---------------------------------
% smooth freq timeseries 
%---------------------------------
freq_new_smoo = smooth(freq_new,freq_smooth_points);

%---------------------------------
% calculate d(Hz)/dt in minutes  point-to-point 
%---------------------------------
for jj=1:length(freq_new_smoo) - 1
  hzpermin_smoo(jj) = -(freq_new_smoo(jj)-freq_new_smoo(jj+1))*...
                                   60/(secs_new(jj+1)-secs_new(jj)); 
end

%---------------------------------
% lookup table for collect_effic versus drop_diam  
%---------------------------------
if drop_diam < 3
  collect_effic = 0.03
elseif drop_diam >= 3 & drop_diam < 4
  collect_effic = 0.16
elseif drop_diam >= 4 & drop_diam < 5
  collect_effic = 0.32
elseif drop_diam >= 5 & drop_diam < 6
  collect_effic = 0.43
elseif drop_diam >= 6 & drop_diam < 7
  collect_effic = 0.52
elseif drop_diam >= 7 & drop_diam < 8
  collect_effic = 0.59
elseif drop_diam >= 8 & drop_diam < 9
  collect_effic = 0.65
elseif drop_diam >= 9 & drop_diam < 10 
  collect_effic = 0.70
elseif drop_diam >= 10 & drop_diam < 11 
  collect_effic = 0.74
elseif drop_diam >= 11 & drop_diam < 12 
  collect_effic = 0.77
elseif drop_diam >= 12 & drop_diam < 13 
  collect_effic = 0.80
elseif drop_diam >= 13 & drop_diam < 13.5 
  collect_effic = 0.82
elseif drop_diam >= 13.5 & drop_diam < 14.0 
  collect_effic = 0.84
elseif drop_diam >= 14.0 & drop_diam < 15.0 
  collect_effic = 0.86
elseif drop_diam >= 15.0 & drop_diam < 16.0 
  collect_effic = 0.87
elseif drop_diam >= 16.0 & drop_diam < 17.0 
  collect_effic = 0.88
elseif drop_diam >= 17.0 & drop_diam < 18.0 
  collect_effic = 0.89
elseif drop_diam >= 18.0 & drop_diam < 19.0 
  collect_effic = 0.90
elseif drop_diam >= 19.0 & drop_diam < 20.0 
  collect_effic = 0.91
elseif drop_diam >= 20.0 & drop_diam < 21.0 
  collect_effic = 0.92
elseif drop_diam >= 21.0 & drop_diam < 22.0 
  collect_effic = 0.93
elseif drop_diam >= 22.0 & drop_diam < 25.0 
  collect_effic = 0.94
elseif drop_diam >= 25.0 & drop_diam < 30.0 
  collect_effic = 0.95
elseif drop_diam >= 30.0 & drop_diam < 35.0 
  collect_effic = 0.96
elseif drop_diam >= 35.0 & drop_diam < 40.0 
  collect_effic = 0.97
elseif drop_diam >= 40.0 & drop_diam < 45.0 
  collect_effic = 0.98
elseif drop_diam >= 45.0 & drop_diam < 55.0 
  collect_effic = 0.99
else 
  collect_effic = 1.00
end

%---------------------------------
% calc LWC for HILL '94 equation #7 
%---------------------------------
for nn=1:length(hzpermin_smoo)

  LWC_HILL94(nn) = -2*100*100*0.0166667*hzpermin_smoo(nn)*freq0^2/(collect_effic*b0*wire_diam*ascent_rate_new(nn)*freq_new(nn)^3);

end

%---------------------------------
% use fit lines for hzpermin to derive LWC 
%---------------------------------
for kk=1:length(hzpermin_smoo)

  % use polyfit line at Anasphere tunnel test values and lower
  if hzpermin_smoo(kk) > -0.38 & hzpermin_smoo(kk) < 0
    lwc(kk) = 1.1032*abs(hzpermin_smoo(kk)).^2 + 1.0275*abs(hzpermin_smoo(kk)); 

  % linear fit line above Anasphere tunnel test values
  elseif hzpermin_smoo(kk) <= -0.38 
    % linear fit line
    lwc(kk) = 1.6917*abs(hzpermin_smoo(kk)) - 0.0952; 

  else
    lwc (kk) = 0;
  end

end

%---------------------------------
% smooth derived slwsonde LWC
%---------------------------------
LWC_smoo              = smooth(LWC_HILL94,lwc_smooth_points);
if flight_portion == 'asc'
  ind_neg_LWC           = find(LWC_smoo < 0); 
  LWC_smoo(ind_neg_LWC) = 0;
  LWC_smoo(end-3:end)   = 0;
  LWC_smoo(1:7)         = 0;
else
  ind_neg_LWC           = find(LWC_smoo > 0); 
  LWC_smoo(ind_neg_LWC) = 0;
  LWC_smoo(end-3:end)   = 0;
  LWC_smoo(1:7)         = 0;
  LWC_smoo              = abs(LWC_smoo);
end

%lwc_smoo              = smooth(lwc,lwc_smooth_points);
%ind_neg_lwc           = find(lwc_smoo < 0); 
%lwc_smoo(ind_neg_lwc) = 0;
%lwc_smoo(end-3:end)   = 0;
%lwc_smoo(1:3)         = 0;

%---------------------------------
% calculate slwsonde ILW
%---------------------------------
sonde_ILW = 0;
for mm = 1:length(LWC_smoo)
  sonde_ILW = sonde_ILW + LWC_smoo(mm)*(secs_new(mm+1)-secs_new(mm))*abs(ascent_rate_new(mm)); 
end

disp(['sonde_ILW= ', num2str(sonde_ILW)]);

%---------------------------------
% plotting 
%---------------------------------
% first figure
figure
subplot(1,3,1);
plot(freq_new,alt_new);
axis tight;
hold on;
plot(freq_new_smoo,alt_new,'r-');
grid on;
ylim([0.0 5.5]);
hline = refline([0 0.803]);
hline.Color = [.7 .5 0];
xlabel('Freq [Hz]');
ylabel('Alt [km MSL]');
legend('Freq','Freq, 11-pt smoothed');

subplot(1,3,2);
plot(hzpermin_smoo,alt_new(1:end-1),'b-');
axis tight;
grid on;
xlabel('dFreq/dt [Hz/min]');
ylim([0.0 5.5]);
hline = refline([0 0.803]);
hline.Color = [.7 .5 0];

subplot(1,3,3);
plot(LWC_smoo,alt_new(1:end-1));
hold on;
%plot(lwc_smoo,alt_new(1:end-1));
axis tight;
xlim([0.0 0.5]);
%xlim([0.0 1.5]);
ylim([0.0 5.5]);
grid on;
xlabel('LWC [gm-3]');
hold on;
hline = refline([0 0.803]);
hline.Color = [.7 .5 0];
if str2num(date_str(1:8)) == 20140109
  if date_flight_num == 3
    text(0.40,1.20,'*')
  elseif date_flight_num == 4
    text(0.25,0.82,'*')
  end
end
%legend('LWC=1.1032*(dfreq/dt)^2+1.0275*dfreq/dt');
%title(title_str);

% second figure
figure

subplot(1,3,1);
%plot(lwc_smoo,alt_new(1:end-1));
plot(LWC_smoo,alt_new(1:end-1));
hold on;
axis tight;
xlim([0.0 0.5]);
%xlim([0.0 1.5]);
ylim([0.0 5.5]);
hline = refline([0 0.803]);
hline.Color = [.7 .5 0];
if radiom_exist == 1
  %RADIOMETER LWC OVERLAY
  %plot(radiom_LWC_data(1:end-1),radiom_alt_data+alt(1),'g-');
  %plot(radiom_LWC_data(1:end-1),radiom_alt_data+alt(1),'r-');
  plot(radiom_LWC_data(1:end),radiom_alt_data+alt(1),'r-');
  xlim([0.0 1.0]);
end
grid on;
xlabel('LWC [gm-3]');
hold on;
%legend('LWC SLW sonde','LWC radiom #3107');
%legend('LWC SLW sonde','LWC radiom #3058');
legend('LWC SLW sonde','LWC radiometer');

subplot(1,3,2);
plot(temp_cor_c(start:stop),alt(start:stop),'b-');
axis tight;
temp_max = max(temp_cor_c);
temp_min = min(temp_frost_c);
xlim([temp_min temp_max] );
ylim([0.0 5.5]);
hline = refline([0 0.803]);
hline.Color = [.7 .5 0];
%xlim([-15 0] );
hold on;
plot(temp_frost_c(start:stop),alt(start:stop),'r-');
grid on;
xlabel('Temp [deg C]');
ylabel('Alt [km MSL]');
hold on;
legend('sonde Ta','sonde Td');

subplot(1,3,3);
plot(humidity(start:stop),alt(start:stop));
axis tight;
grid on;
ylim([0.0 5.5]);
hline = refline([0 0.803]);
hline.Color = [.7 .5 0];
xlabel('Humidity [%]');
%title(title_str);

figure 
plot(ascent_rate(start:stop),alt(start:stop));
xlabel('Ascent Rate [ms^-^1]');
ylim([0.0 5.5]);
hline = refline([0 0.803]);
hline.Color = [.7 .5 0];
grid on;


