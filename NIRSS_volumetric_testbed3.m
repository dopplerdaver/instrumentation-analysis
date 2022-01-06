%-------------------------------------------------------------------
% Name:    NIRSS_volumetric_testbed3.m
%
% Purpose: Refactored developmental testbed code for converting to
%          VP to volumetric NIRSS product
%
% Inputs:  1. text file outputted from CTREC mdv file with
%             lat, lon, mrefl, uvec and vvec fields
%             -or-
%             text file from local ASOS with
%             wind speed and direction
%          2. NSSL 3D mosaic S-band reflectivity file 
%          3. daily NIRSS text log file
%
% Usage:   matlab -nodesktop
%          > NIRSS_volumetric_testbed2
%          > enter start date yyyymmddhhmmss
%          > enter stop date yyyymmddhhmmss
%
% Function calls:
%          plot_3085_radiom                - to read in scanning slant angle radiom
%          forecast_discussion_test_parser - reads text forecast discussion in order to
%                                            decide whether to interpolate hazard to 
%                                            all return not just advect packets
%
% Created: 1.11.2013 dserke 
%
% Coding tasks list:
% Task Status 
% 1.   0      update hardcopy flow chart of volumetric algo
% 2.   x      dynamic date input and processing
%      x        a. user inputted date/time range
%      x        b. auto processing and read of date/time input files 
% 3.   0      first guess haz for outer zones with S-band features
%      0        a. check for VP ILW with no S REFL 
%      0        b. meld first guess with advected packets
% 5.   0      change advection from KCLE CTREC vectors to current location vectors
% 6.   0      overlay S composite on CTREC vector array to QC feature tracking
% 7.   0      advect from 12 km slant haz as well as from VP?
% 8.   /      import other fields for NIRSSvol
%      0        a. haz quality? based on age and origin?
% 9.   /      multi-plot: red/yellow/green warning boxes for each approach/departure
% 10.  0      parallel plots for SLD warning 
% 11.  0      auto index retrieve for 3085 radiometer      
% 12.  x      add ASOS advection capability 
% 13.  /      integrate call and output of 'forecast_discussion_text_parser.m'
%               to testbed to determine when not to interpolate over whole 
%               radar domain
%
% Additional tasks:
%   why box lines giving errors?
%
%-------------------------------------------------------------------

[stat_mnt, res] = unix('mount /d2');

%-------------------------------------------------------------------
% hardcoded constants 
%-------------------------------------------------------------------

% CTREC
kcle_ctrec_lon    = 102;
kcle_ctrec_lat    = 165;

% NIRSS
nirss_radar_nhgts = 526;
nirss_rdiom_nhgts = 334;
nirss_ilw_vp_bias = 0.06;       % gm-3, from clear air comparisons to slant ILW
delta_kft         = 0.0984252;

% time conversions
sec_per_day       = 86400;

% metric conversions
kft_to_km         = 0.3048;
kts_to_mpersec    = 0.51444;

% radar S-band
kcle_lat_max      = 42.4;
kcle_lat_min      = 40.4;
kcle_lon_max      = -80.5;
kcle_lon_min      = -83.5;
min_betw_vols     = 5;
sec_betw_vols     = min_betw_vols * 60;

% hazard output
haz_grid_dim      = 50;
kcle_haz_x        = round(haz_grid_dim/2);
kcle_haz_y        = round(haz_grid_dim/2);
num_haz_zones     = 9;
hgt_min_zone      = [0 0 0 0 0 1700 2000 1500 2000];;                % all hgt mins and maxes in ft
hgt_max_zone      = [32874 3700 3000 3240 3500 4000 7000 5000 4000]; % all hgt mins and maxes in ft

% hazard plot specs
num_rows          = 5;
num_columns       = 5;
start_left        = 0.70;
start_bott        = 0.70;
size_plot         = 0.05;

%-------------------------------------------------------------------
% Prompt user for start and end yyyy, mm, dd
%------
% NOTE: TEST CASE TO BE USED 2/24/2012 22:29 MOD ICING PIREP OV CLE FL 6000 TA M3
%-------------------------------------------------------------------
disp('NOTE: Use hour long time blocks to not overrun multiple days, months or years');
disp('Enter start value        [example: 20120224215500]');
disp('Enter stop  value        [example: 20120224225500]');
input1       = num2str(input('  Enter start year       [yyyymmddhhmmss]:'));
input2       = num2str(input('  Enter  stop year       [yyyymmddhhmmss]:'));
disp('Enter advection method   [choices: 1 for "use only ASOS, use mosaic CTREC if no ASOS data"]');
disp('                         [         2 for "use only mosaic CTREC, use ASOS if no CTREC data"]');
disp('                         [         3 for "use ASOS if no echo in domain, when echo use CTREC"]');
input_advect = num2str(input('  Enter advection method [1,2 or 3]:'));

try
  yyyy1  = input1(1:4);
  mon1   = input1(5:6);
  day1   = input1(7:8);
  hour1  = input1(9:10);
  min1   = input1(11:12);
  sec1   = input1(13:14);
  yyyy2  = input2(1:4);
  mon2   = input2(5:6);
  day2   = input2(7:8);
  hour2  = input2(9:10);
  min2   = input2(11:12);
  sec2   = input2(13:14);
catch
  disp('  Length of input date is not correct ... ERROR!')
end

% determine how many 5 minute time blocks user has bracketed
date1           = datenum(str2num(yyyy1),str2num(mon1),str2num(day1),str2num(hour1),str2num(min1),str2num(sec1));
date2           = datenum(str2num(yyyy2),str2num(mon2),str2num(day2),str2num(hour2),str2num(min2),str2num(sec2));

time_diff       = (date2 - date1) * sec_per_day;
five_min_blocks = ceil(time_diff/300);
hour            = hour1;
minute          = min1;

%-------------------------------------------------------------------
% Define relevant dirs and filenames
%-------------------------------------------------------------------
case_date        = [yyyy1 mon1 day1];
ctrec_data_path  = '/d2/ctrec/data/';
ctrec_date       = [case_date '/'];
ctrec_file       = [];
mos3d_file       = [];

for zz = 1:five_min_blocks+1
  mos3d_file = [mos3d_file; yyyy1 mon1 day1 '_' num2str(hour) minute '00_tile3mosaic.nc'];
  desired_hhmm{zz} = [num2str(hour) minute];
  minute           = str2num(minute) + 5;
  if minute >= 60
    minute = '00';
    try
      hour   = str2num(hour) + 1;
    catch
      hour   = hour + 1;
    end
  else
    minute = num2str(minute);
    if length(minute) == 1 
      minute = ['0' minute];
    end
  end
  ctrec_file = [ctrec_file; num2str(hour) num2str(minute)]; 
end

mos3d_file = cellstr(mos3d_file);
ctrec_file = cellstr(ctrec_file);

mos3d_data_path  = '/d2/3dmosaic/';
nirss_data_path  = '/d1/serke/projects/NIRSS_NASA/data/volumetric_prod_case_data/';
nirss_date       = [case_date '_'];

%-------------------------------------------------------------------
% Test number of times for each input type 
%-------------------------------------------------------------------
inputs = [length(ctrec_file) length(mos3d_file)];

if inputs(1) ~= inputs(2) 
  disp('  *WARNING: Number of input files not equal');
  ind_inputs = find(inputs == min(inputs));
  if ind_inputs == 1
    input_type = 'ctrec_file';
  elseif ind_inputs == 2
    input_type = 'mos3d_file';
  end
  disp(['  Minimum files are from input source: ' input_type(1:5)]);
  num_times = length(eval(input_type));
else
  disp('  Equal number of CTREC and 3Dmosaic inputs available');
  num_times = length(ctrec_file);
end

%-------------------------------------------------------------------
% Check if daily 3085 exists locally, if not then get it
%-------------------------------------------------------------------
try
  disp('  Check if desired 3085 radiometer file exists locally ...');
  radiom3085_yyyymmdd_file_exists = 0;
  cd /d1/serke/projects/NIRSS_NASA/data/volumetric_prod_case_data/;
  string_3085 = ['!ls *.3085 > 3085_files.txt'];
  eval(string_3085);
  fid = fopen('3085_files.txt');
  existing_date_3085files = textscan(fid,'%d');
  fclose(fid);
  for ff = 1:length(existing_date_3085files{1})
    existing_3085file  = num2str(existing_date_3085files{1}(ff));
    existing_3085_mmdd = existing_3085file(5:8);
    existing_3085_yyyy = existing_3085file(1:4);
    if yyyy1 == existing_3085_yyyy & [mon1 day1] == existing_3085_mmdd 
      radiom3085_yyyymmdd_file_exists = 1;
    end
  end

  if radiom3085_yyyymmdd_file_exists == 1
    disp('    Desired 3085 text file does exist locally...');
  else
    disp('    Desired 3085 text file does NOT exist locally...');

    try
      disp('    Retrieve desired 3085 file ...');
      string8 = ['!scp nirss@jlab-nima1.rap.ucar.edu:/d2/nirss/nasa_glen/archive/' yyyy1 mon1 day1 '.3085 /d1/serke/projects/NIRSS_NASA/data/volumetric_prod_case_data' ];
      eval(string8);
      radiom3085_yyyymmdd_file_exists = 1;

    catch
      disp('    Retrieve desired 3085 file ... FAILED');
      disp('    ');
    end

  end
catch
end

%keyboard;

%-------------------------------------------------------------------
% Read in daily 3085 file 
%-------------------------------------------------------------------
%try
  disp('  Read in desired 3085 radiometer daily file...');
  cd /d1/serke/projects/NIRSS_NASA/code/mat;
  plot_3085radiom;
  disp('    Read in desired 3085 radiometer daily file...COMPLETE');
%catch
%  disp('    Read in desired 3085 radiometer daily file...FAILED');
%end

%keyboard;

%-------------------------------------------------------------------
% Read in daily nirss output file 
%-------------------------------------------------------------------
try
  disp('  Check if desired NIRSS daily file exists locally ...');
  nirss_yyyymmdd_file_exists = 0;
  cd /d1/serke/projects/NIRSS_NASA/data/volumetric_prod_case_data/;
  string09 = ['!ls ' yyyy1 mon1 day1 '.nirss > ' yyyy1 mon1 day1 '_nirss_files.txt'];
  eval(string09);
  existing_date_nirssfiles       = importdata([yyyy1 mon1 day1 '_nirss_files.txt']);
  for dd = 1:length(existing_date_nirssfiles)
    existing_nirss_yyyy{dd} = existing_date_nirssfiles{dd}(1:4);
    existing_nirss_mmdd{dd} = existing_date_nirssfiles{dd}(5:8);
    if yyyy1 == existing_nirss_yyyy{dd} & [mon1 day1] == existing_nirss_mmdd{dd} 
      nirss_yyyymmdd_file_exists = 1;
      disp('    Desired NIRSS text file does exist locally...');
    end
  end
catch
  disp('  Check if desired NIRSS daily file exists locally ...FAILED');
  dbstop if error;
end

if nirss_yyyymmdd_file_exists == 0

  try
    disp('  Retrieve desired NIRSS daily file...');
    string10 = ['!scp nirss@jlab-nima1.rap.ucar.edu:/d2/nirss/nasa_glen/nirss_output/' yyyy1 mon1 day1 '.txt /d2/nirss_output' ];
    eval(string10);
    disp('  Retrieve desired NIRSS daily file...COMPLETE');
  catch
    disp('  Retrieve desired NIRSS daily file...FAILED');
    dbstop if error;
  end   % end of try/catch

end     % end of if nirss_yyyymmdd_file_exists == 0

try
  disp('  Read in desired NIRSS daily file...');
  cd /d2/nirss_output;
  nirss_data = importdata([nirss_date(1:8) '.txt']);
  cd /d1/serke/projects/NIRSS_NASA/code/mat;
  disp('    Read in desired NIRSS daily file...COMPLETE');
catch
  disp('    Read in desired NIRSS daily file...FAILED');
  dbstop if error;
end

%-------------------------------------------------------------------
% Check if daily ASOS file exists locally, if not then get it
%-------------------------------------------------------------------
disp('  Check if desired ASOS file exists locally ...');
cd /d1/serke/projects/NIRSS_NASA/data/ASOS
string_asos = ['!ls *ASOS_2* > ASOS_files.txt'];
eval(string_asos);
fid = fopen('ASOS_files.txt');
existing_date_asosfiles = textscan(fid,'%s');
fclose(fid);
asos_yyyymmdd_file_exists = 0;
asos_yyyymmdd_file_new    = 0;
for ff = 1:length(existing_date_asosfiles{1})
  existing_asosfile  = char(existing_date_asosfiles{1}(ff));
  existing_asos_mmdd = existing_asosfile(20:23);
  existing_asos_yyyy = existing_asosfile(16:19);
  if yyyy1 == existing_asos_yyyy & [mon1 day1] == existing_asos_mmdd 
    asos_yyyymmdd_file_exists = 1;
  end
end

if asos_yyyymmdd_file_exists == 1
  disp('    Desired ASOS text file does exist locally...');
else
  disp('    Desired ASOS text file does NOT exist locally...');

  try
    disp('    Retrieve desired ASOS file ...');
    cd /d1/serke/projects/NIRSS_NASA/code;
    system(['./winter_weather_data_query "' yyyy1 '-' mon1 '-' day1 ' 00:00:00" "' yyyy1 '-' mon1 '-' day1 ' 23:59:59" asos KCLE > /d1/serke/projects/NIRSS_NASA/data/ASOS/NIRSS_vol_ASOS_' yyyy1 mon1 day1 ]);
    cd /d1/serke/projects/NIRSS_NASA/code/mat;
    disp('    Retrieve desired ASOS file ... COMPLETE');
    asos_yyyymmdd_file_exists = 1;
    asos_yyyymmdd_file_new    = 1;

  catch
    disp('    Retrieve desired ASOS file ... FAILED');
    disp('    ');
  end

end

%-------------------------------------------------------------------
% Read in daily ASOS file 
%-------------------------------------------------------------------
try
  disp('  Read in desired ASOS daily file...');
  cd /d1/serke/projects/NIRSS_NASA/data/ASOS;
  if asos_yyyymmdd_file_new == 1
    [stat_rmhdr, res]  = unix(['gvim -s /d2/ctrec/data/rm_header.vim NIRSS_vol_ASOS_' yyyy1 mon1 day1]);
  end
  fid_asos = fopen(['NIRSS_vol_ASOS_' yyyy1 mon1 day1]);
  asos     = textscan(fid_asos, '%s%f%s%f%s%f%s%d%d%d%d%s%s%d%f%d%f%f%f%d%d','delimiter',',','treatasempty','NULL');
  fclose(fid_asos);
  disp('    Read in desired ASOS daily file...COMPLETE');
catch
  disp('    Read in desired ASOS daily file...FAILED');
end

cd /d1/serke/projects/NIRSS_NASA/code/mat;

%-------------------------------------------------------------------
% loop over number of 5 minute 3D mosaic and CTREC times 
%-------------------------------------------------------------------
for ii = 1:num_times

  %----------------------------------------------------------
  % extract desired 3085 data from desired time
  %----------------------------------------------------------
  try
    disp(['  Search for 3085 fields at desired time...' desired_hhmm{ii}]);
    found           = 0;
    time_index_3085 = 0;
    for yy=1:length(date_time)

        % test if 2-digit hour is the same
        if date_time{yy}(10:11) == desired_hhmm{ii}(1:2) & found == 0 

          % test if 2-digit minute is the same or very close
          if str2num(date_time{yy}(13:14)) == str2num(desired_hhmm{ii}(3:4)) & found == 0
            ilw_slant{ii}   = [ilw_w(yy) ilw_n(yy) ilw_e(yy) ilw_s(yy)];
            found           = 1;
            time_index_3085 = yy; 
            disp('    Search for 3085 fields at desired time...exact match COMPLETE');
          elseif str2num(date_time{yy}(13:14)) ~= str2num(desired_hhmm{ii}(3:4)) & found == 0 & ...
                 abs(str2num(date_time{yy}(13:14))-str2num(desired_hhmm{ii}(3:4))) <=3
            ilw_slant{ii}   = [ilw_w(yy) ilw_n(yy) ilw_e(yy) ilw_s(yy)];
            found           = 1;
            time_index_3085 = yy; 
            disp('    Search for 3085 fields at desired time...proximate match COMPLETE');
          end

        end   % end of date_time{yy}(10:11) == ... 
 
    end       % end of for yy=1:length(date_time)

  catch
    disp(  'Search for 3085 fields at desired time...FAILED');
  end

  %----------------------------------------------------------
  % Read in iith NIRSS hazard fields from text files 
  %----------------------------------------------------------
  for kk = 1: nirss_radar_nhgts;
    nirss_radar_hgts_kft(kk) = delta_kft * kk;
  end
  for ll = 1: nirss_rdiom_nhgts;
    nirss_rdiom_hgts_kft(ll) = delta_kft * ll;
  end

  nirss_radar_hgts_km = nirss_radar_hgts_kft * kft_to_km;
  nirss_rdiom_hgts_km = nirss_radar_hgts_kft * kft_to_km;

  %----------------------------------------------------------
  % Check if iith 3D mosaic nc and mdv files exists 
  %----------------------------------------------------------
  disp('  Check if desired mosaic nc file exists...');

  cd /d2/3dmosaic;

  nc_yyyymmddhhmm_file_exists  = 0;
  mdv_yyyymmddhhmm_file_exists = 0;

  % get existing ncfiles list
  string20 = ['!ls ' yyyy1 mon1 day1 '* > ' yyyy1 mon1 day1 'files.txt'];
  eval(string20);
  existing_date_ncfiles       = importdata([yyyy1 mon1 day1 'files.txt']);

  % loop over all existing ncfiles
  for ee = 1:length(existing_date_ncfiles)

    existing_nc_yyyy{ee} = existing_date_ncfiles{ee}(1:4);
    existing_nc_mmdd{ee} = existing_date_ncfiles{ee}(5:8);
    existing_nc_hhmm{ee} = existing_date_ncfiles{ee}(10:13);
   
    % check if desired mosaic ncfile exists locally
    if yyyy1 == existing_nc_yyyy{ee} & [mon1 day1] == existing_nc_mmdd{ee} & desired_hhmm{ii} == existing_nc_hhmm{ee} 
      nc_yyyymmddhhmm_file_exists = 1;
    end

  end    % end of for ee=1:length(existing_date_ncfiles_

  if nc_yyyymmddhhmm_file_exists == 1
    disp('    Desired mosaic nc file does exist...');
  else
    disp('    Desired mosaic nc file does NOT exist...');

    % check if mdv dir for desired date exists locally
    try
      string40 = ['cd mdv/nsslMosaic3D/tile3/' num2str(yyyy1) num2str(mon1) num2str(day1) ];
      eval(string40);
      disp('    Desired mosaic mdv directory does exist');
    catch
      disp('    Desired mosaic mdv directory does NOT exist');

      % Retrieve desired mosaic file from hsi and manipulate it 
      try 
        disp('  Getting mosaic file from hsi...');
        !kinit;
        string2 = ['!hsi "cd /RAPDMG/CWX/nsslMosaics/' num2str(yyyy1) '/' num2str(mon1) num2str(day1) '; get ' num2str(yyyy1) num2str(mon1) num2str(day1) '_tile3_nsslMosaic3D.tar"']; 
        eval(string2);
        disp('  Getting mosaic file from hsi...COMPLETE');
      catch
        disp('  Getting mosaic file from hsi...FAILED');
      end
  
      try
        disp('  Untarring tar file...');
        [stat,res] = unix(['tar -xvf ' num2str(yyyy1) num2str(mon1) num2str(day1) '_tile3_nsslMosaic3D.tar']);
        disp('  Untarring tar file...COMPLETE');
      catch
        disp('  Untarring tar file...FAILED');
      end
    
      try
        disp('  Remove tar file...');
        string3 = ['!rm ' num2str(yyyy1) num2str(mon1) num2str(day1) '_tile3_nsslMosaic3D.tar'];
        eval(string3);
        disp('  Remove tar file...COMPLETE');
      catch
        disp('  Remove tar file...FAILED');
      end
  
      try
        disp('  Unzipping files...');
        string4 = ['cd mdv/nsslMosaic3D/tile3/' num2str(yyyy1) num2str(mon1) num2str(day1) ];
        eval(string4);
        hour_diff = str2num(hour2) - str2num(hour1);
        for ww = 0:hour_diff
          string5 = ['gunzip ' num2str(str2num(hour1)+ww) '*'];
          eval(string5);
        end 
        disp('  Unzipping files...COMPLETE');
      catch
        disp('  Unzipping files...FAILED');
      end
    end

    % make a text file with all unzipped mdvfiles listed
    string40 = ['!ls *.mdv > uzip_mdvfiles.txt'];
    eval(string40);

    existing_date_uzip_mdvfiles = importdata(['uzip_mdvfiles.txt']);

    % loop over all unzipped mdvfiles 
    for ff = 1:length(existing_date_uzip_mdvfiles)
      existing_mdv_uzip_hhmmcell{ff} = existing_date_uzip_mdvfiles{ff}(1:4);
    end

    existing_uzip_mdv_hhmm = cell2mat(existing_mdv_uzip_hhmmcell);
    existing_uzip_mdv_hhmm = str2num(strread(existing_uzip_mdv_hhmm,'%4c'));

    % check if desired uzipped mdv file exists
    ind_desired_uzip_mdv = find(existing_uzip_mdv_hhmm == str2num(desired_hhmm{ii}));
    if length(ind_desired_uzip_mdv) >= 1
      disp('    Desired mosaic unzipped mdv file does exist...');
      mdv_hhmm = existing_uzip_mdv_hhmm(ind_desired_uzip_mdv); 
    else
      % find closest unzipped mdv time if exact match does not exist
      disp('    Desired mosaic unzipped file does NOT exist...');
      disp('       Finding and using next earliest mdv file...');
      ind_previousuzipmdv = find(existing_uzip_mdv_hhmm < str2num(desired_hhmm{ii}));
      ind_previousuzipmdv = ind_previousuzipmdv(end);
      mdv_hhmm            = existing_uzip_mdv_hhmm(ind_previousuzipmdv); 
    end  % end of length(ind_desiredmdv) >= 1

    disp('  Convert desired mosaic mdv to netcdf ...');
    try
      cd /d2/3dmosaic;
      %cd ../../../..;
      minutes = num2str(str2num(min1)+(ii*5)-5);
      if length(minutes) == 1; minutes = ['0' minutes]; end;
      string6 = ['!mdv2netCDF -params /d2/3dmosaic/mdv2netCDF.params -if /d2/3dmosaic/mdv/nsslMosaic3D/tile3/' yyyy1 mon1 day1 '/' num2str(mdv_hhmm) '00.mdv'];
      eval(string6);
      disp('  Convert desired mosaic mdv to netcdf...COMPLETE');
    catch
      disp('  Convert desired mosaic mdv to netcdf...FAILED');
    end
          
  end    % end of if yyyy1 = ... 
  
  %----------------------------------------------------------
  % Read in iith 3D mosaic file 
  %----------------------------------------------------------
  try
    disp('  Read in desired mosaic nc file...');
    cd /d1/serke/projects/NIRSS_NASA/code/mat;
    nc = NetcdfRead([mos3d_data_path mos3d_file{ii}]);
    disp('    Read in desired mosaic nc file...COMPLETE');
  catch
    disp('    Read in desired mosaic nc file...FAILED');
    dbstop if error;
  end

  %----------------------------------------------------------
  % strip out refl profiles at several azimuths at 10 km and 50 km ranges
  %----------------------------------------------------------
  mosaic_hgt               = nc.vars.height.data;
  % profiles at 10 km range
  mosaic_refl_prof_NW_10km = nc.vars.mrefl_mosaic.data(:,:,150,807); 
  mosaic_refl_prof_NE_10km = nc.vars.mrefl_mosaic.data(:,:,150,824); 
  mosaic_refl_prof_SE_10km = nc.vars.mrefl_mosaic.data(:,:,137,824); 
  mosaic_refl_prof_SW_10km = nc.vars.mrefl_mosaic.data(:,:,137,807); 
  % profiles at 10 km range
  
  %----------------------------------------------------------
  % compare 10 km proximate mosaic profiles to NIRSS Ka profile 
  %----------------------------------------------------------

  %----------------------------------------------------------
  % some pre-testing before choosing either ASOS or 3D mosaic CTREC advection method 
  %----------------------------------------------------------
  if str2num(input_advect) == 3 & size(asos,1) == 1 

    disp('  Determining whether there is significant trackable echo in domain...');
    echo_in_domain = 0;
    refl_domain_within50km_KCLE = nc.vars.mrefl_mosaic.data(:,3,94:194,765:865);
    ind_refl                    = find(refl_domain_within50km_KCLE > -30);
    if size(ind_refl,1) > 0.1 * size(refl_domain_within50km_KCLE,3)*size(refl_domain_within50km_KCLE,4) 
      echo_in_domain = 1;
      disp('    Determining whether there is significant trackable echo in domain...AFFIRMATIVE');
    else
      disp('    Determining whether there is significant trackable echo in domain...NEGATIVE');
    end
    clear nc;

  elseif (str2num(input_advect) == 1 | str2num(input_advect) == 3) & size(asos,1) == 0 

    disp('WARNING: User chose ASOS advection method but ASOS file does not exist');
    disp('         Attempting to use only CTREC 2D Mosaic advection...');

  end

  %----------------------------------------------------------
  % choose either ASOS or 3D mosaic CTREC advection method 
  %----------------------------------------------------------
  if (str2num(input_advect) == 1 & size(asos,1) == 1) | (str2num(input_advect) == 3 & echo_in_domain == 0)

    try
      disp('  Attempting to process time with ASOS advection...');

      % find index of asos times that matches desired time of current 5-min block
      for lll=1:length(asos{1})
        date_time_asos    = char(asos{1}(lll));
        date_time_asos_hr = date_time_asos(12:13);
        date_time_asos_mi = date_time_asos(15:16);
        if [date_time_asos_hr date_time_asos_mi] == char(desired_hhmm{ii})
          asos_ind = lll; 
        end
      end    % end of for lll=1:length(asos{1})
    
      % convert ASOS native units knots to ms-1
      asos_mpersec = asos{9}(lll) * kts_to_mpersec;
      % convert ASOS native polar coords to cartesian coords
      [vvec,uvec]  = pol2cart(double(asos{8}(lll)) * pi/180, double(asos_mpersec));
      disp('Attempting to process time with ASOS advection...COMPLETE');
    catch
      disp('Attempting to process time with ASOS advection...FAILED');
    end

  elseif (str2num(input_advect) == 2) | ((str2num(input_advect) == 1 | str2num(input_advect) == 3) & size(asos,1) == 0) | (str2num(input_advect) == 3 & size(asos,1) == 1 & echo_in_domain == 1) 

    disp('  Attempting to process time with CTREC 3D Mosaic advection...');
   
    %----------------------------------------------------------
    % Check if iith CTREC mdv and txt files exist 
    %----------------------------------------------------------
    disp('  Check if desired CTREC directory exists...');
    string68 = ['cd ' ctrec_data_path];
    eval(string68);
    dir_to_find = [yyyy1 mon1 day1];
    if isequal(exist(dir_to_find, 'dir'),7)
      disp('     Desired CTREC directory does exist')
    else
      disp('     Desired CTREC directory does NOT exist')
      disp('     Making desired CTREC directory')
      string69=['!mkdir ' yyyy1 mon1 day1];
      eval(string69);
    end
    
    try
      disp('  Check if desired CTREC file exists...');
  
      string70 = ['cd ' ctrec_data_path ctrec_date ];
      string71 = ['!ls *.mdv > ctrec_mdv_files.txt'];
      string72 = ['!ls *.txt > ctrec_nohdrtxt_files.txt'];
      eval(string70);
      eval(string71);
      eval(string72);
  
      ctrecmdv_yyyymmddhhmm_file_exists = 0;
      existing_date_ctrecmdvfiles       = importdata(['ctrec_mdv_files.txt']);
      for gg = 1:length(existing_date_ctrecmdvfiles)
        existing_ctrecmdv_hhmm{gg} = existing_date_ctrecmdvfiles{gg}(1:4);
        if existing_ctrecmdv_hhmm{gg} == ctrec_file{ii}(1:4) 
          ctrecmdv_yyyymmddhhmm_file_exists = 1;
          disp('    Desired CTREC mdv file does exist...');
        end
      end
    catch
      disp('    Desired CTREC mdv file does NOT exist...');
    end
  
    ctrectxt_yyyymmddhhmm_file_exists = 0;
    existing_date_ctrectxtfiles       = importdata(['ctrec_nohdrtxt_files.txt']);
    for gg = 1:length(existing_date_ctrectxtfiles)
      existing_ctrectxt_hhmm = existing_date_ctrectxtfiles{gg}(1:4);
      %if ii+1 <= length(existing_date_ctrectxtfiles)
        if existing_ctrectxt_hhmm == ctrec_file{ii}(1:4) 
          ctrectxt_yyyymmddhhmm_file_exists = 1;
          disp('    Desired CTREC text file does exist...');
        end
      %end
    end
  
    %-------------------------------------------------------------------
    % Process iith CTREC file if it does not already exist 
    %-------------------------------------------------------------------
    cd ~/cvs/apps/trec/src/ctrec;
  
    if ctrecmdv_yyyymmddhhmm_file_exists == 0 & ctrectxt_yyyymmddhhmm_file_exists ==0 
  
      disp('    CTREC mdv file does not exist.  Attempting to process...');
  
      try
        % run ctrec starting 10 minutes before desired hhmm
  
        if str2num(ctrec_file{ii}(3:4)) < 10
          start_hour = num2str(str2num(ctrec_file{ii}(1:2))-1);
          start_min  = num2str(str2num(ctrec_file{ii}(3:4)) + 60 - 10);
        elseif num2str(ctrec_file{ii}(3:4)) == 10  
          start_hour = ctrec_file{ii}(1:2);
          start_min  = '00';
        else
          start_hour = ctrec_file{ii}(1:2);
          start_min  = num2str(str2num(ctrec_file{ii}(3:4)) - 10);
        end
  
        [stat_ctrec, res]  = unix(['./ctrec -params ctrec_params_serke -mode ARCHIVE -debug verbose -start "' yyyy1 ' ' mon1 ' ' day1 ' ' start_hour ' ' start_min ' ' sec1 '" -end "' yyyy2 ' ' mon2 ' ' day2 ' ' ctrec_file{ii}(1:2) ' ' ctrec_file{ii}(3:4) ' ' sec2 '"']); 
        disp('  Processing of desired CTREC file...COMPLETE');
      catch
        disp('  Processing of desired CTREC file...FAILED');
        disp('    Try processing file manually!');
        dbstop if error;
      end   
  
      try
        disp('  Converting CTREC mdv file to text...');
        eval(string70);
        string80           = ['!PrintMdv -f ' ctrec_file{ii}(1:4) '00.mdv -table > ' ctrec_file{ii}(1:4) '.txt'];
        eval(string80);
        [stat_rmhdr, res]  = unix(['gvim -s /d2/ctrec/data/rm_header.vim ' ctrec_file{ii} '.txt']);
        disp('  Converting CTREC mdv file to text...COMPLETE');
      catch
        disp('  Converting CTREC mdv file to text...FAILED');
        dbstop if error;
      end
  
    elseif ctrecmdv_yyyymmddhhmm_file_exists == 1 & ctrectxt_yyyymmddhhmm_file_exists ==0
  
      try
        disp('  Converting CTREC mdv file to text...');
        string80 = ['!PrintMdv -f ' ctrec_file{ii}(1:4) '00.mdv -table > ' ctrec_file{ii}(1:4) '.txt'];
        eval(string80);
        [stat_rmhdr, res]  = unix(['gvim -s /d2/ctrec/data/rm_header.vim ' ctrec_file{ii} '.txt']);
        disp('  Converting CTREC mdv file to text...COMPLETE');
      catch
        disp('  Converting CTREC mdv file to text...FAILED');
        dbstop if error;
      end
  
    elseif ctrecmdv_yyyymmddhhmm_file_exists == 1 & ctrectxt_yyyymmddhhmm_file_exists ==1
  
      eval(string70);
      [stat_rmhdr, res]  = unix(['gvim -s /d2/ctrec/data/rm_header.vim ' ctrec_file{ii} '.txt']);
  
    end     % end of if ctrec_yyyymmddhhmm_file_exists == 0
  
    %-------------------------------------------------------------------
    % Attempt to read in iith CTREC output text file 
    %-------------------------------------------------------------------
    clear ctrec_data;
    eval(string70);
    pause(20);
    %keyboard;
    try
      disp('  Read in desired CTREC file...');
      ctrec_data = load([ctrec_data_path ctrec_date ctrec_file{ii} '.txt']);
      disp('    Read in desired CTREC file...COMPLETE');
    catch
      disp('    Read in desired CTREC file...FAILED');
      dbstop if error;
    end      
  
    cd /d1/serke/projects/NIRSS_NASA/code/mat;
  
    %----------------------------------------------------------
    % strip out CTREC fields 
    %----------------------------------------------------------
    hgt        = ctrec_data(:,1);
    lat        = ctrec_data(:,2);
    lon        = ctrec_data(:,3);
    refl       = ctrec_data(:,4);
    uvec       = ctrec_data(:,5);
    vvec       = ctrec_data(:,6);
    clear ctrec_data;

    %-------------------------------------------------------------------
    % Crop CTREC lat/lon ranges to KCLE area of interest 
    %-------------------------------------------------------------------
    count  = 1;
    count2 = 1;
    
    for j = 1:size(uvec,1)
      if lat(j) >= kcle_lat_min & lat(j) <= kcle_lat_max & ...
         lon(j) <= kcle_lon_max & lon(j) >= kcle_lon_min
        hgt_cle(count)  = hgt(j);
        lat_cle(count)  = lat(j);
        lon_cle(count)  = lon(j);
        refl_cle(count) = refl(j);
        uvec_cle(count) = uvec(j);
        vvec_cle(count) = vvec(j);
        count           = count + 1;
      end
    end
    clear ctrec_data hgt lat lon refl uvec vvec;
  
    % reorganize arrays
    for k = 1:201
      for l = 1:301
        lat(k,l)  = lat_cle(count2);
        lon(k,l)  = lon_cle(count2);
        refl(k,l) = refl_cle(count2);
        uvec(k,l) = uvec_cle(count2);
        vvec(k,l) = vvec_cle(count2);
        count2    = count2 + 1;
      end
    end
  
  end      % end of if input_advect == 1

  %----------------------------------------------------------
  % Based on Zhang's AMS conf paper (Fig 5), a precip VCP at 10km 
  %    horiz from radar will have a 19.5 degree tilt up to 3.5 km 
  %    alt AGL and a clear air VCP at 10km horiz 4.5 degree tilt up to 1km alt AGL.
  % Based on Zhang's 1997 JAOT paper (Fig 4), a precip VCP at 50 km 
  %    horiz from radar will have a 19.5 degree tilt up to 18 km 
  %    alt AGL and a clear air VCP at 50 km horiz (4.5 degree tilt) up to 5 km alt AGL. 
  %
  % GOAL: divide plan view into 90 degree azimuth quadrants within ranges:
  %    radial dist:
  %    0   to 2.5 km - covered by VPR NIRSS
  %    2.5 to 15 km  - covered by interpolating VPR NIRSS to 10 km mosaic profile and 3085 ILW
  %    15 to 50 km   - covered by interpolating VPR NIRSS to 50 km mosaic profiles
  %----------------------------------------------------------

  %----------------------------------------------------------
  % set up first guess arrays for NIRSS volumetric hazards 
  %
  % Hazard scale:                  0-8, -1 is unknown
  %
  % NIRSS volumetric hazard array: [haz_ZE haz_NW_10km haz_NE_10km haz_SE_10km ...
  %                                 haz_SW_10km haz_NW_50km haz_NE_50km ...
  %                                 haz_SE_50km haz_SW_50km];
  %----------------------------------------------------------
  
  %----------------------------------------------------------
  % extract desired NIRSS data from desired time
  %----------------------------------------------------------
  drops_nirss = [];

  try
    disp(  '  Search for NIRSS fields at desired time');
    for yy=1:length(nirss_data)
      if nirss_data{yy}(1:5) == 'RADAR' 
        if nirss_data{yy}(7:8) == desired_hhmm{ii}(1:2) 
          if nirss_data{yy}(10:11) == desired_hhmm{ii}(3:4)
            line1          = textscan(nirss_data{yy}(43:end),'%7.4f',nirss_radar_nhgts*2);
            radar_nirss    = line1{1}(1:nirss_radar_nhgts);
            cld_mask_nirss = line1{1}(nirss_radar_nhgts+1:end);
            nirss_ilw      = str2num(nirss_data{yy+1}(32:35));
            line3          = textscan(nirss_data{yy+3}(41:end),'%10.9f',nirss_radar_nhgts*2);
            lwc_nirss      = line3{1}(1:nirss_radar_nhgts);
            haz_nirss{ii}  = line3{1}(nirss_radar_nhgts+1:end);

            if nirss_data{yy+4}(1:5) == 'DROPS'
              line4       = textscan(nirss_data{yy+4}(41:end),'%7.4f',nirss_radar_nhgts*2);
              drops_nirss = line4{1}(1:nirss_radar_nhgts);
              disp(  '    Search for NIRSS drops field at desired time...COMPLETE');
            else
              disp(  '    Search for NIRSS drops field at desired time...FAILED');
            end

            disp(  '    Search for NIRSS fields at desired time...COMPLETE');

          end
        end 
      end

    end

  catch
    disp(  '    Search for NIRSS fields at desired time...FAILED');
  end

  % remove high bias in VP ILW
  nirss_ilw(ii)                     = nirss_ilw - nirss_ilw_vp_bias;

  NIRSS_volumetric_lwc_array{ii}    = ones(1,num_haz_zones) * -1;
  NIRSS_volumetric_lwc_array{ii}(1) = max(lwc_nirss);
  NIRSS_volumetric_hazcol_array{ii} = ones(1,num_haz_zones) * -1;
  ind1                              = find(cld_mask_nirss == 1);
  fiveper                           = ceil(length(ind1) * 0.05);
  haz_nirss_max                     = max(haz_nirss{ii});
  ind2                              = find(haz_nirss{ii} == haz_nirss_max);
  clear ind1;
  
  % check if max of NIRSS hazard is an outlier in profile
  %   and set max haz in col for zone 1 accordingly
  if length(ind2) >= fiveper
    NIRSS_volumetric_hazcol_array{ii}(1) = max(haz_nirss{ii}); 
  else
    NIRSS_volumetric_hazcol_array{ii}(1) = max(haz_nirss{ii}) - 1;
  end
  
  NIRSS_volumetric_hazbox_array{ii} = ones(1,num_haz_zones) * -1;
  haz_nirss_plus1                   = haz_nirss{ii} + 1;
  ind_plus                          = find(haz_nirss_plus1 == 1);
  haz_nirss_plus1(ind_plus)         = 0;
  haz_nirss_minu1                   = haz_nirss{ii} - 1;
  ind_minu                          = find(haz_nirss_minu1 == -1);
  haz_nirss_minu1(ind_minu)         = 0;

  % this field will help determine if should search through 
  %   normal, increased or decreased nirss hazard profile
  compare_zone_ilw                  = 0;

  for bb = 1:num_haz_zones

    % do this chunk only for zones 2-5
    if bb >= 2 & bb <= 5

      % if slant within +- 20%, set max hazard value equal to VP
      if ilw_slant{ii}(bb-1) < nirss_ilw(ii)*1.2 & ilw_slant{ii}(bb-1) > nirss_ilw(ii)*0.8
        NIRSS_volumetric_hazcol_array{ii}(bb) = NIRSS_volumetric_hazcol_array{ii}(1);
        compare_zone_ilw                      = 0;
      % if slant much different than VP, adjust slant hazard up or down
      elseif ilw_slant{ii}(bb-1) > nirss_ilw(ii)*1.2
        NIRSS_volumetric_hazcol_array{ii}(bb) = NIRSS_volumetric_hazcol_array{ii}(1) + 1;
        compare_zone_ilw                      = 1;
      else 
        NIRSS_volumetric_hazcol_array{ii}(bb) = NIRSS_volumetric_hazcol_array{ii}(1) - 1;
        compare_zone_ilw                      = -1;
      end
      % if slant near noise level, set slant hazard to 0
      if ilw_slant{ii}(bb-1) < 0.06
        NIRSS_volumetric_hazcol_array{ii}(bb) = 0; 
        compare_zone_ilw                      = 0;
      end;

    end   % end of if bb>=2&bb<=5

    % find indices of heights within bbth zone box
    ind1 = find(nirss_rdiom_hgts_kft >= hgt_min_zone(bb)/1000 & ...
                nirss_rdiom_hgts_kft <= hgt_max_zone(bb)/1000);

    % set max nirss haz in each zone box (horiz and alt) to normal, increased or decreased haz prof
    if NIRSS_volumetric_hazcol_array{ii}(bb) >= 0
      if compare_zone_ilw == 0
        NIRSS_volumetric_hazbox_array{ii}(bb) = max(haz_nirss{ii}(ind1));
      elseif compare_zone_ilw == 1
        NIRSS_volumetric_hazbox_array{ii}(bb) = max(haz_nirss_plus1(ind1));
      else
        NIRSS_volumetric_hazbox_array{ii}(bb) = max(haz_nirss_minu1(ind1));
      end
    else
      NIRSS_volumetric_hazbox_array{ii}(bb) = NIRSS_volumetric_hazcol_array{ii}(bb);
    end
    clear ind1;

  end     % end of for bb=1:num_haz_zones

  %----------------------------------------------------------
  % Create hazard data packet structure and fields 
  %
  % key: 
  %     source: n=nirss,        k=ka-band,   s=s-band, c=ctrec
  %     type:   max=max in col, prf=profile, flg=flag, val=value
  %----------------------------------------------------------
  NIRSSvol(ii).date_time       = date_time(time_index_3085);
  NIRSSvol(ii).LIWCON_n_max    = NIRSS_volumetric_lwc_array{ii};
  NIRSSvol(ii).LIWCON_n_prf    = lwc_nirss;
  NIRSSvol(ii).HAZARD_n_maxcol = NIRSS_volumetric_hazcol_array{ii};
  NIRSSvol(ii).HAZARD_n_maxbox = NIRSS_volumetric_hazbox_array{ii};
  NIRSSvol(ii).HAZARD_n_prf    = haz_nirss{ii};
  NIRSSvol(ii).DROPSZ_n_prf    = drops_nirss;
  NIRSSvol(ii).PRECIP_n_flg    = [];
  NIRSSvol(ii).BIMODE_k_prf    = [];
  NIRSSvol(ii).REFLEC_k_prf    = radar_nirss;
  NIRSSvol(ii).CLOUDS_k_prf    = cld_mask_nirss;
  NIRSSvol(ii).REFLEC_s_prf    = [];
  if str2num(input_advect) == 1 | (str2num(input_advect) == 3 & echo_in_domain == 0)
    NIRSSvol(ii).UVECTR_c_val    = uvec;
    NIRSSvol(ii).VVECTR_c_val    = vvec;
  elseif str2num(input_advect) == 2 | (str2num(input_advect) == 3 & echo_in_domain == 1)
    NIRSSvol(ii).UVECTR_c_val    = uvec(kcle_ctrec_lon,kcle_ctrec_lat);
    NIRSSvol(ii).VVECTR_c_val    = vvec(kcle_ctrec_lon,kcle_ctrec_lat);
  end
  
  %-------------------------------------------------------------------
  % Define a grid in the airport environment where a colorcoded NIRSS 
  % hazard will later be overlayed
  %-------------------------------------------------------------------
  grid_airport               = load('NIRSS_volumetric_KCLEapprdepart_map.txt', '-ascii');
  grid_airport_maxbox        = grid_airport;

  grid_airport2              = ones(haz_grid_dim,haz_grid_dim);
  grid_airport2(1:25,1:26)   = 6;
  grid_airport2(1:26,27:51)  = 7;
  grid_airport2(27:51,26:51) = 8;
  grid_airport2(26:51,1:25)  = 9;
  grid_airport2(14:26,14:26) = 2;
  grid_airport2(14:26,27:39) = 3;
  grid_airport2(27:39,26:39) = 4;
  grid_airport2(26:39,14:25) = 5;
  grid_airport2(24:28,24:28) = 1;

  %----------------------------------------------------------
  % edit array for tracking of packets advected by CTREC vectors 
  %----------------------------------------------------------

  % define hazard origin (KCLE) 
  %haz_advect_array(kcle_haz_y,kcle_haz_x,ii) = -1; 

  % find the change in u and v CTREC values at KCLE for this time
  if NIRSSvol(ii).UVECTR_c_val >= 0
    delta_u_km(ii) = ceil(NIRSSvol(ii).UVECTR_c_val*sec_betw_vols/1000);
  else
    delta_u_km(ii) = floor(NIRSSvol(ii).UVECTR_c_val*sec_betw_vols/1000);
  end

  if  NIRSSvol(ii).VVECTR_c_val >= 0
    delta_v_km(ii) = ceil(NIRSSvol(ii).VVECTR_c_val*sec_betw_vols/1000);
  else 
    delta_v_km(ii) = floor(NIRSSvol(ii).VVECTR_c_val*sec_betw_vols/1000);
  end

  % initialize ctrec vectors for iith NIRSSvol
  NIRSSvol(ii).UVECTR_c_val = [];
  NIRSSvol(ii).VVECTR_c_val = [];

  % check to see if delta u and v values from CTREC at KCLE are within reason
  if abs(delta_u_km(ii)) < 15 & abs(delta_v_km(ii)) < 15

    % determine the x and y coordinates for the latest hazard packet 
    haz_array_xind(ii) = kcle_haz_x + delta_u_km(ii); 
    haz_array_yind(ii) = kcle_haz_y - delta_v_km(ii); 
  
    %-----
    % move packets previous to iith packet 
    %   set up loop to find numbered index on grid and 
    %   move it from location in ii-1 
    %   to new advected location in iith array
    %-----
    nn = ii - 1;
  
    if nn > 0
      num_packets = length(NIRSSvol(nn).UVECTR_c_val);
  
      for mm = 1:num_packets
  
        % test to see if new coords are off the hazard grid
        if NIRSSvol(nn).VVECTR_c_val(ii-mm)-delta_v_km(ii) > haz_grid_dim | ...
           NIRSSvol(nn).UVECTR_c_val(ii-mm)+delta_u_km(ii) > haz_grid_dim | ...
           NIRSSvol(nn).VVECTR_c_val(ii-mm)-delta_v_km(ii) < 0 | ...
           NIRSSvol(nn).UVECTR_c_val(ii-mm)+delta_u_km(ii) < 0 
          display('Warning: Packet coordinate is beyond the hazard grid');
          NIRSSvol(ii).UVECTR_c_val = [NaN NIRSSvol(ii).UVECTR_c_val];  
          NIRSSvol(ii).VVECTR_c_val = [NaN NIRSSvol(ii).VVECTR_c_val];  
  
        else
          % find uvec and vvec at the locations of previous packets and
          % move them to their new locations within haz_advect_array(ii) 
          %haz_advect_array(NIRSSvol(nn).VVECTR_c_val(ii-mm)-delta_v_km(ii), ...
          %                 NIRSSvol(nn).UVECTR_c_val(ii-mm)+delta_u_km(ii),ii) = ii-mm; 
  
          if mm ~= ii
            NIRSSvol(ii).UVECTR_c_val = [NIRSSvol(nn).UVECTR_c_val(ii-mm)+delta_u_km(ii) ...
                                         NIRSSvol(ii).UVECTR_c_val];  
            NIRSSvol(ii).VVECTR_c_val = [NIRSSvol(nn).VVECTR_c_val(ii-mm)-delta_v_km(ii)  ...
                                         NIRSSvol(ii).VVECTR_c_val];  
          end  % end of if mm ~=ii
  
        end    % end of if NIRSSvol etc
  
      end      % end of for mm=1:num_packets
  
    end        % end of if nn>0

    %-----
    % move iith packet 
    %-----
    NIRSSvol(ii).UVECTR_c_val = [NIRSSvol(ii).UVECTR_c_val kcle_haz_x+delta_u_km(ii)];
    NIRSSvol(ii).VVECTR_c_val = [NIRSSvol(ii).VVECTR_c_val kcle_haz_y-delta_v_km(ii)];
  
    % change appropriate gridbox value to the time index (ii) for the newest packet
    %haz_advect_array(haz_array_yind(ii),haz_array_xind(ii),ii) = ii;

  else

    % since delta values are unreasonable, advect packets at previous rates

    % determine the x and y coordinates for the latest hazard packet 
    haz_array_xind(ii) = kcle_haz_x + delta_u_km(ii-1); 
    haz_array_yind(ii) = kcle_haz_y - delta_v_km(ii-1); 

    % set current deltas to previous delta values
    delta_u_km(ii)     = delta_u_km(ii-1);
    delta_v_km(ii)     = delta_v_km(ii-1);

    %-----
    % move packets previous to iith packet 
    %   set up loop to find numbered index on grid and 
    %   move it from location in ii-1 
    %   to new advected location in iith array
    %-----
    nn = ii - 1;
  
    if nn > 0
      num_packets = length(NIRSSvol(nn).UVECTR_c_val);
  
      for mm = 1:num_packets
  
        % test to see if new coords are off the hazard grid
        if NIRSSvol(nn).VVECTR_c_val(ii-mm)-delta_v_km(nn) > haz_grid_dim | ...
           NIRSSvol(nn).UVECTR_c_val(ii-mm)+delta_u_km(nn) > haz_grid_dim | ...
           NIRSSvol(nn).VVECTR_c_val(ii-mm)-delta_v_km(nn) < 0 | ...
           NIRSSvol(nn).UVECTR_c_val(ii-mm)+delta_u_km(nn) < 0 | ...
           isnan(NIRSSvol(nn).UVECTR_c_val(ii-mm)) | ...
           isnan(NIRSSvol(nn).VVECTR_c_val(ii-mm))  
           display('Warning: Packet coordinate is beyond the hazard grid');
           NIRSSvol(ii).UVECTR_c_val = [NaN NIRSSvol(ii).UVECTR_c_val];  
           NIRSSvol(ii).VVECTR_c_val = [NaN NIRSSvol(ii).VVECTR_c_val];  
  
        else
           % find uvec and vvec at the locations of previous packets and
           %   move them to their new locations within haz_advect_array(ii) 
           %haz_advect_array(NIRSSvol(nn).VVECTR_c_val(ii-mm)-delta_v_km(nn), ...
           %                 NIRSSvol(nn).UVECTR_c_val(ii-mm)+delta_u_km(nn),ii) = ii-mm; 

          if mm ~= ii
            NIRSSvol(ii).UVECTR_c_val = [NIRSSvol(nn).UVECTR_c_val(ii-mm)+delta_u_km(nn) ...
                                         NIRSSvol(ii).UVECTR_c_val];  
            NIRSSvol(ii).VVECTR_c_val = [NIRSSvol(nn).VVECTR_c_val(ii-mm)-delta_v_km(nn)  ...
                                         NIRSSvol(ii).VVECTR_c_val];  
          end  % end of if mm ~=ii

        end    % end of if NIRSSvol etc
      end      % end of for mm=1:num_packets
    end        % end of if nn>0

    %-----
    % move iith packet 
    %-----
    NIRSSvol(ii).UVECTR_c_val = [NIRSSvol(ii).UVECTR_c_val kcle_haz_x+delta_u_km(ii-1)];
    NIRSSvol(ii).VVECTR_c_val = [NIRSSvol(ii).VVECTR_c_val kcle_haz_y-delta_v_km(ii-1)];

    % change appropriate gridbox value to the time index (ii) for the newest packet
    %haz_advect_array(haz_array_yind(ii),haz_array_xind(ii),ii) = ii;
  
  end          % end of if abs(delta_u_km) < 20 etc

  %----------------------------------------------------------
  % use input from function 'forecast_discussion_text_parser.m' to  
  %   decide how to apply hazard values to terminal domain
  %----------------------------------------------------------
  if exist(num2str('interp_haz_to_radardomain_now')) == 1
  end

  if exist(num2str('interp_haz_to_radardomain_fut')) == 1
  end

  %----------------------------------------------------------
  % Check if any advected packets have advected into zones 6-9 
  %  if so, insert appropriate hazard values into zones
  %----------------------------------------------------------
  num_packets2 = length(NIRSSvol(ii).UVECTR_c_val);
  
  % oldest packet first (lowest #) so that oldest can be overwritten with newer
  for rr = 1:num_packets2

    for ss = 6:9

      if ~isnan(NIRSSvol(ii).VVECTR_c_val(rr)) | ~isnan(NIRSSvol(ii).UVECTR_c_val(rr))

        if grid_airport2(NIRSSvol(ii).VVECTR_c_val(rr), NIRSSvol(ii).UVECTR_c_val(rr)) == ss 

          NIRSSvol(ii).HAZARD_n_maxcol(ss) = NIRSSvol(rr).HAZARD_n_maxcol(1);
          NIRSSvol(ii).HAZARD_n_maxbox(ss) = NIRSSvol(rr).HAZARD_n_maxbox(1);
          ind_zone                         = find(grid_airport == ss);
          grid_airport(ind_zone)           = NIRSSvol(ii).HAZARD_n_maxcol(ss)/10;
          grid_airport_maxbox(ind_zone)    = NIRSSvol(ii).HAZARD_n_maxbox(ss)/10;
          clear ind_zone;
          disp(['haz added for zone ', num2str(ss), ' for packet number ', num2str(rr)]);

        end   % end of if grid_airport2==ss 

      end     % if ~isnan(etc)

    end       % end of for ss=6:9

  end         % end of for rr=1:num_packets2
  
  %----------------------------------------------------------
  % insert hazards into zones 1-5 of gridded airport space 
  %----------------------------------------------------------

  for pp = num_haz_zones:-1:1

    ind_zone = find(grid_airport == pp);
    if NIRSSvol(ii).HAZARD_n_maxcol(pp) == -1
      grid_airport(ind_zone) = -0.1;
    else
      grid_airport(ind_zone) = NIRSSvol(ii).HAZARD_n_maxcol(pp)/10;
    end
    if NIRSSvol(ii).HAZARD_n_maxbox(pp) == -1
      grid_airport_maxbox(ind_zone) = -0.1;
    else
      grid_airport_maxbox(ind_zone) = NIRSSvol(ii).HAZARD_n_maxbox(pp)/10;
    end
    clear ind_zone;
  end

  grid_airport        = grid_airport * 10;
  grid_airport_maxbox = grid_airport_maxbox * 10;

  %----------------------------------------------------------
  % Display some data to screen from iith volume 
  %----------------------------------------------------------
  if str2num(input_advect) == 1
    disp(['KCLE UVECTR and VVECTR = ', num2str(uvec), ', ' num2str(vvec)]); 
  elseif str2num(input_advect) == 2
    disp(['KCLE UVECTR and VVECTR = ', num2str(uvec(kcle_ctrec_lon,kcle_ctrec_lat)), ', ' ...
                                             num2str(vvec(kcle_ctrec_lon,kcle_ctrec_lat))]);
  end
  disp(['delta_u_km = ', num2str(delta_u_km)]);
  disp(['delta_v_km = ', num2str(delta_v_km)]);
  disp('NIRSS volumetric hazard structure packet ');
  NIRSSvol(ii)
  disp('----------')
  
  %----------------------------------------------------------
  % Current product plots for iith volume 
  %----------------------------------------------------------
  figure;
  subplot('position',[0.0 0.0 0.33 1.0])
  %usamap('Ohio');
  hold off;
  if str2num(input_advect) == 2
    contour(lon,lat,refl,[3,10]);
  end
  lat_new = 41.16:0.01:41.66;
  lon_new = -82.15:0.012:-81.55;
  ohiohi = shaperead('usastatehi', 'UseGeoCoords', true,...   
                     'Selector',{@(name) strcmpi(name,'Ohio'), 'Name'});
  geoshow(ohiohi, 'FaceColor','none');
  hold on;
  imagesc(lon_new, lat_new, flipud(grid_airport),[-1 8]);
  if str2num(input_advect) == 2
    contour(lon,lat,refl,[3,10]);
  end

  % location labels
  text(-82.4,41.5,'Lake Erie');
  text(-82.4,41.3,'Ohio');

  % NIRSS
  text(-81.855,41.409,'*');

  %% interior box to 50 km lines
  %%line('XData',[-82.150 -81.883],'YData',[41.415 41.415]);
  %%line('XData',[-81.820 -81.550],'YData',[41.405 41.405]);
  %%line('XData',[-81.855 -81.855],'YData',[41.160 41.383]);
  %%line('XData',[-81.845 -81.845],'YData',[41.660 41.437]);

  % interior box
  %line('XData',[-81.820 -81.820],'YData',[41.385 41.435]);
  %line('XData',[-81.880 -81.880],'YData',[41.385 41.435]);
  %line('XData',[-81.820 -81.880],'YData',[41.385 41.385]);
  %line('XData',[-81.820 -81.880],'YData',[41.435 41.435]);

  % runways
  %line('XData',[-81.830 -81.850],'YData',[41.410 41.398]); % 6R and 24L 
  %line('XData',[-81.836 -81.856],'YData',[41.412 41.401]); % 6L and 24R
  %line('XData',[-81.840 -81.822],'YData',[41.414 41.413]); % 28 and 10

  % 10 km lines
  %line('XData',[-82.000 -81.685],'YData',[41.275 41.275]);
  %line('XData',[-82.000 -81.685],'YData',[41.535 41.535]);
  %line('XData',[-82.000 -82.000],'YData',[41.275 41.535]);
  %line('XData',[-81.685 -81.685],'YData',[41.275 41.535]);

  % way points
  for zz=2:num_haz_zones 

    ind_zone = find(grid_airport2 == zz & grid_airport ~= 0);

    if mode(grid_airport(ind_zone)) == -1
      colorval = [1 1 1];
    else
      colorval = [0 0 0];
    end

    clear ind_zone;

    % use white marker if no hazard, black if yes hazard
    if zz == 2
      text(-81.938, 41.426, '*', 'Color', colorval);  % JUUBA 
    elseif zz == 3
    elseif zz == 4
      text(-81.720, 41.403, '*', 'Color', colorval);  % CINSO
    elseif zz == 5
      text(-81.890, 41.380, '*', 'Color', colorval);  % ZAPAN 
      text(-81.972, 41.328, '*', 'Color', colorval);  % SASCO
    elseif zz == 6
      text(-82.055, 41.443, '*', 'Color', colorval);  % FIMOL 
    elseif zz == 7
      text(-81.702, 41.485, '*', 'Color', colorval);  % JEDSU 
      text(-81.646, 41.517, '*', 'Color', colorval);  % ????? 
      text(-81.593, 41.547, '*', 'Color', colorval);  % ????? 
      text(-81.548, 41.580, '*', 'Color', colorval);  % PUDSE 
    elseif zz == 8
      text(-81.650, 41.400, '*', 'Color', colorval);  % ZANGI 
      text(-81.548, 41.396, '*', 'Color', colorval);  % EKUME 
    elseif zz == 9
      text(-82.079, 41.265, '*', 'Color', colorval);  % DUPAY 
      text(-82.125, 41.374, '*', 'Color', colorval);  % DRYER 
    end

  end
  
  title([NIRSSvol(ii).date_time{:} ' NIRSS volume max-in-col hazard']);

  subplot('position',[0.33 0.0 0.33 1.0])
  %usamap('Ohio');
  hold off;
  if str2num(input_advect) == 2
    contour(lon,lat,refl,[3,10]);
  end
  lat_new = 41.16:0.01:41.66;
  lon_new = -82.15:0.012:-81.55;
  ohiohi  = shaperead('usastatehi', 'UseGeoCoords', true,...   
                      'Selector',{@(name) strcmpi(name,'Ohio'), 'Name'});
  geoshow(ohiohi, 'FaceColor','none');
  hold on;
  imagesc(lon_new, lat_new, flipud(grid_airport_maxbox),[-1 8]);
  if str2num(input_advect) == 2
    contour(lon,lat,refl,[3,10]);
  end
  title([NIRSSvol(ii).date_time{:} ' NIRSS volume max-in-box hazard']);

  % interior box
  %line('XData',[-81.820 -81.820],'YData',[41.385 41.435]);
  %line('XData',[-81.880 -81.880],'YData',[41.385 41.435]);
  %line('XData',[-81.820 -81.880],'YData',[41.385 41.385]);
  %line('XData',[-81.820 -81.880],'YData',[41.435 41.435]);
  % runways
  %line('XData',[-81.830 -81.850],'YData',[41.410 41.398]); % 6R and 24L 
  %line('XData',[-81.836 -81.856],'YData',[41.412 41.401]); % 6L and 24R
  %line('XData',[-81.840 -81.822],'YData',[41.414 41.413]); % 28 and 10
  % NIRSS
  text(-81.855,41.409,'*');
  
  % way points
  for zz=2:num_haz_zones 

    ind_zone = find(grid_airport2 == zz & grid_airport ~= 0);

    if mode(grid_airport(ind_zone)) == -1
      colorval = [1 1 1];
    else
      colorval = [0 0 0];
    end

    clear ind_zone;

    % use white marker if no hazard, black if yes hazard
    if zz == 2
      text(-81.938, 41.426, '*', 'Color', colorval);  % JUUBA 
    elseif zz == 3
    elseif zz == 4
      text(-81.720, 41.403, '*', 'Color', colorval);  % CINSO
    elseif zz == 5
      text(-81.890, 41.380, '*', 'Color', colorval);  % ZAPAN 
      text(-81.972, 41.328, '*', 'Color', colorval);  % SASCO
    elseif zz == 6
      text(-82.055, 41.443, '*', 'Color', colorval);  % FIMOL 
    elseif zz == 7
      text(-81.702, 41.485, '*', 'Color', colorval);  % JEDSU 
      text(-81.646, 41.517, '*', 'Color', colorval);  % ????? 
      text(-81.593, 41.547, '*', 'Color', colorval);  % ????? 
      text(-81.548, 41.580, '*', 'Color', colorval);  % PUDSE 
    elseif zz == 8
      text(-81.650, 41.400, '*', 'Color', colorval);  % ZANGI 
      text(-81.548, 41.396, '*', 'Color', colorval);  % EKUME 
    elseif zz == 9
      text(-82.079, 41.265, '*', 'Color', colorval);  % DUPAY 
      text(-82.125, 41.374, '*', 'Color', colorval);  % DRYER 
    end

  end
  
  % create sample hazard plot for control tower operations
  %grid_airport grid_airport2;

  subplot('position',[0.8 0.5 0.05 0.05])
      axis square; axis off;
      %ind = find(grid_airport2 == 1);
      %if grid_airport(ind) == 
      %
      test=ones(2,2)*0.5;
      imagesc(test,[0 1]);
      axis square; axis off;
      title('KCLE');
  subplot('position',[0.75 0.51 0.05 0.05])
      axis square; axis off;
      test=ones(2,2)*0.5;
      imagesc(test,[0 1]);
      axis square; axis off;
      title('JUUBA');
  subplot('position',[0.70 0.52 0.05 0.05])
      axis square; axis off;
      test=ones(2,2)*0.5;
      imagesc(test,[0 1]);
      axis square; axis off;
      title('FIMOL');
  subplot('position',[0.85 0.6 0.05 0.05])
      axis square; axis off;
      test=ones(2,2)*0.5;
      imagesc(test,[0 1]);
      axis square; axis off;
      title('24 APRO');
  subplot('position',[0.9 0.7 0.05 0.05])
      axis square; axis off;
      test=ones(2,2)*0.5;
      imagesc(test,[0 1]);
      axis square; axis off;
      title('JEDSU');
  subplot('position',[0.85 0.9 0.05 0.05])
      axis square; axis off;
      test=ones(2,2)*0.5;
      imagesc(test,[0 1]);
      axis square; axis off;
      title('PASLE-H');
  subplot('position',[0.95 0.6 0.05 0.05])
      axis square; axis off;
      test=ones(2,2)*0.5;
      imagesc(test,[0 1]);
      axis square; axis off;
      title('LEBRN-H');
  subplot('position',[0.85 0.49 0.05 0.05])
      axis square; axis off;
      test=ones(2,2)*0.5;
      imagesc(test,[0 1]);
      axis square; axis off;
      title('CINSO');
  subplot('position',[0.90 0.48 0.05 0.05])
      axis square; axis off;
      test=ones(2,2)*0.5;
      imagesc(test,[0 1]);
      axis square; axis off;
      title('ZANGI');
  subplot('position',[0.95 0.47 0.05 0.05])
      axis square; axis off;
      test=ones(2,2)*0.5;
      imagesc(test,[0 1]);
      axis square; axis off;
      title('EKUME');
  subplot('position',[0.95 0.3 0.05 0.05])
      axis square; axis off;
      test=ones(2,2)*0.5;
      imagesc(test,[0 1]);
      axis square; axis off;
      title('IXORE-H');
  subplot('position',[0.85 0.2 0.05 0.05])
      axis square; axis off;
      test=ones(2,2)*0.5;
      imagesc(test,[0 1]);
      axis square; axis off;
      title('WOMGO-H');
  subplot('position',[0.68 0.45 0.05 0.05])
      axis square; axis off;
      test=ones(2,2)*0.5;
      imagesc(test,[0 1]);
      axis square; axis off;
      title('DRYER-H');
  subplot('position',[0.75 0.4 0.05 0.05])
      axis square; axis off;
      test=ones(2,2)*0.5;
      imagesc(test,[0 1]);
      axis square; axis off;
      title('ZAPAN/SASCO');
  subplot('position',[0.7 0.3 0.05 0.05])
      axis square; axis off;
      test=ones(2,2)*0.5;
      imagesc(test,[0 1]);
      axis square; axis off;
      title('DUPAY/SERLE');

  %----------------------------------------------------------
  % output plot to daily case directory 
  %----------------------------------------------------------
  cd /d1/serke/projects/NIRSS_NASA/image/volumetric_product/;
  try
    string100 = ['cd ' num2str(yyyy1) num2str(mon1) num2str(day1)];
    eval(string100);
  catch
    string101= ['mkdir ' num2str(yyyy1) num2str(mon1) num2str(day1)];
    eval(string101);
    string100 = ['cd ' num2str(yyyy1) num2str(mon1) num2str(day1)];
    eval(string100);
  end
  saveas(gcf,[yyyy1 mon1 day1 '_' desired_hhmm{ii} '_maxinbox.png'],'png');
  cd /d1/serke/projects/NIRSS_NASA/code/mat;

  %----------------------------------------------------------
  % Older or commented out plots 
  %----------------------------------------------------------

    
end   % end of for i=1:num_times

    
