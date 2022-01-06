%------------------------------------
% NAME:    compare_NIRSSfusion_to_CIPPIREP_v7up
%
% PURPOSE: loop thru NIRSS fusion day files
%            read day mat file
%            load severity fields into array
%            plot colocated NIRSS, CIP and PIREP severities
%
% INPUTS: 
%
% OUTPUTS: 
%
% CREATED: 4.20.2010 dserke
%
% MODIFICATIONS:
%
% EXAMPLE: matlab -nodesktop
%          compare_NIRSSfusion_to_CIPPIREP
%
%------------------------------------

%------------------------------------
% update paths 
%------------------------------------
addpath ~/cvs/matlab/raptoolboxes/String/mfiles;

%------------------------------------
% keywords 
%------------------------------------
compute_res = 0;

%------------------------------------
% define looping dates 
%------------------------------------
start_date_yyyymmdd = num2str(20081104);
%start_date_yyyymmdd = num2str(20100321);
end_date_yyyymmdd   = num2str(20100330);
start_year          = start_date_yyyymmdd(1:4);
start_mon           = start_date_yyyymmdd(5:6);
start_day           = start_date_yyyymmdd(7:8);
end_year            = end_date_yyyymmdd(1:4);
end_mon             = end_date_yyyymmdd(5:6);
end_day             = end_date_yyyymmdd(7:8); 
start_date          = datenum(str2num(start_year),str2num(start_mon),str2num(start_day));
end_date            = datenum(str2num(end_year),str2num(end_mon),str2num(end_day));

%keyboard;

%------------------------------------
% define working paths
%------------------------------------
NIRSS_CIP_PIREP_input_path = '/d2/nirss/nasa_glen/NIRSS_CIP_PIREP_output/data';
%NIRSS_CIP_PIREP_input_path = '/d2/nirss/nasa_glen/NIRSS_CIP_PIREP_output/data/old_sevmatchup';

%------------------------------------
% get list of matchup data files  
%------------------------------------
cd_string = ['cd ' NIRSS_CIP_PIREP_input_path];
eval(cd_string);
unix('ls *sevmatchup.11.mat > matchup_files.txt');  
load('matchup_files.txt');

%keyboard;

%-----------------------------------------------------------
% loop through all desired dates (NIRSS, PIREP and CIP data
%-----------------------------------------------------------
date_dry_array                       = [];
date_wet_array                       = [];
pirep_hh_dry_array                   = [];
pirep_hh_wet_array                   = [];
pirep_mm_dry_array                   = [];
pirep_mm_wet_array                   = [];
pirep_alt_dry_array                  = [];
pirep_alt_wet_array                  = [];
pirep_sev_dry_array                  = [];
pirep_sev_wet_array                  = [];
cipice_sev_PIREP_dry_array           = [];
cipice_sev_PIREP_wet_array           = [];
nirss_hazard_PIREP_max_array         = [];
nirss_hazard_PIREP_mean_array        = [];
nirss_hazard_PIREP_dry_closest_array = [];
nirss_hazard_PIREP_wet_closest_array = [];
cip_warn_vol_array                   = [];
nirss_warn_vol_array                 = [];

day_num = 0;  % start a count of total days processed
for i = start_date:end_date

  day_num = day_num + 1;
  matchup_filename = datestr(i,'yyyymmddTHHMMSS');
  matchup_filedate = str2num(matchup_filename(1:8));
  disp('************************************');
  disp(['looking for file:     date= ' num2str(matchup_filedate)]);

  % check to see if matchup file exists
  matchup_file_ind = find(matchup_files == matchup_filedate); 
  if size(matchup_file_ind,1) == 0
    disp('  error: no matchup file found!');
  else
    % if exists then load mat file
    load([num2str(matchup_filedate) '_sevmatchup.11.mat']);
    
    % load daily arrays into larger array
    %disp(date);
    %disp(size(date));
    %disp(size(pirep_sev));
    %disp(size(nirss_hazard_PIREP_closest));
    %disp(size(nirss_hazard_PIREP_mean));
    date_dry_array                      = [date_dry_array date_dry];
    pirep_hh_dry_array                  = [pirep_hh_dry_array pirep_hh_dry];
    pirep_mm_dry_array                  = [pirep_mm_dry_array pirep_mm_dry];
    pirep_alt_dry_array                 = [pirep_alt_dry_array pirep_alt_dry];
    pirep_sev_dry_array                 = [pirep_sev_dry_array pirep_sev_dry];
    cipice_sev_PIREP_dry_array          = [cipice_sev_PIREP_dry_array cipice_sev_PIREP_dry];
    nirss_hazard_PIREP_max_array        = [nirss_hazard_PIREP_max_array nirss_hazard_PIREP_max];
    nirss_hazard_PIREP_mean_array       = [nirss_hazard_PIREP_mean_array nirss_hazard_PIREP_mean];
    nirss_hazard_PIREP_dry_closest_array = [nirss_hazard_PIREP_dry_closest_array nirss_hazard_PIREP_dry_closest];

    date_wet_array                      = [date_wet_array date_wet];
    pirep_hh_wet_array                  = [pirep_hh_wet_array pirep_hh_wet];
    pirep_mm_wet_array                  = [pirep_mm_wet_array pirep_mm_wet];
    pirep_alt_wet_array                 = [pirep_alt_wet_array pirep_alt_wet];
    pirep_sev_wet_array                 = [pirep_sev_wet_array pirep_sev_wet];
    cipice_sev_PIREP_wet_array          = [cipice_sev_PIREP_wet_array cipice_sev_PIREP_wet];
    nirss_hazard_PIREP_max_array        = [nirss_hazard_PIREP_max_array nirss_hazard_PIREP_max];
    nirss_hazard_PIREP_mean_array       = [nirss_hazard_PIREP_mean_array nirss_hazard_PIREP_mean];
    nirss_hazard_PIREP_wet_closest_array = [nirss_hazard_PIREP_wet_closest_array nirss_hazard_PIREP_wet_closest];
    cip_warn_vol_array                  = [cip_warn_vol_array cip_warn_vol_day];
    nirss_warn_vol_array                = [nirss_warn_vol_array nirss_warn_vol_day];

  end

  %keyboard;

end      % end of for i=start_date:end_date

test=find(pirep_sev_dry_array == -1);
pirep_sev_dry_array(test) = 0;

%--------------------------------------------------------------------------
% compute mean product warning volumes for all days run through 'compare.m 
%--------------------------------------------------------------------------
cip_warn_vol_average = mean(cip_warn_vol_array);
cip_warn_vol_average
nirss_warn_vol_average = mean(nirss_warn_vol_array);
nirss_warn_vol_average

%---------------------------------
% loop through accumulated fields and bin the values 
%---------------------------------
%bin_nirss_clo_0 = 0;
%bin_nirss_clo_1 = 0;
%bin_nirss_clo_2 = 0;
%bin_nirss_clo_3 = 0;
%bin_nirss_clo_4 = 0;
%bin_nirss_clo_5 = 0;
%bin_nirss_clo_6 = 0;
%bin_nirss_clo_7 = 0;
%bin_nirss_clo_8 = 0;
%for jj = 1:length(nirss_hazard_PIREP_closest_array)
%  if nirss_hazard_PIREP_closest_array(jj) < 0.5
%    bin_nirss_clo_0 = bin_nirss_clo_0 + 1;
%  elseif nirss_hazard_PIREP_closest_array(jj) < 1.5 & nirss_hazard_PIREP_closest_array(jj) >= 0.5
%    bin_nirss_clo_1 = bin_nirss_clo_1 + 1;
%  elseif nirss_hazard_PIREP_closest_array(jj) < 2.5 & nirss_hazard_PIREP_closest_array(jj) >= 1.5
%    bin_nirss_clo_2 = bin_nirss_clo_2 + 1;
%  elseif nirss_hazard_PIREP_closest_array(jj) < 3.5 & nirss_hazard_PIREP_closest_array(jj) >= 2.5
%    bin_nirss_clo_3 = bin_nirss_clo_3 + 1;
%  elseif nirss_hazard_PIREP_closest_array(jj) < 4.5 & nirss_hazard_PIREP_closest_array(jj) >= 3.5
%    bin_nirss_clo_4 = bin_nirss_clo_4 + 1;
%  elseif nirss_hazard_PIREP_closest_array(jj) < 5.5 & nirss_hazard_PIREP_closest_array(jj) >= 4.5
%    bin_nirss_clo_5 = bin_nirss_clo_5 + 1;
%  elseif nirss_hazard_PIREP_closest_array(jj) < 6.5 & nirss_hazard_PIREP_closest_array(jj) >= 5.5
%    bin_nirss_clo_6 = bin_nirss_clo_6 + 1;
%  elseif nirss_hazard_PIREP_closest_array(jj) < 7.5 & nirss_hazard_PIREP_closest_array(jj) >= 6.5
%    bin_nirss_clo_7 = bin_nirss_clo_7 + 1;
%  elseif nirss_hazard_PIREP_closest_array(jj) >= 7.5
%    bin_nirss_clo_8 = bin_nirss_clo_8 + 1;
%  end
%end
%bins_nirss_closest_summary = [bin_nirss_clo_0 bin_nirss_clo_1 bin_nirss_clo_2 bin_nirss_clo_3 bin_nirss_clo_4 bin_nirss_clo_5 bin_nirss_clo_6 bin_nirss_clo_7 bin_nirss_clo_8];

%bin_nirss_max_0 = 0;
%bin_nirss_max_1 = 0;
%bin_nirss_max_2 = 0;
%bin_nirss_max_3 = 0;
%bin_nirss_max_4 = 0;
%bin_nirss_max_5 = 0;
%bin_nirss_max_6 = 0;
%bin_nirss_max_7 = 0;
%bin_nirss_max_8 = 0;
%for jj = 1:length(nirss_hazard_PIREP_max_array)
%  if nirss_hazard_PIREP_max_array(jj) < 0.5
%    bin_nirss_max_0 = bin_nirss_max_0 + 1;
%  elseif nirss_hazard_PIREP_max_array(jj) < 1.5 & nirss_hazard_PIREP_max_array(jj) >= 0.5
%    bin_nirss_max_1 = bin_nirss_max_1 + 1;
%  elseif nirss_hazard_PIREP_max_array(jj) < 2.5 & nirss_hazard_PIREP_max_array(jj) >= 1.5
%    bin_nirss_max_2 = bin_nirss_max_2 + 1;
%  elseif nirss_hazard_PIREP_max_array(jj) < 3.5 & nirss_hazard_PIREP_max_array(jj) >= 2.5
%    bin_nirss_max_3 = bin_nirss_max_3 + 1;
%  elseif nirss_hazard_PIREP_max_array(jj) < 4.5 & nirss_hazard_PIREP_max_array(jj) >= 3.5
%    bin_nirss_max_4 = bin_nirss_max_4 + 1;
%  elseif nirss_hazard_PIREP_max_array(jj) < 5.5 & nirss_hazard_PIREP_max_array(jj) >= 4.5
%    bin_nirss_max_5 = bin_nirss_max_5 + 1;
%  elseif nirss_hazard_PIREP_max_array(jj) < 6.5 & nirss_hazard_PIREP_max_array(jj) >= 5.5
%    bin_nirss_max_6 = bin_nirss_max_6 + 1;
%  elseif nirss_hazard_PIREP_max_array(jj) < 7.5 & nirss_hazard_PIREP_max_array(jj) >= 6.5
%    bin_nirss_max_7 = bin_nirss_max_7 + 1;
%  elseif nirss_hazard_PIREP_max_array(jj) >= 7.5
%    bin_nirss_max_8 = bin_nirss_max_8 + 1;
%  end
%end
%bins_nirss_max_summary = [bin_nirss_max_0 bin_nirss_max_1 bin_nirss_max_2 bin_nirss_max_3 bin_nirss_max_4 bin_nirss_max_5 bin_nirss_max_6 bin_nirss_max_7 bin_nirss_max_8];

%bin_pirep_0 = 0;
%bin_pirep_1 = 0;
%bin_pirep_2 = 0;
%bin_pirep_3 = 0;
%bin_pirep_4 = 0;
%bin_pirep_5 = 0;
%bin_pirep_6 = 0;
%bin_pirep_7 = 0;
%bin_pirep_8 = 0;
%for jj = 1:length(pirep_sev_array)
%  if pirep_sev_array(jj) < 0.5
%    bin_pirep_0 = bin_pirep_0 + 1;
%  elseif pirep_sev_array(jj) < 1.5 & pirep_sev_array(jj) >= 0.5
%    bin_pirep_1 = bin_pirep_1 + 1;
%  elseif pirep_sev_array(jj) < 2.5 & pirep_sev_array(jj) >= 1.5
%    bin_pirep_2 = bin_pirep_2 + 1;
%  elseif pirep_sev_array(jj) < 3.5 & pirep_sev_array(jj) >= 2.5
%    bin_pirep_3 = bin_pirep_3 + 1;
%  elseif pirep_sev_array(jj) < 4.5 & pirep_sev_array(jj) >= 3.5
%    bin_pirep_4 = bin_pirep_4 + 1;
%  elseif pirep_sev_array(jj) < 5.5 & pirep_sev_array(jj) >= 4.5
%    bin_pirep_5 = bin_pirep_5 + 1;
%  elseif pirep_sev_array(jj) < 6.5 & pirep_sev_array(jj) >= 5.5
%    bin_pirep_6 = bin_pirep_6 + 1;
%  elseif pirep_sev_array(jj) < 7.5 & pirep_sev_array(jj) >= 6.5
%    bin_pirep_7 = bin_pirep_7 + 1;
%  elseif pirep_sev_array(jj) >= 7.5
%    bin_pirep_8 = bin_pirep_8 + 1;
%  end
%end
%bins_pirep = [bin_pirep_0 bin_pirep_1 bin_pirep_2 bin_pirep_3 bin_pirep_4 bin_pirep_5 bin_pirep_6 bin_pirep_7 bin_pirep_8];

%---------------------------------------------
% sum up bins for nirss column MEAN versus pireps
%---------------------------------------------
bin_p0_n0 = 0;bin_p1_n0 = 0;bin_p2_n0 = 0;bin_p3_n0 = 0;bin_p4_n0 = 0;bin_p5_n0 = 0;bin_p6_n0 = 0;bin_p7_n0 = 0;bin_p8_n0 = 0;
bin_p0_n1 = 0;bin_p1_n1 = 0;bin_p2_n1 = 0;bin_p3_n1 = 0;bin_p4_n1 = 0;bin_p5_n1 = 0;bin_p6_n1 = 0;bin_p7_n1 = 0;bin_p8_n1 = 0;
bin_p0_n2 = 0;bin_p1_n2 = 0;bin_p2_n2 = 0;bin_p3_n2 = 0;bin_p4_n2 = 0;bin_p5_n2 = 0;bin_p6_n2 = 0;bin_p7_n2 = 0;bin_p8_n2 = 0;
bin_p0_n3 = 0;bin_p1_n3 = 0;bin_p2_n3 = 0;bin_p3_n3 = 0;bin_p4_n3 = 0;bin_p5_n3 = 0;bin_p6_n3 = 0;bin_p7_n3 = 0;bin_p8_n3 = 0;
bin_p0_n4 = 0;bin_p1_n4 = 0;bin_p2_n4 = 0;bin_p3_n4 = 0;bin_p4_n4 = 0;bin_p5_n4 = 0;bin_p6_n4 = 0;bin_p7_n4 = 0;bin_p8_n4 = 0;
bin_p0_n5 = 0;bin_p1_n5 = 0;bin_p2_n5 = 0;bin_p3_n5 = 0;bin_p4_n5 = 0;bin_p5_n5 = 0;bin_p6_n5 = 0;bin_p7_n5 = 0;bin_p8_n5 = 0;
bin_p0_n6 = 0;bin_p1_n6 = 0;bin_p2_n6 = 0;bin_p3_n6 = 0;bin_p4_n6 = 0;bin_p5_n6 = 0;bin_p6_n6 = 0;bin_p7_n6 = 0;bin_p8_n6 = 0;
bin_p0_n7 = 0;bin_p1_n7 = 0;bin_p2_n7 = 0;bin_p3_n7 = 0;bin_p4_n7 = 0;bin_p5_n7 = 0;bin_p6_n7 = 0;bin_p7_n7 = 0;bin_p8_n7 = 0;
bin_p0_n8 = 0;bin_p1_n8 = 0;bin_p2_n8 = 0;bin_p3_n8 = 0;bin_p4_n8 = 0;bin_p5_n8 = 0;bin_p6_n8 = 0;bin_p7_n8 = 0;bin_p8_n8 = 0;


for jj = 1:length(pirep_sev_dry_array)
  if pirep_sev_dry_array(jj) < 0.5
    if nirss_hazard_PIREP_mean_array(jj) < 0.5
      bin_p0_n0 = bin_p0_n0 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 1.5 & nirss_hazard_PIREP_mean_array(jj) >= 0.5
      bin_p0_n1 = bin_p0_n1 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 2.5 & nirss_hazard_PIREP_mean_array(jj) >= 1.5
      bin_p0_n2 = bin_p0_n2 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 3.5 & nirss_hazard_PIREP_mean_array(jj) >= 2.5
      bin_p0_n3 = bin_p0_n3 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 4.5 & nirss_hazard_PIREP_mean_array(jj) >= 3.5
      bin_p0_n4 = bin_p0_n4 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 5.5 & nirss_hazard_PIREP_mean_array(jj) >= 4.5
      bin_p0_n5 = bin_p0_n5 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 6.5 & nirss_hazard_PIREP_mean_array(jj) >= 5.5
      bin_p0_n6 = bin_p0_n6 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 7.5 & nirss_hazard_PIREP_mean_array(jj) >= 6.5
      bin_p0_n7 = bin_p0_n7 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) >= 7.5
      bin_p0_n8 = bin_p0_n8 + 1;
    end
  elseif pirep_sev_dry_array(jj) < 1.5 & pirep_sev_dry_array(jj) >= 0.5
    if nirss_hazard_PIREP_mean_array(jj) < 0.5
      bin_p1_n0 = bin_p1_n0 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 1.5 & nirss_hazard_PIREP_mean_array(jj) >= 0.5
      bin_p1_n1 = bin_p1_n1 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 2.5 & nirss_hazard_PIREP_mean_array(jj) >= 1.5
      bin_p1_n2 = bin_p1_n2 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 3.5 & nirss_hazard_PIREP_mean_array(jj) >= 2.5
      bin_p1_n3 = bin_p1_n3 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 4.5 & nirss_hazard_PIREP_mean_array(jj) >= 3.5
      bin_p1_n4 = bin_p1_n4 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 5.5 & nirss_hazard_PIREP_mean_array(jj) >= 4.5
      bin_p1_n5 = bin_p1_n5 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 6.5 & nirss_hazard_PIREP_mean_array(jj) >= 5.5
      bin_p1_n6 = bin_p1_n6 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 7.5 & nirss_hazard_PIREP_mean_array(jj) >= 6.5
      bin_p1_n7 = bin_p1_n7 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) >= 7.5
      bin_p1_n8 = bin_p1_n8 + 1;
    end
  elseif pirep_sev_dry_array(jj) < 2.5 & pirep_sev_dry_array(jj) >= 1.5
    if nirss_hazard_PIREP_mean_array(jj) < 0.5
      bin_p2_n0 = bin_p2_n0 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 1.5 & nirss_hazard_PIREP_mean_array(jj) >= 0.5
      bin_p2_n1 = bin_p2_n1 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 2.5 & nirss_hazard_PIREP_mean_array(jj) >= 1.5
      bin_p2_n2 = bin_p2_n2 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 3.5 & nirss_hazard_PIREP_mean_array(jj) >= 2.5
      bin_p2_n3 = bin_p2_n3 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 4.5 & nirss_hazard_PIREP_mean_array(jj) >= 3.5
      bin_p2_n4 = bin_p2_n4 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 5.5 & nirss_hazard_PIREP_mean_array(jj) >= 4.5
      bin_p2_n5 = bin_p2_n5 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 6.5 & nirss_hazard_PIREP_mean_array(jj) >= 5.5
      bin_p2_n6 = bin_p2_n6 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 7.5 & nirss_hazard_PIREP_mean_array(jj) >= 6.5
      bin_p2_n7 = bin_p2_n7 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) >= 7.5
      bin_p2_n8 = bin_p2_n8 + 1;
    end
  elseif pirep_sev_dry_array(jj) < 3.5 & pirep_sev_dry_array(jj) >= 2.5
    if nirss_hazard_PIREP_mean_array(jj) < 0.5
      bin_p3_n0 = bin_p3_n0 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 1.5 & nirss_hazard_PIREP_mean_array(jj) >= 0.5
      bin_p3_n1 = bin_p3_n1 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 2.5 & nirss_hazard_PIREP_mean_array(jj) >= 1.5
      bin_p3_n2 = bin_p3_n2 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 3.5 & nirss_hazard_PIREP_mean_array(jj) >= 2.5
      bin_p3_n3 = bin_p3_n3 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 4.5 & nirss_hazard_PIREP_mean_array(jj) >= 3.5
      bin_p3_n4 = bin_p3_n4 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 5.5 & nirss_hazard_PIREP_mean_array(jj) >= 4.5
      bin_p3_n5 = bin_p3_n5 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 6.5 & nirss_hazard_PIREP_mean_array(jj) >= 5.5
      bin_p3_n6 = bin_p3_n6 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 7.5 & nirss_hazard_PIREP_mean_array(jj) >= 6.5
      bin_p3_n7 = bin_p3_n7 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) >= 7.5
      bin_p3_n8 = bin_p3_n8 + 1;
    end
  elseif pirep_sev_dry_array(jj) < 4.5 & pirep_sev_dry_array(jj) >= 3.5
    if nirss_hazard_PIREP_mean_array(jj) < 0.5
      bin_p1_n0 = bin_p4_n0 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 1.5 & nirss_hazard_PIREP_mean_array(jj) >= 0.5
      bin_p4_n1 = bin_p4_n1 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 2.5 & nirss_hazard_PIREP_mean_array(jj) >= 1.5
      bin_p4_n2 = bin_p4_n2 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 3.5 & nirss_hazard_PIREP_mean_array(jj) >= 2.5
      bin_p4_n3 = bin_p4_n3 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 4.5 & nirss_hazard_PIREP_mean_array(jj) >= 3.5
      bin_p4_n4 = bin_p4_n4 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 5.5 & nirss_hazard_PIREP_mean_array(jj) >= 4.5
      bin_p4_n5 = bin_p4_n5 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 6.5 & nirss_hazard_PIREP_mean_array(jj) >= 5.5
      bin_p4_n6 = bin_p4_n6 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 7.5 & nirss_hazard_PIREP_mean_array(jj) >= 6.5
      bin_p4_n7 = bin_p4_n7 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) >= 7.5
      bin_p4_n8 = bin_p4_n8 + 1;
    end
  elseif pirep_sev_dry_array(jj) < 5.5 & pirep_sev_dry_array(jj) >= 4.5
    if nirss_hazard_PIREP_mean_array(jj) < 0.5
      bin_p5_n0 = bin_p5_n0 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 1.5 & nirss_hazard_PIREP_mean_array(jj) >= 0.5
      bin_p5_n1 = bin_p5_n1 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 2.5 & nirss_hazard_PIREP_mean_array(jj) >= 1.5
      bin_p5_n2 = bin_p5_n2 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 3.5 & nirss_hazard_PIREP_mean_array(jj) >= 2.5
      bin_p5_n3 = bin_p5_n3 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 4.5 & nirss_hazard_PIREP_mean_array(jj) >= 3.5
      bin_p5_n4 = bin_p5_n4 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 5.5 & nirss_hazard_PIREP_mean_array(jj) >= 4.5
      bin_p5_n5 = bin_p5_n5 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 6.5 & nirss_hazard_PIREP_mean_array(jj) >= 5.5
      bin_p5_n6 = bin_p5_n6 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 7.5 & nirss_hazard_PIREP_mean_array(jj) >= 6.5
      bin_p5_n7 = bin_p5_n7 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) >= 7.5
      bin_p5_n8 = bin_p5_n8 + 1;
    end
  elseif pirep_sev_dry_array(jj) < 6.5 & pirep_sev_dry_array(jj) >= 5.5
    if nirss_hazard_PIREP_mean_array(jj) < 0.5
      bin_p6_n0 = bin_p6_n0 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 1.5 & nirss_hazard_PIREP_mean_array(jj) >= 0.5
      bin_p6_n1 = bin_p6_n1 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 2.5 & nirss_hazard_PIREP_mean_array(jj) >= 1.5
      bin_p6_n2 = bin_p6_n2 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 3.5 & nirss_hazard_PIREP_mean_array(jj) >= 2.5
      bin_p6_n3 = bin_p6_n3 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 4.5 & nirss_hazard_PIREP_mean_array(jj) >= 3.5
      bin_p6_n4 = bin_p6_n4 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 5.5 & nirss_hazard_PIREP_mean_array(jj) >= 4.5
      bin_p6_n5 = bin_p6_n5 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 6.5 & nirss_hazard_PIREP_mean_array(jj) >= 5.5
      bin_p6_n6 = bin_p6_n6 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 7.5 & nirss_hazard_PIREP_mean_array(jj) >= 6.5
      bin_p6_n7 = bin_p6_n7 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) >= 7.5
      bin_p6_n8 = bin_p6_n8 + 1;
    end
  elseif pirep_sev_dry_array(jj) < 7.5 & pirep_sev_dry_array(jj) >= 6.5
    if nirss_hazard_PIREP_mean_array(jj) < 0.5
      bin_p7_n0 = bin_p7_n0 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 1.5 & nirss_hazard_PIREP_mean_array(jj) >= 0.5
      bin_p7_n1 = bin_p7_n1 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 2.5 & nirss_hazard_PIREP_mean_array(jj) >= 1.5
      bin_p7_n2 = bin_p7_n2 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 3.5 & nirss_hazard_PIREP_mean_array(jj) >= 2.5
      bin_p7_n3 = bin_p7_n3 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 4.5 & nirss_hazard_PIREP_mean_array(jj) >= 3.5
      bin_p7_n4 = bin_p7_n4 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 5.5 & nirss_hazard_PIREP_mean_array(jj) >= 4.5
      bin_p7_n5 = bin_p7_n5 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 6.5 & nirss_hazard_PIREP_mean_arrray(jj) >= 5.5
      bin_p7_n6 = bin_p7_n6 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 7.5 & nirss_hazard_PIREP_mean_array(jj) >= 6.5
      bin_p7_n7 = bin_p7_n7 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) >= 7.5
      bin_p7_n8 = bin_p7_n8 + 1;
    end
  elseif pirep_sev_dry_array(jj) < 8.5 & pirep_sev_dry_array(jj) >= 7.5
    if nirss_hazard_PIREP_mean_array(jj) < 0.5
      bin_p8_n0 = bin_p8_n0 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 1.5 & nirss_hazard_PIREP_mean_array(jj) >= 0.5
      bin_p8_n1 = bin_p8_n1 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 2.5 & nirss_hazard_PIREP_mean_array(jj) >= 1.5
      bin_p8_n2 = bin_p8_n2 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 3.5 & nirss_hazard_PIREP_mean_array(jj) >= 2.5
      bin_p8_n3 = bin_p8_n3 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 4.5 & nirss_hazard_PIREP_mean_array(jj) >= 3.5
      bin_p8_n4 = bin_p8_n4 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 5.5 & nirss_hazard_PIREP_mean_array(jj) >= 4.5
      bin_p8_n5 = bin_p8_n5 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 6.5 & nirss_hazard_PIREP_mean_array(jj) >= 5.5
      bin_p8_n6 = bin_p8_n6 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) < 7.5 & nirss_hazard_PIREP_mean_array(jj) >= 6.5
      bin_p8_n7 = bin_p8_n7 + 1;
    elseif nirss_hazard_PIREP_mean_array(jj) >= 7.5
      bin_p8_n8 = bin_p8_n8 + 1;
    end
  end
end
bins_p_vs_nmean = [bin_p0_n8 bin_p1_n8 bin_p2_n8 bin_p3_n8 bin_p4_n8 bin_p5_n8 bin_p6_n8 bin_p7_n8 bin_p8_n8; bin_p0_n7 bin_p1_n7 bin_p2_n7 bin_p3_n7 bin_p4_n7 bin_p5_n7 bin_p6_n7 bin_p7_n7 bin_p8_n7; bin_p0_n6 bin_p1_n6 bin_p2_n6 bin_p3_n6 bin_p4_n6 bin_p5_n6 bin_p6_n6 bin_p7_n6 bin_p8_n6; bin_p0_n5 bin_p1_n5 bin_p2_n5 bin_p3_n5 bin_p4_n5 bin_p5_n5 bin_p6_n5 bin_p7_n5 bin_p8_n5; bin_p0_n4 bin_p1_n4 bin_p2_n4 bin_p3_n4 bin_p4_n4 bin_p5_n4 bin_p6_n4 bin_p7_n4 bin_p8_n4; bin_p0_n3 bin_p1_n3 bin_p2_n3 bin_p3_n3 bin_p4_n3 bin_p5_n3 bin_p6_n3 bin_p7_n3 bin_p8_n3; bin_p0_n2 bin_p1_n2 bin_p2_n2 bin_p3_n2 bin_p4_n2 bin_p5_n2 bin_p6_n2 bin_p7_n2 bin_p8_n2; bin_p0_n1 bin_p1_n1 bin_p2_n1 bin_p3_n1 bin_p4_n1 bin_p5_n1 bin_p6_n1 bin_p7_n1 bin_p8_n1; bin_p0_n0 bin_p1_n0 bin_p2_n0 bin_p3_n0 bin_p4_n0 bin_p5_n0 bin_p6_n0 bin_p7_n0 bin_p8_n0];

%---------------------------------------------
% sum up bins for nirss dry closest versus pireps
%---------------------------------------------
bin_p0_n0c = 0;bin_p1_n0c = 0;bin_p2_n0c = 0;bin_p3_n0c = 0;bin_p4_n0c = 0;bin_p5_n0c = 0;bin_p6_n0c = 0;bin_p7_n0c = 0;bin_p8_n0c = 0;
bin_p0_n1c = 0;bin_p1_n1c = 0;bin_p2_n1c = 0;bin_p3_n1c = 0;bin_p4_n1c = 0;bin_p5_n1c = 0;bin_p6_n1c = 0;bin_p7_n1c = 0;bin_p8_n1c = 0;
bin_p0_n2c = 0;bin_p1_n2c = 0;bin_p2_n2c = 0;bin_p3_n2c = 0;bin_p4_n2c = 0;bin_p5_n2c = 0;bin_p6_n2c = 0;bin_p7_n2c = 0;bin_p8_n2c = 0;
bin_p0_n3c = 0;bin_p1_n3c = 0;bin_p2_n3c = 0;bin_p3_n3c = 0;bin_p4_n3c = 0;bin_p5_n3c = 0;bin_p6_n3c = 0;bin_p7_n3c = 0;bin_p8_n3c = 0;
bin_p0_n4c = 0;bin_p1_n4c = 0;bin_p2_n4c = 0;bin_p3_n4c = 0;bin_p4_n4c = 0;bin_p5_n4c = 0;bin_p6_n4c = 0;bin_p7_n4c = 0;bin_p8_n4c = 0;
bin_p0_n5c = 0;bin_p1_n5c = 0;bin_p2_n5c = 0;bin_p3_n5c = 0;bin_p4_n5c = 0;bin_p5_n5c = 0;bin_p6_n5c = 0;bin_p7_n5c = 0;bin_p8_n5c = 0;
bin_p0_n6c = 0;bin_p1_n6c = 0;bin_p2_n6c = 0;bin_p3_n6c = 0;bin_p4_n6c = 0;bin_p5_n6c = 0;bin_p6_n6c = 0;bin_p7_n6c = 0;bin_p8_n6c = 0;
bin_p0_n7c = 0;bin_p1_n7c = 0;bin_p2_n7c = 0;bin_p3_n7c = 0;bin_p4_n7c = 0;bin_p5_n7c = 0;bin_p6_n7c = 0;bin_p7_n7c = 0;bin_p8_n7c = 0;
bin_p0_n8c = 0;bin_p1_n8c = 0;bin_p2_n8c = 0;bin_p3_n8c = 0;bin_p4_n8c = 0;bin_p5_n8c = 0;bin_p6_n8c = 0;bin_p7_n8c = 0;bin_p8_n8c = 0;

for jj = 1:length(pirep_sev_dry_array)
  if pirep_sev_dry_array(jj) < 0.5
    if nirss_hazard_PIREP_dry_closest_array(jj) < 0.5
      bin_p0_n0c = bin_p0_n0c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 1.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 0.5
      bin_p0_n1c = bin_p0_n1c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 2.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 1.5
      bin_p0_n2c = bin_p0_n2c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 3.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 2.5
      bin_p0_n3c = bin_p0_n3c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 4.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 3.5
      bin_p0_n4c = bin_p0_n4c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 5.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 4.5
      bin_p0_n5c = bin_p0_n5c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 6.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 5.5
      bin_p0_n6c = bin_p0_n6c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 7.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 6.5
      bin_p0_n7c = bin_p0_n7c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) >= 7.5
      bin_p0_n8c = bin_p0_n8c + 1;
    end
  elseif pirep_sev_dry_array(jj) < 1.5 & pirep_sev_dry_array(jj) >= 0.5
    if nirss_hazard_PIREP_dry_closest_array(jj) < 0.5
      bin_p1_n0c = bin_p1_n0c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 1.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 0.5
      bin_p1_n1c = bin_p1_n1c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 2.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 1.5
      bin_p1_n2c = bin_p1_n2c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 3.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 2.5
      bin_p1_n3c = bin_p1_n3c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 4.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 3.5
      bin_p1_n4c = bin_p1_n4c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 5.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 4.5
      bin_p1_n5c = bin_p1_n5c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 6.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 5.5
      bin_p1_n6c = bin_p1_n6c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 7.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 6.5
      bin_p1_n7c = bin_p1_n7c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) >= 7.5
      bin_p1_n8c = bin_p1_n8c + 1;
    end
  elseif pirep_sev_dry_array(jj) < 2.5 & pirep_sev_dry_array(jj) >= 1.5
    if nirss_hazard_PIREP_dry_closest_array(jj) < 0.5
      bin_p2_n0c = bin_p2_n0c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 1.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 0.5
      bin_p2_n1c = bin_p2_n1c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 2.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 1.5
      bin_p2_n2c = bin_p2_n2c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 3.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 2.5
      bin_p2_n3c = bin_p2_n3c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 4.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 3.5
      bin_p2_n4c = bin_p2_n4c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 5.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 4.5
      bin_p2_n5c = bin_p2_n5c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 6.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 5.5
      bin_p2_n6c = bin_p2_n6c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 7.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 6.5
      bin_p2_n7c = bin_p2_n7c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) >= 7.5
      bin_p2_n8c = bin_p2_n8c + 1;
    end
  elseif pirep_sev_dry_array(jj) < 3.5 & pirep_sev_dry_array(jj) >= 2.5
    if nirss_hazard_PIREP_dry_closest_array(jj) < 0.5
      bin_p3_n0c                      = bin_p3_n0c + 1;
      counter_p3_n0c                  = bin_p3_n0c;
      bin_p3_n0c_date(counter_p3_n0c) = date_dry_array(jj);
      bin_p3_n0c_hh(counter_p3_n0c)   = pirep_hh_dry_array(jj);
      bin_p3_n0c_mm(counter_p3_n0c)   = pirep_mm_dry_array(jj);
      %bin_p3_n0c_datehhmm(counter_p3_n0c)   = [num2str(bin_p3_n0c_date(counter_p3_n0c)) ':' ...
      %                                         num2str(bin_p3_n0c_hh(counter_p3_n0c)) ...
      %                                         num2str(bin_p3_n0c_mm(counter_p3_n0c))];
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 1.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 0.5
      bin_p3_n1c = bin_p3_n1c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 2.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 1.5
      bin_p3_n2c = bin_p3_n2c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 3.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 2.5
      bin_p3_n3c = bin_p3_n3c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 4.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 3.5
      bin_p3_n4c = bin_p3_n4c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 5.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 4.5
      bin_p3_n5c = bin_p3_n5c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 6.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 5.5
      bin_p3_n6c = bin_p3_n6c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 7.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 6.5
      bin_p3_n7c = bin_p3_n7c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) >= 7.5
      bin_p3_n8c = bin_p3_n8c + 1;
    end
  elseif pirep_sev_dry_array(jj) < 4.5 & pirep_sev_dry_array(jj) >= 3.5
    if nirss_hazard_PIREP_dry_closest_array(jj) < 0.5
      bin_p4_n0c = bin_p4_n0c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 1.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 0.5
      bin_p4_n1c = bin_p4_n1c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 2.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 1.5
      bin_p4_n2c = bin_p4_n2c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 3.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 2.5
      bin_p4_n3c = bin_p4_n3c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 4.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 3.5
      bin_p4_n4c = bin_p4_n4c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 5.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 4.5
      bin_p4_n5c = bin_p4_n5c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 6.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 5.5
      bin_p4_n6c = bin_p4_n6c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 7.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 6.5
      bin_p4_n7c = bin_p4_n7c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) >= 7.5
      bin_p4_n8c = bin_p4_n8c + 1;
    end
  elseif pirep_sev_dry_array(jj) < 5.5 & pirep_sev_dry_array(jj) >= 4.5
    if nirss_hazard_PIREP_dry_closest_array(jj) < 0.5
      bin_p5_n0c                      = bin_p5_n0c + 1;
      counter_p5_n0c                  = bin_p5_n0c;
      bin_p5_n0c_date(counter_p5_n0c) = date_dry_array(jj);
      bin_p5_n0c_hh(counter_p5_n0c)   = pirep_hh_dry_array(jj);
      bin_p5_n0c_mm(counter_p5_n0c)   = pirep_mm_dry_array(jj);
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 1.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 0.5
      bin_p5_n1c = bin_p5_n1c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 2.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 1.5
      bin_p5_n2c = bin_p5_n2c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 3.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 2.5
      bin_p5_n3c = bin_p5_n3c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 4.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 3.5
      bin_p5_n4c = bin_p5_n4c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 5.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 4.5
      bin_p5_n5c = bin_p5_n5c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 6.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 5.5
      bin_p5_n6c = bin_p5_n6c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 7.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 6.5
      bin_p5_n7c = bin_p5_n7c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) >= 7.5
      bin_p5_n8c = bin_p5_n8c + 1;
    end
  elseif pirep_sev_dry_array(jj) < 6.5 & pirep_sev_dry_array(jj) >= 5.5
    if nirss_hazard_PIREP_dry_closest_array(jj) < 0.5
      bin_p6_n0c = bin_p6_n0c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 1.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 0.5
      bin_p6_n1c = bin_p6_n1c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 2.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 1.5
      bin_p6_n2c = bin_p6_n2c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 3.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 2.5
      bin_p6_n3c = bin_p6_n3c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 4.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 3.5
      bin_p6_n4c = bin_p6_n4c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 5.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 4.5
      bin_p6_n5c = bin_p6_n5c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 6.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 5.5
      bin_p6_n6c = bin_p6_n6c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 7.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 6.5
      bin_p6_n7c = bin_p6_n7c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) >= 7.5
      bin_p6_n8c = bin_p6_n8c + 1;
    end
  elseif pirep_sev_dry_array(jj) < 7.5 & pirep_sev_dry_array(jj) >= 6.5
    if nirss_hazard_PIREP_dry_closest_array(jj) < 0.5
      bin_p7_n0c = bin_p7_n0c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 1.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 0.5
      bin_p7_n1c = bin_p7_n1c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 2.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 1.5
      bin_p7_n2c = bin_p7_n2c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 3.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 2.5
      bin_p7_n3c = bin_p7_n3c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 4.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 3.5
      bin_p7_n4c = bin_p7_n4c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 5.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 4.5
      bin_p7_n5c = bin_p7_n5c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 6.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 5.5
      bin_p7_n6c = bin_p7_n6c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 7.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 6.5
      bin_p7_n7c = bin_p7_n7c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) >= 7.5
      bin_p7_n8c = bin_p7_n8c + 1;
    end
  elseif pirep_sev_dry_array(jj) < 8.5 & pirep_sev_dry_array(jj) >= 7.5
    if nirss_hazard_PIREP_dry_closest_array(jj) < 0.5
      bin_p8_n0c = bin_p8_n0c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 1.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 0.5
      bin_p8_n1c = bin_p8_n1c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 2.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 1.5
      bin_p8_n2c = bin_p8_n2c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 3.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 2.5
      bin_p8_n3c = bin_p8_n3c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 4.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 3.5
      bin_p8_n4c = bin_p8_n4c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 5.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 4.5
      bin_p8_n5c = bin_p8_n5c + 1;
      counter_p8_n5c                  = bin_p8_n5c;
      bin_p8_n5c_date(counter_p8_n5c) = date_dry_array(jj);
      bin_p8_n5c_hh(counter_p8_n5c)   = pirep_hh_dry_array(jj);
      bin_p8_n5c_mm(counter_p8_n5c)   = pirep_mm_dry_array(jj);
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 6.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 5.5
      bin_p8_n6c = bin_p8_n6c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 7.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 6.5
      bin_p8_n7c = bin_p8_n7c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) >= 7.5
      bin_p8_n8c = bin_p8_n8c + 1;
    end
  end
end

bins_p_vs_nclose_dry = [bin_p0_n8c bin_p1_n8c bin_p2_n8c bin_p3_n8c bin_p4_n8c bin_p5_n8c bin_p6_n8c bin_p7_n8c bin_p8_n8c; bin_p0_n7c bin_p1_n7c bin_p2_n7c bin_p3_n7c bin_p4_n7c bin_p5_n7c bin_p6_n7c bin_p7_n7c bin_p8_n7c; bin_p0_n6c bin_p1_n6c bin_p2_n6c bin_p3_n6c bin_p4_n6c bin_p5_n6c bin_p6_n6c bin_p7_n6c bin_p8_n6c; bin_p0_n5c bin_p1_n5c bin_p2_n5c bin_p3_n5c bin_p4_n5c bin_p5_n5c bin_p6_n5c bin_p7_n5c bin_p8_n5c; bin_p0_n4c bin_p1_n4c bin_p2_n4c bin_p3_n4c bin_p4_n4c bin_p5_n4c bin_p6_n4c bin_p7_n4c bin_p8_n4c; bin_p0_n3c bin_p1_n3c bin_p2_n3c bin_p3_n3c bin_p4_n3c bin_p5_n3c bin_p6_n3c bin_p7_n3c bin_p8_n3c; bin_p0_n2c bin_p1_n2c bin_p2_n2c bin_p3_n2c bin_p4_n2c bin_p5_n2c bin_p6_n2c bin_p7_n2c bin_p8_n2c; bin_p0_n1c bin_p1_n1c bin_p2_n1c bin_p3_n1c bin_p4_n1c bin_p5_n1c bin_p6_n1c bin_p7_n1c bin_p8_n1c; bin_p0_n0c bin_p1_n0c bin_p2_n0c bin_p3_n0c bin_p4_n0c bin_p5_n0c bin_p6_n0c bin_p7_n0c bin_p8_n0c];

%---------------------------------------------
% sum up bins for nirss dry closest versus cip
%---------------------------------------------
bin_c0_n0c = 0;bin_c1_n0c = 0;bin_c2_n0c = 0;bin_c3_n0c = 0;bin_c4_n0c = 0;bin_c5_n0c = 0;bin_c6_n0c = 0;bin_c7_n0c = 0;bin_c8_n0c = 0;
bin_c0_n1c = 0;bin_c1_n1c = 0;bin_c2_n1c = 0;bin_c3_n1c = 0;bin_c4_n1c = 0;bin_c5_n1c = 0;bin_c6_n1c = 0;bin_c7_n1c = 0;bin_c8_n1c = 0;
bin_c0_n2c = 0;bin_c1_n2c = 0;bin_c2_n2c = 0;bin_c3_n2c = 0;bin_c4_n2c = 0;bin_c5_n2c = 0;bin_c6_n2c = 0;bin_c7_n2c = 0;bin_c8_n2c = 0;
bin_c0_n3c = 0;bin_c1_n3c = 0;bin_c2_n3c = 0;bin_c3_n3c = 0;bin_c4_n3c = 0;bin_c5_n3c = 0;bin_c6_n3c = 0;bin_c7_n3c = 0;bin_c8_n3c = 0;
bin_c0_n4c = 0;bin_c1_n4c = 0;bin_c2_n4c = 0;bin_c3_n4c = 0;bin_c4_n4c = 0;bin_c5_n4c = 0;bin_c6_n4c = 0;bin_c7_n4c = 0;bin_c8_n4c = 0;
bin_c0_n5c = 0;bin_c1_n5c = 0;bin_c2_n5c = 0;bin_c3_n5c = 0;bin_c4_n5c = 0;bin_c5_n5c = 0;bin_c6_n5c = 0;bin_c7_n5c = 0;bin_c8_n5c = 0;
bin_c0_n6c = 0;bin_c1_n6c = 0;bin_c2_n6c = 0;bin_c3_n6c = 0;bin_c4_n6c = 0;bin_c5_n6c = 0;bin_c6_n6c = 0;bin_c7_n6c = 0;bin_c8_n6c = 0;
bin_c0_n7c = 0;bin_c1_n7c = 0;bin_c2_n7c = 0;bin_c3_n7c = 0;bin_c4_n7c = 0;bin_c5_n7c = 0;bin_c6_n7c = 0;bin_c7_n7c = 0;bin_c8_n7c = 0;
bin_c0_n8c = 0;bin_c1_n8c = 0;bin_c2_n8c = 0;bin_c3_n8c = 0;bin_c4_n8c = 0;bin_c5_n8c = 0;bin_c6_n8c = 0;bin_c7_n8c = 0;bin_c8_n8c = 0;

for jj = 1:length(cipice_sev_PIREP_dry_array)

  if cipice_sev_PIREP_dry_array(jj) < 0.02
    cipice_sev_PIREP_dry_array2(jj) = 0;
  elseif cipice_sev_PIREP_dry_array(jj) < 0.175 & cipice_sev_PIREP_dry_array(jj) >= 0.02
    cipice_sev_PIREP_dry_array2(jj) = 1;
  elseif cipice_sev_PIREP_dry_array(jj) < 0.275 & cipice_sev_PIREP_dry_array(jj) >= 0.175
    cipice_sev_PIREP_dry_array2(jj) = 2;
  elseif cipice_sev_PIREP_dry_array(jj) < 0.375 & cipice_sev_PIREP_dry_array(jj) >= 0.275
    cipice_sev_PIREP_dry_array2(jj) = 3;
  elseif cipice_sev_PIREP_dry_array(jj) < 0.538 & cipice_sev_PIREP_dry_array(jj) >= 0.375
    cipice_sev_PIREP_dry_array2(jj) = 4;
  elseif cipice_sev_PIREP_dry_array(jj) < 0.700 & cipice_sev_PIREP_dry_array(jj) >= 0.538
    cipice_sev_PIREP_dry_array2(jj) = 5;
  elseif cipice_sev_PIREP_dry_array(jj) < 0.800 & cipice_sev_PIREP_dry_array(jj) >= 0.700
    cipice_sev_PIREP_dry_array2(jj) = 6;
  elseif cipice_sev_PIREP_dry_array(jj) < 0.900 & cipice_sev_PIREP_dry_array(jj) >= 0.800
    cipice_sev_PIREP_dry_array2(jj) = 7;
  elseif cipice_sev_PIREP_dry_array(jj) <= 1.000 & cipice_sev_PIREP_dry_array(jj) >= 0.900
    cipice_sev_PIREP_dry_array2(jj) = 8;
  end

  if cipice_sev_PIREP_dry_array2(jj) < 0.5
    if nirss_hazard_PIREP_dry_closest_array(jj) < 0.5
      bin_c0_n0c = bin_c0_n0c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 1.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 0.5
      bin_c0_n1c = bin_c0_n1c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 2.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 1.5
      bin_c0_n2c = bin_c0_n2c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 3.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 2.5
      bin_c0_n3c = bin_c0_n3c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 4.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 3.5
      bin_c0_n4c = bin_c0_n4c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 5.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 4.5
      bin_c0_n5c = bin_c0_n5c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 6.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 5.5
      bin_c0_n6c = bin_c0_n6c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 7.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 6.5
      bin_c0_n7c = bin_c0_n7c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) >= 7.5
      bin_c0_n8c = bin_c0_n8c + 1;
    end
  elseif cipice_sev_PIREP_dry_array2(jj) < 1.5 & cipice_sev_PIREP_dry_array2(jj) >= 0.5
    if nirss_hazard_PIREP_dry_closest_array(jj) < 0.5
      bin_c1_n0c = bin_c1_n0c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 1.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 0.5
      bin_c1_n1c = bin_c1_n1c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 2.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 1.5
      bin_c1_n2c = bin_c1_n2c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 3.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 2.5
      bin_c1_n3c = bin_c1_n3c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 4.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 3.5
      bin_c1_n4c = bin_c1_n4c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 5.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 4.5
      bin_c1_n5c = bin_c1_n5c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 6.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 5.5
      bin_c1_n6c = bin_c1_n6c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 7.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 6.5
      bin_c1_n7c = bin_c1_n7c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) >= 7.5
      bin_c1_n8c = bin_c1_n8c + 1;
    end
  elseif cipice_sev_PIREP_dry_array2(jj) < 2.5 & cipice_sev_PIREP_dry_array2(jj) >= 1.5
    if nirss_hazard_PIREP_dry_closest_array(jj) < 0.5
      bin_c2_n0c = bin_c2_n0c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 1.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 0.5
      bin_c2_n1c = bin_c2_n1c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 2.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 1.5
      bin_c2_n2c = bin_c2_n2c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 3.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 2.5
      bin_c2_n3c = bin_c2_n3c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 4.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 3.5
      bin_c2_n4c = bin_c2_n4c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 5.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 4.5
      bin_c2_n5c = bin_c2_n5c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 6.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 5.5
      bin_c2_n6c = bin_c2_n6c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 7.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 6.5
      bin_c2_n7c = bin_c2_n7c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) >= 7.5
      bin_c2_n8c = bin_c2_n8c + 1;
    end
  elseif cipice_sev_PIREP_dry_array2(jj) < 3.5 & cipice_sev_PIREP_dry_array2(jj) >= 2.5
    if nirss_hazard_PIREP_dry_closest_array(jj) < 0.5
      %bin_c3_n0c                      = bin_c3_n0c + 1;
      %counter_c3_n0c                  = bin_c3_n0c;
      %bin_c3_n0c_date(counter_c3_n0c) = date_dry_array(jj);
      %bin_c3_n0c_hh(counter_c3_n0c)   = pirep_hh_dry_array(jj);
      %bin_c3_n0c_mm(counter_c3_n0c)   = pirep_mm_dry_array(jj);
      %bin_p3_n0c_datehhmm(counter_p3_n0c)   = [num2str(bin_p3_n0c_date(counter_p3_n0c)) ':' ...
      %                                         num2str(bin_p3_n0c_hh(counter_p3_n0c)) ...
      %                                         num2str(bin_p3_n0c_mm(counter_p3_n0c))];
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 1.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 0.5
      bin_c3_n1c = bin_c3_n1c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 2.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 1.5
      bin_c3_n2c = bin_c3_n2c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 3.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 2.5
      bin_c3_n3c = bin_c3_n3c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 4.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 3.5
      bin_c3_n4c = bin_c3_n4c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 5.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 4.5
      bin_c3_n5c = bin_c3_n5c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 6.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 5.5
      bin_c3_n6c = bin_c3_n6c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 7.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 6.5
      bin_c3_n7c = bin_c3_n7c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) >= 7.5
      bin_c3_n8c = bin_c3_n8c + 1;
    end
  elseif cipice_sev_PIREP_dry_array2(jj) < 4.5 & cipice_sev_PIREP_dry_array2(jj) >= 3.5
    if nirss_hazard_PIREP_dry_closest_array(jj) < 0.5
      bin_c4_n0c = bin_c4_n0c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 1.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 0.5
      bin_c4_n1c = bin_c4_n1c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 2.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 1.5
      bin_c4_n2c = bin_c4_n2c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 3.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 2.5
      bin_c4_n3c = bin_c4_n3c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 4.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 3.5
      bin_c4_n4c = bin_c4_n4c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 5.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 4.5
      bin_c4_n5c = bin_c4_n5c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 6.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 5.5
      bin_c4_n6c = bin_c4_n6c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 7.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 6.5
      bin_c4_n7c = bin_c4_n7c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) >= 7.5
      bin_c4_n8c = bin_c4_n8c + 1;
    end
  elseif cipice_sev_PIREP_dry_array2(jj) < 5.5 & cipice_sev_PIREP_dry_array2(jj) >= 4.5
    if nirss_hazard_PIREP_dry_closest_array(jj) < 0.5
      bin_c5_n0c                      = bin_c5_n0c + 1;
      counter_c5_n0c                  = bin_c5_n0c;
      bin_c5_n0c_date(counter_c5_n0c) = date_dry_array(jj);
      bin_c5_n0c_hh(counter_c5_n0c)   = pirep_hh_dry_array(jj);
      bin_c5_n0c_mm(counter_c5_n0c)   = pirep_mm_dry_array(jj);
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 1.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 0.5
      bin_c5_n1c = bin_c5_n1c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 2.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 1.5
      bin_c5_n2c = bin_c5_n2c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 3.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 2.5
      bin_c5_n3c = bin_c5_n3c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 4.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 3.5
      bin_c5_n4c = bin_c5_n4c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 5.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 4.5
      bin_c5_n5c = bin_c5_n5c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 6.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 5.5
      bin_c5_n6c = bin_c5_n6c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 7.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 6.5
      bin_c5_n7c = bin_c5_n7c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) >= 7.5
      bin_c5_n8c = bin_c5_n8c + 1;
    end
  elseif cipice_sev_PIREP_dry_array2(jj) < 6.5 & cipice_sev_PIREP_dry_array2(jj) >= 5.5
    if nirss_hazard_PIREP_dry_closest_array(jj) < 0.5
      bin_c6_n0c = bin_c6_n0c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 1.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 0.5
      bin_c6_n1c = bin_c6_n1c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 2.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 1.5
      bin_c6_n2c = bin_c6_n2c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 3.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 2.5
      bin_c6_n3c = bin_c6_n3c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 4.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 3.5
      bin_c6_n4c = bin_c6_n4c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 5.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 4.5
      bin_c6_n5c = bin_c6_n5c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 6.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 5.5
      bin_c6_n6c = bin_c6_n6c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 7.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 6.5
      bin_c6_n7c = bin_c6_n7c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) >= 7.5
      bin_c6_n8c = bin_c6_n8c + 1;
    end
  elseif cipice_sev_PIREP_dry_array2(jj) < 7.5 & cipice_sev_PIREP_dry_array2(jj) >= 6.5
    if nirss_hazard_PIREP_dry_closest_array(jj) < 0.5
      bin_c7_n0c = bin_c7_n0c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 1.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 0.5
      bin_c7_n1c = bin_c7_n1c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 2.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 1.5
      bin_c7_n2c = bin_c7_n2c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 3.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 2.5
      bin_c7_n3c = bin_c7_n3c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 4.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 3.5
      bin_c7_n4c = bin_c7_n4c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 5.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 4.5
      bin_c7_n5c = bin_c7_n5c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 6.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 5.5
      bin_c7_n6c = bin_c7_n6c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 7.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 6.5
      bin_c7_n7c = bin_c7_n7c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) >= 7.5
      bin_c7_n8c = bin_c7_n8c + 1;
    end
  elseif cipice_sev_PIREP_dry_array2(jj) < 8.5 & cipice_sev_PIREP_dry_array2(jj) >= 7.5
    if nirss_hazard_PIREP_dry_closest_array(jj) < 0.5
      bin_c8_n0c = bin_c8_n0c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 1.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 0.5
      bin_c8_n1c = bin_c8_n1c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 2.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 1.5
      bin_c8_n2c = bin_c8_n2c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 3.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 2.5
      bin_c8_n3c = bin_c8_n3c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 4.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 3.5
      bin_c8_n4c = bin_c8_n4c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 5.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 4.5
      bin_c8_n5c = bin_c8_n5c + 1;
      counter_c8_n5c                  = bin_c8_n5c;
      bin_c8_n5c_date(counter_c8_n5c) = date_dry_array(jj);
      bin_c8_n5c_hh(counter_c8_n5c)   = pirep_hh_dry_array(jj);
      bin_c8_n5c_mm(counter_c8_n5c)   = pirep_mm_dry_array(jj);
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 6.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 5.5
      bin_c8_n6c = bin_c8_n6c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) < 7.5 & nirss_hazard_PIREP_dry_closest_array(jj) >= 6.5
      bin_c8_n7c = bin_c8_n7c + 1;
    elseif nirss_hazard_PIREP_dry_closest_array(jj) >= 7.5
      bin_c8_n8c = bin_c8_n8c + 1;
    end
  end
end

bins_c_vs_nclose_dry = [bin_c0_n8c bin_c1_n8c bin_c2_n8c bin_c3_n8c bin_c4_n8c bin_c5_n8c bin_c6_n8c bin_c7_n8c bin_c8_n8c; bin_c0_n7c bin_c1_n7c bin_c2_n7c bin_c3_n7c bin_c4_n7c bin_c5_n7c bin_c6_n7c bin_c7_n7c bin_c8_n7c; bin_c0_n6c bin_c1_n6c bin_c2_n6c bin_c3_n6c bin_c4_n6c bin_c5_n6c bin_c6_n6c bin_c7_n6c bin_c8_n6c; bin_c0_n5c bin_c1_n5c bin_c2_n5c bin_c3_n5c bin_c4_n5c bin_c5_n5c bin_c6_n5c bin_c7_n5c bin_c8_n5c; bin_c0_n4c bin_c1_n4c bin_c2_n4c bin_c3_n4c bin_c4_n4c bin_c5_n4c bin_c6_n4c bin_c7_n4c bin_c8_n4c; bin_c0_n3c bin_c1_n3c bin_c2_n3c bin_c3_n3c bin_c4_n3c bin_c5_n3c bin_c6_n3c bin_c7_n3c bin_c8_n3c; bin_c0_n2c bin_c1_n2c bin_c2_n2c bin_c3_n2c bin_c4_n2c bin_c5_n2c bin_c6_n2c bin_c7_n2c bin_c8_n2c; bin_c0_n1c bin_c1_n1c bin_c2_n1c bin_c3_n1c bin_c4_n1c bin_c5_n1c bin_c6_n1c bin_c7_n1c bin_c8_n1c; bin_c0_n0c bin_c1_n0c bin_c2_n0c bin_c3_n0c bin_c4_n0c bin_c5_n0c bin_c6_n0c bin_c7_n0c bin_c8_n0c];

%---------------------------------------------
% sum up bins for CIP  versus pireps
%---------------------------------------------

bin_p0_c0 = 0;bin_p1_c0 = 0;bin_p2_c0 = 0;bin_p3_c0 = 0;bin_p4_c0 = 0;bin_p5_c0 = 0;bin_p6_c0 = 0;bin_p7_c0 = 0;bin_p8_c0 = 0;
bin_p0_c1 = 0;bin_p1_c1 = 0;bin_p2_c1 = 0;bin_p3_c1 = 0;bin_p4_c1 = 0;bin_p5_c1 = 0;bin_p6_c1 = 0;bin_p7_c1 = 0;bin_p8_c1 = 0;
bin_p0_c2 = 0;bin_p1_c2 = 0;bin_p2_c2 = 0;bin_p3_c2 = 0;bin_p4_c2 = 0;bin_p5_c2 = 0;bin_p6_c2 = 0;bin_p7_c2 = 0;bin_p8_c2 = 0;
bin_p0_c3 = 0;bin_p1_c3 = 0;bin_p2_c3 = 0;bin_p3_c3 = 0;bin_p4_c3 = 0;bin_p5_c3 = 0;bin_p6_c3 = 0;bin_p7_c3 = 0;bin_p8_c3 = 0;
bin_p0_c4 = 0;bin_p1_c4 = 0;bin_p2_c4 = 0;bin_p3_c4 = 0;bin_p4_c4 = 0;bin_p5_c4 = 0;bin_p6_c4 = 0;bin_p7_c4 = 0;bin_p8_c4 = 0;
bin_p0_c5 = 0;bin_p1_c5 = 0;bin_p2_c5 = 0;bin_p3_c5 = 0;bin_p4_c5 = 0;bin_p5_c5 = 0;bin_p6_c5 = 0;bin_p7_c5 = 0;bin_p8_c5 = 0;
bin_p0_c6 = 0;bin_p1_c6 = 0;bin_p2_c6 = 0;bin_p3_c6 = 0;bin_p4_c6 = 0;bin_p5_c6 = 0;bin_p6_c6 = 0;bin_p7_c6 = 0;bin_p8_c6 = 0;
bin_p0_c7 = 0;bin_p1_c7 = 0;bin_p2_c7 = 0;bin_p3_c7 = 0;bin_p4_c7 = 0;bin_p5_c7 = 0;bin_p6_c7 = 0;bin_p7_c7 = 0;bin_p8_c7 = 0;
bin_p0_c8 = 0;bin_p1_c8 = 0;bin_p2_c8 = 0;bin_p3_c8 = 0;bin_p4_c8 = 0;bin_p5_c8 = 0;bin_p6_c8 = 0;bin_p7_c8 = 0;bin_p8_c8 = 0;

for jj = 1:length(cipice_sev_PIREP_dry_array)

  if cipice_sev_PIREP_dry_array(jj) < 0.02
    cipice_sev_PIREP_dry_array(jj) = 0;
  elseif cipice_sev_PIREP_dry_array(jj) < 0.175 & cipice_sev_PIREP_dry_array(jj) >= 0.02
    cipice_sev_PIREP_dry_array(jj) = 1;
  elseif cipice_sev_PIREP_dry_array(jj) < 0.275 & cipice_sev_PIREP_dry_array(jj) >= 0.175
    cipice_sev_PIREP_dry_array(jj) = 2;
  elseif cipice_sev_PIREP_dry_array(jj) < 0.375 & cipice_sev_PIREP_dry_array(jj) >= 0.275
    cipice_sev_PIREP_dry_array(jj) = 3;
  elseif cipice_sev_PIREP_dry_array(jj) < 0.538 & cipice_sev_PIREP_dry_array(jj) >= 0.375
    cipice_sev_PIREP_dry_array(jj) = 4;
  elseif cipice_sev_PIREP_dry_array(jj) < 0.700 & cipice_sev_PIREP_dry_array(jj) >= 0.538
    cipice_sev_PIREP_dry_array(jj) = 5;
  elseif cipice_sev_PIREP_dry_array(jj) < 0.800 & cipice_sev_PIREP_dry_array(jj) >= 0.700
    cipice_sev_PIREP_dry_array(jj) = 6;
  elseif cipice_sev_PIREP_dry_array(jj) < 0.900 & cipice_sev_PIREP_dry_array(jj) >= 0.800
    cipice_sev_PIREP_dry_array(jj) = 7;
  elseif cipice_sev_PIREP_dry_array(jj) <= 1.000 & cipice_sev_PIREP_dry_array(jj) >= 0.900
    cipice_sev_PIREP_dry_array(jj) = 8;
  end

  if pirep_sev_dry_array(jj) < 0.5
    if cipice_sev_PIREP_dry_array(jj) < 0.5
      bin_p0_c0 = bin_p0_c0 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 1.5 & cipice_sev_PIREP_dry_array(jj) >= 0.5
      bin_p0_c1 = bin_p0_c1 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 2.5 & cipice_sev_PIREP_dry_array(jj) >= 1.5
      bin_p0_c2 = bin_p0_c2 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 3.5 & cipice_sev_PIREP_dry_array(jj) >= 2.5
      bin_p0_c3 = bin_p0_c3 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 4.5 & cipice_sev_PIREP_dry_array(jj) >= 3.5
      bin_p0_c4 = bin_p0_c4 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 5.5 & cipice_sev_PIREP_dry_array(jj) >= 4.5
      bin_p0_c5 = bin_p0_c5 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 6.5 & cipice_sev_PIREP_dry_array(jj) >= 5.5
      bin_p0_c6 = bin_p0_c6 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 7.5 & cipice_sev_PIREP_dry_array(jj) >= 6.5
      bin_p0_c7 = bin_p0_c7 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) >= 7.5
      bin_p0_c8 = bin_p0_c8 + 1;
    end
  elseif pirep_sev_dry_array(jj) < 1.5 & pirep_sev_dry_array(jj) >= 0.5
    if cipice_sev_PIREP_dry_array(jj) < 0.5
      bin_p1_c0 = bin_p1_c0 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 1.5 & cipice_sev_PIREP_dry_array(jj) >= 0.5
      bin_p1_c1 = bin_p1_c1 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 2.5 & cipice_sev_PIREP_dry_array(jj) >= 1.5
      bin_p1_c2 = bin_p1_c2 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 3.5 & cipice_sev_PIREP_dry_array(jj) >= 2.5
      bin_p1_c3 = bin_p1_c3 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 4.5 & cipice_sev_PIREP_dry_array(jj) >= 3.5
      bin_p1_c4 = bin_p1_c4 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 5.5 & cipice_sev_PIREP_dry_array(jj) >= 4.5
      bin_p1_c5 = bin_p1_c5 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 6.5 & cipice_sev_PIREP_dry_array(jj) >= 5.5
      bin_p1_c6 = bin_p1_c6 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 7.5 & cipice_sev_PIREP_dry_array(jj) >= 6.5
      bin_p1_c7 = bin_p1_c7 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) >= 7.5
      bin_p1_c8 = bin_p1_c8 + 1;
    end
  elseif pirep_sev_dry_array(jj) < 2.5 & pirep_sev_dry_array(jj) >= 1.5
    if cipice_sev_PIREP_dry_array(jj) < 0.5
      bin_p2_c0 = bin_p2_c0 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 1.5 & cipice_sev_PIREP_dry_array(jj) >= 0.5
      bin_p2_c1 = bin_p2_c1 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 2.5 & cipice_sev_PIREP_dry_array(jj) >= 1.5
      bin_p2_c2 = bin_p2_c2 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 3.5 & cipice_sev_PIREP_dry_array(jj) >= 2.5
      bin_p2_c3 = bin_p2_c3 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 4.5 & cipice_sev_PIREP_dry_array(jj) >= 3.5
      bin_p2_c4 = bin_p2_c4 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 5.5 & cipice_sev_PIREP_dry_array(jj) >= 4.5
      bin_p2_c5 = bin_p2_c5 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 6.5 & cipice_sev_PIREP_dry_array(jj) >= 5.5
      bin_p2_c6 = bin_p2_c6 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 7.5 & cipice_sev_PIREP_dry_array(jj) >= 6.5
      bin_p2_c7 = bin_p2_c7 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) >= 7.5
      bin_p2_c8 = bin_p2_c8 + 1;
    end
  elseif pirep_sev_dry_array(jj) < 3.5 & pirep_sev_dry_array(jj) >= 2.5
    if cipice_sev_PIREP_dry_array(jj) < 0.5
      bin_p3_c0                      = bin_p3_c0 + 1;
      counter_p3_c0                  = bin_p3_c0;
      bin_p3_c0_date(counter_p3_c0)  = date_dry_array(jj);
      bin_p3_c0_hh(counter_p3_c0)    = pirep_hh_dry_array(jj);
      bin_p3_c0_mm(counter_p3_c0)    = pirep_mm_dry_array(jj);
    elseif cipice_sev_PIREP_dry_array(jj) < 1.5 & cipice_sev_PIREP_dry_array(jj) >= 0.5
      bin_p3_c1 = bin_p3_c1 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 2.5 & cipice_sev_PIREP_dry_array(jj) >= 1.5
      bin_p3_c2 = bin_p3_c2 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 3.5 & cipice_sev_PIREP_dry_array(jj) >= 2.5
      bin_p3_c3 = bin_p3_c3 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 4.5 & cipice_sev_PIREP_dry_array(jj) >= 3.5
      bin_p3_c4 = bin_p3_c4 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 5.5 & cipice_sev_PIREP_dry_array(jj) >= 4.5
      bin_p3_c5 = bin_p3_c5 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 6.5 & cipice_sev_PIREP_dry_array(jj) >= 5.5
      bin_p3_c6 = bin_p3_c6 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 7.5 & cipice_sev_PIREP_dry_array(jj) >= 6.5
      bin_p3_c7 = bin_p3_c7 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) >= 7.5
      bin_p3_c8 = bin_p3_c8 + 1;
    end
  elseif pirep_sev_dry_array(jj) < 4.5 & pirep_sev_dry_array(jj) >= 3.5
    if cipice_sev_PIREP_dry_array(jj) < 0.5
      bin_p1_c0 = bin_p4_c0 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 1.5 & cipice_sev_PIREP_dry_array(jj) >= 0.5
      bin_p4_c1 = bin_p4_c1 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 2.5 & cipice_sev_PIREP_dry_array(jj) >= 1.5
      bin_p4_c2 = bin_p4_c2 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 3.5 & cipice_sev_PIREP_dry_array(jj) >= 2.5
      bin_p4_c3 = bin_p4_c3 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 4.5 & cipice_sev_PIREP_dry_array(jj) >= 3.5
      bin_p4_c4 = bin_p4_c4 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 5.5 & cipice_sev_PIREP_dry_array(jj) >= 4.5
      bin_p4_c5 = bin_p4_c5 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 6.5 & cipice_sev_PIREP_dry_array(jj) >= 5.5
      bin_p4_c6 = bin_p4_c6 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 7.5 & cipice_sev_PIREP_dry_array(jj) >= 6.5
      bin_p4_c7 = bin_p4_c7 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) >= 7.5
      bin_p4_c8 = bin_p4_c8 + 1;
    end
  elseif pirep_sev_dry_array(jj) < 5.5 & pirep_sev_dry_array(jj) >= 4.5
    if cipice_sev_PIREP_dry_array(jj) < 0.5
      bin_p5_c0                      = bin_p5_c0 + 1;
      counter_p5_c0                  = bin_p5_c0;
      bin_p5_c0_date(counter_p5_c0)  = date_dry_array(jj);
      bin_p5_c0_hh(counter_p5_c0)    = pirep_hh_dry_array(jj);
      bin_p5_c0_mm(counter_p5_c0)    = pirep_mm_dry_array(jj);
    elseif cipice_sev_PIREP_dry_array(jj) < 1.5 & cipice_sev_PIREP_dry_array(jj) >= 0.5
      bin_p5_c1 = bin_p5_c1 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 2.5 & cipice_sev_PIREP_dry_array(jj) >= 1.5
      bin_p5_c2 = bin_p5_c2 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 3.5 & cipice_sev_PIREP_dry_array(jj) >= 2.5
      bin_p5_c3 = bin_p5_c3 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 4.5 & cipice_sev_PIREP_dry_array(jj) >= 3.5
      bin_p5_c4 = bin_p5_c4 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 5.5 & cipice_sev_PIREP_dry_array(jj) >= 4.5
      bin_p5_c5 = bin_p5_c5 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 6.5 & cipice_sev_PIREP_dry_array(jj) >= 5.5
      bin_p5_c6 = bin_p5_c6 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 7.5 & cipice_sev_PIREP_dry_array(jj) >= 6.5
      bin_p5_c7 = bin_p5_c7 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) >= 7.5
      bin_p5_c8 = bin_p5_c8 + 1;
    end
  elseif pirep_sev_dry_array(jj) < 6.5 & pirep_sev_dry_array(jj) >= 5.5
    if cipice_sev_PIREP_dry_array(jj) < 0.5
      bin_p6_c0 = bin_p6_c0 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 1.5 & cipice_sev_PIREP_dry_array(jj) >= 0.5
      bin_p6_c1 = bin_p6_c1 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 2.5 & cipice_sev_PIREP_dry_array(jj) >= 1.5
      bin_p6_c2 = bin_p6_c2 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 3.5 & cipice_sev_PIREP_dry_array(jj) >= 2.5
      bin_p6_c3 = bin_p6_c3 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 4.5 & cipice_sev_PIREP_dry_array(jj) >= 3.5
      bin_p6_c4 = bin_p6_c4 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 5.5 & cipice_sev_PIREP_dry_array(jj) >= 4.5
      bin_p6_c5 = bin_p6_c5 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 6.5 & cipice_sev_PIREP_dry_array(jj) >= 5.5
      bin_p6_c6 = bin_p6_c6 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 7.5 & cipice_sev_PIREP_dry_array(jj) >= 6.5
      bin_p6_c7 = bin_p6_c7 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) >= 7.5
      bin_p6_c8 = bin_p6_c8 + 1;
    end
  elseif pirep_sev_dry_array(jj) < 7.5 & pirep_sev_dry_array(jj) >= 6.5
    if cipice_sev_PIREP_dry_array(jj) < 0.5
      bin_p7_n0 = bin_p7_n0 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 1.5 & cipice_sev_PIREP_dry_array(jj) >= 0.5
      bin_p7_n1 = bin_p7_n1 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 2.5 & cipice_sev_PIREP_dry_array(jj) >= 1.5
      bin_p7_n2 = bin_p7_n2 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 3.5 & cipice_sev_PIREP_dry_array(jj) >= 2.5
      bin_p7_n3 = bin_p7_n3 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 4.5 & cipice_sev_PIREP_dry_array(jj) >= 3.5
      bin_p7_n4 = bin_p7_n4 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 5.5 & cipice_sev_PIREP_dry_array(jj) >= 4.5
      bin_p7_n5 = bin_p7_n5 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 6.5 & cipice_sev_PIREP_dry_array(jj) >= 5.5
      bin_p7_n6 = bin_p7_n6 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 7.5 & cipice_sev_PIREP_dry_array(jj) >= 6.5
      bin_p7_n7 = bin_p7_n7 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) >= 7.5
      bin_p7_n8 = bin_p7_n8 + 1;
    end
  elseif pirep_sev_dry_array(jj) < 8.5 & pirep_sev_dry_array(jj) >= 7.5
    if cipice_sev_PIREP_dry_array(jj) < 0.5
      bin_p8_c0 = bin_p8_c0 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 1.5 & cipice_sev_PIREP_dry_array(jj) >= 0.5
      bin_p8_c1 = bin_p8_c1 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 2.5 & cipice_sev_PIREP_dry_array(jj) >= 1.5
      bin_p8_c2 = bin_p8_c2 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 3.5 & cipice_sev_PIREP_dry_array(jj) >= 2.5
      bin_p8_c3 = bin_p8_c3 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 4.5 & cipice_sev_PIREP_dry_array(jj) >= 3.5
      bin_p8_c4 = bin_p8_c4 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 5.5 & cipice_sev_PIREP_dry_array(jj) >= 4.5
      bin_p8_c5 = bin_p8_c5 + 1;
      counter_p8_c5                  = bin_p8_c5;
      bin_p8_c5_date(counter_p8_c5) = date_dry_array(jj);
      bin_p8_c5_hh(counter_p8_c5)   = pirep_hh_dry_array(jj);
      bin_p8_c5_mm(counter_p8_c5)   = pirep_mm_dry_array(jj);
    elseif cipice_sev_PIREP_dry_array(jj) < 6.5 & cipice_sev_PIREP_dry_array(jj) >= 5.5
      bin_p8_c6 = bin_p8_c6 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) < 7.5 & cipice_sev_PIREP_dry_array(jj) >= 6.5
      bin_p8_c7 = bin_p8_c7 + 1;
    elseif cipice_sev_PIREP_dry_array(jj) >= 7.5
      bin_p8_c8 = bin_p8_c8 + 1;
    end
  end
end

bins_p_vs_c = [bin_p0_c8 bin_p1_c8 bin_p2_c8 bin_p3_c8 bin_p4_c8 bin_p5_c8 bin_p6_c8 bin_p7_c8 bin_p8_c8; bin_p0_c7 bin_p1_c7 bin_p2_c7 bin_p3_c7 bin_p4_c7 bin_p5_c7 bin_p6_c7 bin_p7_c7 bin_p8_c7; bin_p0_c6 bin_p1_c6 bin_p2_c6 bin_p3_c6 bin_p4_c6 bin_p5_c6 bin_p6_c6 bin_p7_c6 bin_p8_c6; bin_p0_c5 bin_p1_c5 bin_p2_c5 bin_p3_c5 bin_p4_c5 bin_p5_c5 bin_p6_c5 bin_p7_c5 bin_p8_c5; bin_p0_c4 bin_p1_c4 bin_p2_c4 bin_p3_c4 bin_p4_c4 bin_p5_c4 bin_p6_c4 bin_p7_c4 bin_p8_c4; bin_p0_c3 bin_p1_c3 bin_p2_c3 bin_p3_c3 bin_p4_c3 bin_p5_c3 bin_p6_c3 bin_p7_c3 bin_p8_c3; bin_p0_c2 bin_p1_c2 bin_p2_c2 bin_p3_c2 bin_p4_c2 bin_p5_c2 bin_p6_c2 bin_p7_c2 bin_p8_c2; bin_p0_c1 bin_p1_c1 bin_p2_c1 bin_p3_c1 bin_p4_c1 bin_p5_c1 bin_p6_c1 bin_p7_c1 bin_p8_c1; bin_p0_c0 bin_p1_c0 bin_p2_c0 bin_p3_c0 bin_p4_c0 bin_p5_c0 bin_p6_c0 bin_p7_c0 bin_p8_c0];

%---------------------------------------------
% sum up bins for nirss closest versus wet pireps
%---------------------------------------------
bin_p0_n0c = 0;bin_p1_n0c = 0;bin_p2_n0c = 0;bin_p3_n0c = 0;bin_p4_n0c = 0;bin_p5_n0c = 0;bin_p6_n0c = 0;bin_p7_n0c = 0;bin_p8_n0c = 0;
bin_p0_n1c = 0;bin_p1_n1c = 0;bin_p2_n1c = 0;bin_p3_n1c = 0;bin_p4_n1c = 0;bin_p5_n1c = 0;bin_p6_n1c = 0;bin_p7_n1c = 0;bin_p8_n1c = 0;
bin_p0_n2c = 0;bin_p1_n2c = 0;bin_p2_n2c = 0;bin_p3_n2c = 0;bin_p4_n2c = 0;bin_p5_n2c = 0;bin_p6_n2c = 0;bin_p7_n2c = 0;bin_p8_n2c = 0;
bin_p0_n3c = 0;bin_p1_n3c = 0;bin_p2_n3c = 0;bin_p3_n3c = 0;bin_p4_n3c = 0;bin_p5_n3c = 0;bin_p6_n3c = 0;bin_p7_n3c = 0;bin_p8_n3c = 0;
bin_p0_n4c = 0;bin_p1_n4c = 0;bin_p2_n4c = 0;bin_p3_n4c = 0;bin_p4_n4c = 0;bin_p5_n4c = 0;bin_p6_n4c = 0;bin_p7_n4c = 0;bin_p8_n4c = 0;
bin_p0_n5c = 0;bin_p1_n5c = 0;bin_p2_n5c = 0;bin_p3_n5c = 0;bin_p4_n5c = 0;bin_p5_n5c = 0;bin_p6_n5c = 0;bin_p7_n5c = 0;bin_p8_n5c = 0;
bin_p0_n6c = 0;bin_p1_n6c = 0;bin_p2_n6c = 0;bin_p3_n6c = 0;bin_p4_n6c = 0;bin_p5_n6c = 0;bin_p6_n6c = 0;bin_p7_n6c = 0;bin_p8_n6c = 0;
bin_p0_n7c = 0;bin_p1_n7c = 0;bin_p2_n7c = 0;bin_p3_n7c = 0;bin_p4_n7c = 0;bin_p5_n7c = 0;bin_p6_n7c = 0;bin_p7_n7c = 0;bin_p8_n7c = 0;
bin_p0_n8c = 0;bin_p1_n8c = 0;bin_p2_n8c = 0;bin_p3_n8c = 0;bin_p4_n8c = 0;bin_p5_n8c = 0;bin_p6_n8c = 0;bin_p7_n8c = 0;bin_p8_n8c = 0;

for jj = 1:length(pirep_sev_wet_array)
  if pirep_sev_dry_array(jj) < 0.5
    if nirss_hazard_PIREP_wet_closest_array(jj) < 0.5
      bin_p0_n0c = bin_p0_n0c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 1.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 0.5
      bin_p0_n1c = bin_p0_n1c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 2.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 1.5
      bin_p0_n2c = bin_p0_n2c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 3.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 2.5
      bin_p0_n3c = bin_p0_n3c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 4.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 3.5
      bin_p0_n4c = bin_p0_n4c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 5.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 4.5
      bin_p0_n5c = bin_p0_n5c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 6.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 5.5
      bin_p0_n6c = bin_p0_n6c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 7.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 6.5
      bin_p0_n7c = bin_p0_n7c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) >= 7.5
      bin_p0_n8c = bin_p0_n8c + 1;
    end
  elseif pirep_sev_wet_array(jj) < 1.5 & pirep_sev_wet_array(jj) >= 0.5
    if nirss_hazard_PIREP_wet_closest_array(jj) < 0.5
      bin_p1_n0c = bin_p1_n0c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 1.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 0.5
      bin_p1_n1c = bin_p1_n1c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 2.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 1.5
      bin_p1_n2c = bin_p1_n2c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 3.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 2.5
      bin_p1_n3c = bin_p1_n3c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 4.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 3.5
      bin_p1_n4c = bin_p1_n4c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 5.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 4.5
      bin_p1_n5c = bin_p1_n5c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 6.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 5.5
      bin_p1_n6c = bin_p1_n6c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 7.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 6.5
      bin_p1_n7c = bin_p1_n7c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) >= 7.5
      bin_p1_n8c = bin_p1_n8c + 1;
    end
  elseif pirep_sev_wet_array(jj) < 2.5 & pirep_sev_wet_array(jj) >= 1.5
    if nirss_hazard_PIREP_wet_closest_array(jj) < 0.5
      bin_p2_n0c = bin_p2_n0c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 1.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 0.5
      bin_p2_n1c = bin_p2_n1c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 2.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 1.5
      bin_p2_n2c = bin_p2_n2c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 3.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 2.5
      bin_p2_n3c = bin_p2_n3c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 4.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 3.5
      bin_p2_n4c = bin_p2_n4c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 5.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 4.5
      bin_p2_n5c = bin_p2_n5c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 6.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 5.5
      bin_p2_n6c = bin_p2_n6c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 7.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 6.5
      bin_p2_n7c = bin_p2_n7c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) >= 7.5
      bin_p2_n8c = bin_p2_n8c + 1;
    end
  elseif pirep_sev_wet_array(jj) < 3.5 & pirep_sev_wet_array(jj) >= 2.5
    if nirss_hazard_PIREP_wet_closest_array(jj) < 0.5
      bin_p3_n0c                      = bin_p3_n0c + 1;
      counter_p3_n0c                  = bin_p3_n0c;
      bin_p3_n0c_date(counter_p3_n0c) = date_wet_array(jj);
      bin_p3_n0c_hh(counter_p3_n0c)   = pirep_hh_wet_array(jj);
      bin_p3_n0c_mm(counter_p3_n0c)   = pirep_mm_wet_array(jj);
      %bin_p3_n0c_datehhmm(counter_p3_n0c)   = [num2str(bin_p3_n0c_date(counter_p3_n0c)) ':' ...
      %                                         num2str(bin_p3_n0c_hh(counter_p3_n0c)) ...
      %                                         num2str(bin_p3_n0c_mm(counter_p3_n0c))];
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 1.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 0.5
      bin_p3_n1c = bin_p3_n1c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 2.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 1.5
      bin_p3_n2c = bin_p3_n2c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 3.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 2.5
      bin_p3_n3c = bin_p3_n3c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 4.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 3.5
      bin_p3_n4c = bin_p3_n4c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 5.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 4.5
      bin_p3_n5c = bin_p3_n5c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 6.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 5.5
      bin_p3_n6c = bin_p3_n6c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 7.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 6.5
      bin_p3_n7c = bin_p3_n7c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) >= 7.5
      bin_p3_n8c = bin_p3_n8c + 1;
    end
  elseif pirep_sev_wet_array(jj) < 4.5 & pirep_sev_wet_array(jj) >= 3.5
    if nirss_hazard_PIREP_wet_closest_array(jj) < 0.5
      bin_p4_n0c = bin_p4_n0c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 1.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 0.5
      bin_p4_n1c = bin_p4_n1c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 2.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 1.5
      bin_p4_n2c = bin_p4_n2c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 3.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 2.5
      bin_p4_n3c = bin_p4_n3c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 4.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 3.5
      bin_p4_n4c = bin_p4_n4c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 5.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 4.5
      bin_p4_n5c = bin_p4_n5c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 6.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 5.5
      bin_p4_n6c = bin_p4_n6c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 7.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 6.5
      bin_p4_n7c = bin_p4_n7c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) >= 7.5
      bin_p4_n8c = bin_p4_n8c + 1;
    end
  elseif pirep_sev_wet_array(jj) < 5.5 & pirep_sev_wet_array(jj) >= 4.5
    if nirss_hazard_PIREP_wet_closest_array(jj) < 0.5
      bin_p5_n0c                      = bin_p5_n0c + 1;
      counter_p5_n0c                  = bin_p5_n0c;
      bin_p5_n0c_date(counter_p5_n0c) = date_wet_array(jj);
      bin_p5_n0c_hh(counter_p5_n0c)   = pirep_hh_wet_array(jj);
      bin_p5_n0c_mm(counter_p5_n0c)   = pirep_mm_wet_array(jj);
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 1.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 0.5
      bin_p5_n1c = bin_p5_n1c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 2.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 1.5
      bin_p5_n2c = bin_p5_n2c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 3.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 2.5
      bin_p5_n3c = bin_p5_n3c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 4.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 3.5
      bin_p5_n4c = bin_p5_n4c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 5.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 4.5
      bin_p5_n5c = bin_p5_n5c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 6.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 5.5
      bin_p5_n6c = bin_p5_n6c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 7.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 6.5
      bin_p5_n7c = bin_p5_n7c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) >= 7.5
      bin_p5_n8c = bin_p5_n8c + 1;
    end
  elseif pirep_sev_wet_array(jj) < 6.5 & pirep_sev_wet_array(jj) >= 5.5
    if nirss_hazard_PIREP_wet_closest_array(jj) < 0.5
      bin_p6_n0c = bin_p6_n0c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 1.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 0.5
      bin_p6_n1c = bin_p6_n1c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 2.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 1.5
      bin_p6_n2c = bin_p6_n2c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 3.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 2.5
      bin_p6_n3c = bin_p6_n3c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 4.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 3.5
      bin_p6_n4c = bin_p6_n4c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 5.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 4.5
      bin_p6_n5c = bin_p6_n5c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 6.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 5.5
      bin_p6_n6c = bin_p6_n6c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 7.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 6.5
      bin_p6_n7c = bin_p6_n7c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) >= 7.5
      bin_p6_n8c = bin_p6_n8c + 1;
    end
  elseif pirep_sev_wet_array(jj) < 7.5 & pirep_sev_wet_array(jj) >= 6.5
    if nirss_hazard_PIREP_wet_closest_array(jj) < 0.5
      bin_p7_n0c = bin_p7_n0c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 1.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 0.5
      bin_p7_n1c = bin_p7_n1c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 2.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 1.5
      bin_p7_n2c = bin_p7_n2c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 3.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 2.5
      bin_p7_n3c = bin_p7_n3c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 4.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 3.5
      bin_p7_n4c = bin_p7_n4c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 5.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 4.5
      bin_p7_n5c = bin_p7_n5c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 6.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 5.5
      bin_p7_n6c = bin_p7_n6c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 7.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 6.5
      bin_p7_n7c = bin_p7_n7c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) >= 7.5
      bin_p7_n8c = bin_p7_n8c + 1;
    end
  elseif pirep_sev_wet_array(jj) < 8.5 & pirep_sev_wet_array(jj) >= 7.5
    if nirss_hazard_PIREP_wet_closest_array(jj) < 0.5
      bin_p8_n0c = bin_p8_n0c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 1.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 0.5
      bin_p8_n1c = bin_p8_n1c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 2.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 1.5
      bin_p8_n2c = bin_p8_n2c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 3.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 2.5
      bin_p8_n3c = bin_p8_n3c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 4.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 3.5
      bin_p8_n4c = bin_p8_n4c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 5.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 4.5
      bin_p8_n5c = bin_p8_n5c + 1;
      counter_p8_n5c                  = bin_p8_n5c;
      bin_p8_n5c_date(counter_p8_n5c) = date_wet_array(jj);
      bin_p8_n5c_hh(counter_p8_n5c)   = pirep_hh_wet_array(jj);
      bin_p8_n5c_mm(counter_p8_n5c)   = pirep_mm_wet_array(jj);
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 6.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 5.5
      bin_p8_n6c = bin_p8_n6c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) < 7.5 & nirss_hazard_PIREP_wet_closest_array(jj) >= 6.5
      bin_p8_n7c = bin_p8_n7c + 1;
    elseif nirss_hazard_PIREP_wet_closest_array(jj) >= 7.5
      bin_p8_n8c = bin_p8_n8c + 1;
    end
  end
end
 
bins_p_vs_nclose_wet = [bin_p0_n8c bin_p1_n8c bin_p2_n8c bin_p3_n8c bin_p4_n8c bin_p5_n8c bin_p6_n8c bin_p7_n8c bin_p8_n8c; bin_p0_n7c bin_p1_n7c bin_p2_n7c bin_p3_n7c bin_p4_n7c bin_p5_n7c bin_p6_n7c bin_p7_n7c bin_p8_n7c; bin_p0_n6c bin_p1_n6c bin_p2_n6c bin_p3_n6c bin_p4_n6c bin_p5_n6c bin_p6_n6c bin_p7_n6c bin_p8_n6c; bin_p0_n5c bin_p1_n5c bin_p2_n5c bin_p3_n5c bin_p4_n5c bin_p5_n5c bin_p6_n5c bin_p7_n5c bin_p8_n5c; bin_p0_n4c bin_p1_n4c bin_p2_n4c bin_p3_n4c bin_p4_n4c bin_p5_n4c bin_p6_n4c bin_p7_n4c bin_p8_n4c; bin_p0_n3c bin_p1_n3c bin_p2_n3c bin_p3_n3c bin_p4_n3c bin_p5_n3c bin_p6_n3c bin_p7_n3c bin_p8_n3c; bin_p0_n2c bin_p1_n2c bin_p2_n2c bin_p3_n2c bin_p4_n2c bin_p5_n2c bin_p6_n2c bin_p7_n2c bin_p8_n2c; bin_p0_n1c bin_p1_n1c bin_p2_n1c bin_p3_n1c bin_p4_n1c bin_p5_n1c bin_p6_n1c bin_p7_n1c bin_p8_n1c; bin_p0_n0c bin_p1_n0c bin_p2_n0c bin_p3_n0c bin_p4_n0c bin_p5_n0c bin_p6_n0c bin_p7_n0c bin_p8_n0c];

%---------------------------------
% plotting commands....
%---------------------------------
disp(['  start plotting commands: ']);

%figure(1);
%title(matchup_filename);

%---------------------------------------------------------
% output daily plots 
%----------------------------------------------------------
%cd_string = ['cd ' NIRSS_CIP_PIREP_output_path 'image'];
%eval(cd_string);
%saveas(1, [NIRSS_filename(1:8) '_NIRSS_CIP_PIREP.fig'], 'fig' );

% clear key variables for next nirss_date
%keyboard;
    
