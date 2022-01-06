%------------------------------------
% NAME:    crunch_NIRSSfusion_to_CIPPIREP
%
% PURPOSE: loop thru NIRSS fusion day files
%            read day file
%            read PIREP file
%              search for PIREPs on date
%            read CLE CIP gridpoint file
%            plot compare fields
%
% INPUTS: 
%
% OUTPUTS: 
%
% CREATED: 4.1.2010 dserke
%
% MODIFICATIONS:
%
% EXAMPLE: matlab -nodesktop
%          crunch_NIRSSfusion_to_CIPPIREP
%
% **********************************************
% ******  RUN IN 1 WEEK TIME CHUNKS ************* 
% **********************************************
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
output_data = 1;

%------------------------------------
% define looping dates 
%------------------------------------
start_date_yyyymmdd = num2str(20100322)
end_date_yyyymmdd   = num2str(20100330)
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
NIRSS_datadir_path          = '/d3/nirss/nasa_glen/fusion_output/';
CIP_datadir_path            = '/d3/nirss/CIP_output/';
PIREP_datadir_path          = '/d3/nirss/nasa_glen/PIREP_output/';
%PIREP_file                  = 'cle.40.nirss.dec';
PIREP_file                  = 'cle.111.nirss.dec';
NIRSS_CIP_PIREP_output_path = '/d3/nirss/nasa_glen/NIRSS_CIP_PIREP_output/';

PIREP_file_path             = [PIREP_datadir_path PIREP_file];

%------------------------------------
% get list of fusion data files  
%------------------------------------
cd_string = ['cd ' NIRSS_datadir_path];
eval(cd_string);
unix('ls 2* > nirss_files.txt');  
load('nirss_files.txt');

%keyboard;

%-----------------------------------------------------------
% loop through all desired dates (NIRSS, PIREP and CIP data
%-----------------------------------------------------------
day_num = 0;  % start a count of total days processed
for i = start_date:end_date

  day_num = day_num + 1;
  NIRSS_filename = datestr(i,'yyyymmddTHHMMSS');
  NIRSS_filedate = str2num(NIRSS_filename(1:8));
  disp('************************************');
  disp(['looking for file:     date= ' num2str(NIRSS_filedate)]);

  % check to see if fusion file exists
  nirss_file_ind = find(nirss_files == NIRSS_filedate); 

  % if fusion file does not exist, do no further processing for this day
  if size(nirss_file_ind,1) == 0
    disp('  error: no fusion file found!');


  % if fusion file does exist, do full processing for this day
  else

    %unzip filename if zipped 
    cd_string      = ['cd ' NIRSS_datadir_path];
    eval(cd_string);
    NIRSS_filename = [NIRSS_filename(1:8) '.nirss.gz'];
    unix(['gunzip ' NIRSS_filename]);
    NIRSS_filename = NIRSS_filename(1:end-3);

    % load fusion file
    disp(['  start loading file: date= ' num2str(NIRSS_filedate)]);
    fid            = fopen([NIRSS_datadir_path NIRSS_filename]);
    % find number of lines in nirss file
    file           = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '','Bufsize',10000);
    lines          = file{1};

    % rewind file to start
    fseek(fid,0,-1);

    % loop for reading each NIRSS line
    kk = 0; ll = 0; mm = 0; nn = 0;
    for jj = 1:length(lines)

      % get a single nirss data line into variable 'tline'
      tline          = fgetl(fid);
      data_field     = tline(1:5);

      % date_offset is char space diff for 1 vs 2 char mon and day
      if str2num(NIRSS_filename(5:6)) < 10
        mon_offset = 1;
      else
        mon_offset = 0;
      end
      if str2num(NIRSS_filename(7:8)) < 10
        day_offset = 1;
      else
        day_offset = 0;
      end
      date_offset    = mon_offset + day_offset;

      % parse data line based on type of field (RADAR,RADIO,TEST,NIRSS) 
      if (data_field == 'RADAR')
        ll          = ll + 1;
        gates       = tline(29-date_offset:31-date_offset);
        if gates == 'NOD'
          radar_refl{ll} = NaN*ones(1,480);
        elseif (gates == 'NOD' & radar_offset == 2)
          radar_refl{ll} = NaN*ones(1,500);
        end
        gates       = str2num(tline(29-date_offset:31-date_offset));
        if gates == 480
          tline2          = tline(43-date_offset:end);
          %radar_cloud{ll} = str2num(tline2(end-960:end));
          %tline2          = tline2(1:end-960);
          radar_reflect = textscan(tline2,'%9.7f','Delimiter',' ');
          radar_refl{ll}  = radar_reflect{1}(end-480+1:end);
        elseif gates == 500
          tline2          = tline(43-date_offset:end);
          radar_cloud{ll} = str2num(tline2(end-1000:end));
          tline2          = tline2(1:end-1000);
          radar_refl{ll}  = textscan(tline2,'%9.7f','Delimiter',' ');
        end

        %keyboard;

        %if ll == 693
        %  keyboard;
        %end

        clear tline{ll} tline2{ll};

      elseif (data_field == 'RADIO')
        nn            = nn + 1;
        rain_flag{nn} = str2num(tline(37-date_offset));
        %disp(rain_flag{nn});
        %keyboard;

      elseif (data_field == 'TEST1')
        mm = mm + 1;

        if length(tline) >= 27 
          if tline(27) == 'N'
            cip_haz{mm}    = NaN*ones(1,37);   
          else
            tline2         = tline(27-date_offset:end);
            cip_haz{mm}    = str2num(tline2);
          end 
        end

      elseif (data_field == 'NIRSS')
        kk              = kk + 1;
        hour{kk}        = str2num(tline(7:8));
        minu{kk}        = str2num(tline(10:11));
        sec{kk}         = str2num(tline(13:14));
        mon{kk}         = str2num(tline(16:17));
        day{kk}         = str2num(tline(19:20));
        year{kk}        = str2num(tline(22:25));
        gates_radar     = gates;
        gates           = tline(27-date_offset:29-date_offset);

        %if hour{kk} == 19 & minu{kk} == 35
        %  keyboard;
        %end

        if gates == 'NOD'
          nirss_haz{kk}    = NaN*ones(1,gates_radar);
          nirss_lwc{kk}{1} = NaN*ones(1,gates_radar);
        end      % end of if gates = 'NOD' 
        gates       = str2num(tline(27-date_offset:29-date_offset));
        %range_delta = str2num(tline(31:39));
        if gates == 480
          tline2          = tline(41-date_offset:end);
          nirss_haz{kk}   = str2num(tline2(end-960:end));
          tline3          = tline2(1:end-960);
          nirss_lwc{kk}   = textscan(tline3,'%11.9f','Delimiter',' ');
        elseif gates == 500
          tline2          = tline(41-date_offset:end);
          nirss_haz{kk}   = str2num(tline2(end-1000:end));
          tline3          = tline2(1:end-1000);
          nirss_lwc{kk}   = textscan(tline3,'%11.9f','Delimiter',' ');
        elseif gates == 420 
          tline2          = tline(39-date_offset:end);
          nirss_haz{kk}   = str2num(tline2(end-840:end));
          tline3          = tline2(1:end-840);
          nirss_lwc{kk}   = textscan(tline3,'%11.9f','Delimiter',' ');
        end      % end of if gates = 480

        %if kk == 693
        %  keyboard;
        %end

        clear tline{kk} tline2{kk};

      end        % end of if data_field

    end          % end of for jj=1:length(lines)

    disp(['  done loading file:  date= ' num2str(NIRSS_filedate)]);

    fclose(fid);

    %keyboard;

    for iii=1:length(rain_flag)
      if isempty(rain_flag{iii})
        rain_flag{iii} = 2;
      end
      %clear wet_flag;
      wet_flag(iii) = rain_flag{iii};
    end

    %keyboard;

    %rezip nirss_filename  
    unix(['gzip ' NIRSS_filename]);
    disp(['  done gzipping NIRSS file:']);

    %keyboard;

    %---------------------------------
    % calculate RES based on LWC and z
    %---------------------------------
    % only calc res if keyword set
    disp('  RES:');
    if compute_res == 1;
      disp(['    start calculating RES: ']);
      for i = 1:length(nirss_lwc)
        for j = 1:gates_radar
          try
            %radar_refl_factor(i,j) = radar_refl{i}{1}(j);
            radar_refl_factor(i,j) = 10^((radar_refl{i}(j))/10);
            %radar_refl_factor(i,j) = 10^((radar_refl{i}{1}(j))/10);
            %radar_refl_factor(i,j) = 10^((radar_refl{i}{1}(j)-(30-j*0.6))/10);
            rho_water              = 0.99812;                     % gm-3 @ -10 degrees C
            res(i,j)               = ((radar_refl_factor(i,j)*pi*rho_water)/(6*flipud(nirss_lwc{i}{1}(j))))^.33333;
            %units:                  [ ] [g
            nirss_hazard(i,j)      = nirss_haz{i}(j);
            time_hour(i,j)         = hour{i}(j);
          catch
          end      % end of try/catch
        end        % end of for j=1:length(nirss_lwc{1}{1})
      end          % end of for i=1:length(nirss_lwc)
      disp(['    done calculating RES: ']);

      % find maximum res in each column
      for pp = 1:length(res)
        try
          ind_res          = find(res(pp,1:50) ~= Inf & ~isnan(res(pp,1:50)));
          res_maxincol(pp) = max(res(pp,ind_res)); 
        catch
          res_maxincol(pp) = 0;
        end      % end of local try/catch
      end        % end of for pp=1:length(res)

    % elseif compute_res == 0, load hour and nirss_haz cells into simple arrays
    else
      disp(['    not calculating RES: ']);
      stop_loop = 0;
  
      for k = 1:length(nirss_lwc)
        time_hour(k)  = hour{k};

        % for each time in the daily nirss data
        for j = 1:length(nirss_lwc{k}{1})
          if length(nirss_haz{k}) == 0
            nirss_hazard(k,j)    = NaN;
          else
            nirss_hazard(k,j)    = nirss_haz{k}(j);
          end      % end of if length(nirss_haz{k}) == 0
        end        % end of for j=1:length(nirss_lwc{1}{1})
      end          % end of for k = 1:length(nirss_lwc)
    end            % end of if compute_res == 'yes'
%      % for each time in the daily nirss data
%      for i = 1:length(nirss_lwc)
%
%        % test if loop should be stopped
%        if stop_loop == 0
%          if length(nirss_lwc{i}) > 0
%            stop_loop = 1;
%  
%            for k = 1:length(nirss_lwc)
%              time_hour(k)  = hour{k};
%
%              % for each time in the daily nirss data
%              for j = 1:length(nirss_lwc{i}{1})
%                if length(nirss_haz{k}) == 0
%                  nirss_hazard(k,j)    = NaN;
%                else
%                  nirss_hazard(k,j)    = nirss_haz{k}(j);
%                end      % end of if length(nirss_haz{k}) == 0
%              end        % end of for j=1:length(nirss_lwc{1}{1})
%            end          % end of for k = 1:length(nirss_lwc)
%            %keyboard;
%            i = length(nirss_lwc) + 1;
%          end      % end of if length(nirss_lwc{i}) > 0
%        end        % end of if stop_loop == 0
%      end          % end of for i=1:length(nirss_lwc)
%    end            % end of if compute_res == 'yes'

    %keyboard;

    %---------------------------------
    % work with matched CIP data 
    %---------------------------------

    % try reformating cip_haz cell to cip_hazard matrix
    %try
    %  for i =1:length(cip_haz)
    %    for j = 1:length(cip_haz{1})
    %      try
    %        cip_hazard(i,j)        = cip_haz{i}(j);
    %      catch
    %      end
    %    end
    %  end
    %catch
    %  disp(['  error: failed to load CIP data from NIRSS fusion file ']);
    %end

    % try loading CIP gridpoint dump files for given day 
    count_no_cip_file = 0;  % add to counter if file does not exist
    count_no_cip_data = 0;  % add to counter if data lines in file do not exist
    for kk = 0:23
  
      % append a leading zero if hour_loop is a single digit month
      if kk < 10 hour_loop = ['0' num2str(kk) '00'];
      else hour_loop = [num2str(kk) '00'];
      end
  
      if exist([CIP_datadir_path NIRSS_filename(1:8) '/' NIRSS_filename(1:8) '_' hour_loop '_i221_j123.txt']) 
        fid3           = fopen([CIP_datadir_path NIRSS_filename(1:8) '/' NIRSS_filename(1:8) ...
                                '_' hour_loop '_i221_j123.txt']); 
  
        num_CIP_lines = 336;
        %num_CIP_lines = 371-1-32-14-38-38;
        count_good_tline3 = 0; count_bad_tline3 = 0;
        for mm = 1:num_CIP_lines
          tline3 = fgetl(fid3);
          if size(tline3,1) > 0 & size(tline3,2) > 2
            count_good_tline3 = count_good_tline3 + 1;
            CIP_data_field    = tline3(1:7);
  
            if (CIP_data_field == 'MODEL D')
              header         = fgetl(fid3);
              %disp(['  ' CIP_data_field]);
  
              % load model data into array: (ll,kk,day_num) ll=each model hgt, kk=each model time, day_num=each day
              for ll = 1:37
                tline3      = fgetl(fid3);
                model_data  = parse_string(tline3,'space');
                count = 1;
                for lll = 1:length(model_data)
                  test = isempty(model_data{lll});
                  if test ~=1
                    model_data2{count} = model_data{lll};
                    count              = count + 1;
                  end
                end
                model_k(ll,kk+1,day_num)     = str2num(model_data2{1}); 
                model_p(ll,kk+1,day_num)     = str2num(model_data2{2}); 
                model_hgt(ll,kk+1,day_num)   = str2num(model_data2{3}); 
                model_t(ll,kk+1,day_num)     = str2num(model_data2{4}); 
                model_rh(ll,kk+1,day_num)    = str2num(model_data2{5}); 
                model_vv(ll,kk+1,day_num)    = str2num(model_data2{6}); 
                model_theta(ll,kk+1,day_num) = str2num(model_data2{7}); 
                %keyboard;
              end
            elseif (CIP_data_field == 'SURFACE')
              
            elseif (CIP_data_field == 'CIP OUT' & tline3(12) == '3')
              header         = fgetl(fid3);
            elseif (CIP_data_field == 'CIP ICI')
              header         = fgetl(fid3);
              %disp(['  ' CIP_data_field]);
              % load model data into array: (ll,kk,day_num) ll=each model hgt, kk=each model time, day_num=each day
              for ll = 1:36
                tline3       = fgetl(fid3); 
                cipice_data  = parse_string(tline3,'space');
                count = 1;
                for lll = 1:length(cipice_data)
                  test = isempty(cipice_data{lll});
                  if test ~=1
                    cipice_data2{count} = cipice_data{lll};
                    count               = count + 1;
                  end
                end
                cipice_klevel(ll,kk+1)  = str2num(cipice_data2{1});
                cipice_ice(ll,kk+1)     = str2num(cipice_data2{2});
                cipice_sev(ll,kk+1)     = str2num(cipice_data2{3});
                cipice_sld(ll,kk+1)     = str2num(cipice_data2{4});
                cipice_prob(ll,kk+1)    = str2num(cipice_data2{5});
                %cipice_klevel(ll,kk+1,day_num)  = str2num(cipice_data2{1});
                %cipice_ice(ll,kk+1,day_num)     = str2num(cipice_data2{2});
                %cipice_sev(ll,kk+1,day_num)     = str2num(cipice_data2{3});
                %cipice_sld(ll,kk+1,day_num)     = str2num(cipice_data2{4});
                %cipice_prob(ll,kk+1,day_num)    = str2num(cipice_data2{5});
                %keyboard;
              end
            end      % end of if(CIP_data_field == ...) 
          else
            count_bad_tline3 = count_bad_tline3 + 1;
          end        % end of if(size(tline3,1) > 0
        end          % end of for mm=1:num_CIP_lines
        if count_bad_tline3 > 326
          disp('  error: CIP data file in CIP_datadir_path ... is not complete');
          count_no_cip_data = count_no_cip_data + 1;
        end

        fclose(fid3);  % close CIP file
      end              % end of if exist() 
    end                % end of for kk=0:23
    if count_no_cip_data >= 23
      disp('  error: No good CIP data on this day');
      day_num = day_num - 1; 
    end
    if count_no_cip_file >= 20
      disp('  error: No good CIP data on this day');
      day_num = day_num - 1; 
    end

    %keyboard;
  
    %---------------------------------
    % retrieve PIREP matchup data
    %---------------------------------
    % load PIREP data file
    fid2       = fopen(PIREP_file_path);
    nirss_date = str2num(NIRSS_filename(1:8));
    aa         = 0;
  
    % read in each PIREP data line
    disp(['  start loading PIREP file: ']);
    for pp = 1:14538  % length of pirep file in lines
    %for pp = 1:4941  % length of pirep file in lines
      tline2          = fgetl(fid2);
      pirep_data{pp}  = textscan(tline2,'%s%s%s%s%s%s%s','Delimiter',' ');
      [pirep_date,pirep_time,pirep_lat,pirep_lon,pirep_alti,pirep_sever,pirep_type] = deal(pirep_data{pp}{:});
  
      % check if there are PIREP data for this nirss date
      if (str2num(pirep_date{1}) == nirss_date)
        aa             = aa + 1;
        pirep_hhmm     = pirep_time{1};
        pirep_hh(aa)   = str2num(pirep_hhmm(1:2));
        pirep_mm(aa)   = str2num(pirep_hhmm(3:4));
        pirep_alt(aa)  = str2num(pirep_alti{1});
        if size(pirep_sever{1},1) == 0
          pirep_sev(aa)  = 0;  
        else
          pirep_sev(aa)  = str2num(pirep_sever{1});
        end
      end
  
    end
  
    % make a date field that is the length of the number of daily PIREPs
    date = str2num(NIRSS_filename(1:8))*ones(1,aa);

    % close PIREP data file
    fclose(fid2);
  
    % create array of nirss data altitudes where each gate is 98 feet apart
    nirss_altitude = ones(1,480);
    for lll = 1:480
      nirss_altitude(lll) = -98*(lll);
    end
  
    %keyboard;

    %---------------------------------------------------
    % calculate % of daily volume each product warns on 
    %---------------------------------------------------

    % for cip
    if exist('cipice_sev')
      for ii = 1:size(cipice_sev,2) 
        count_gt_0 = 0;
        for qq = 1:30
          if cipice_sev(qq,ii) > 0
            count_gt_0 = count_gt_0 + 1;
          end 
        end
        cip_warn_vol(ii) = count_gt_0/30;
        clear count_gt_0;
      end
      cip_warn_vol_day = sum(cip_warn_vol)/size(cipice_sev,2);
      cip_warn_vol_day
    end

    %keyboard;

    clear ii qq;

    % for nirss 
    for ii = 1:size(nirss_hazard,1)
      count_gt_0 = 0;
      for qq = 1:306
        if nirss_hazard(ii,qq) > 0
          count_gt_0 = count_gt_0 + 1;
        end 
      end
      nirss_warn_vol(ii) = count_gt_0/306;
      clear count_gt_0;
    end
    nirss_warn_vol_day = sum(nirss_warn_vol)/size(nirss_hazard,1);
    nirss_warn_vol_day

    %keyboard;


    %---------------------------------
    % plotting commands....
    %---------------------------------
    disp(['  start plotting commands: ']);
  
    %figure;
    %plot(res_maxincol')
    %mean_max_resincol = mean(res_maxincol(1:end));
    %%%mean(res_maxincol(1:210));
    %%mean(res_maxincol(211:810));
    %%mean(res_maxincol(811:end));
    %disp(['  mean max res in column= ' mean_max_resincol]);
    %disp(['  start plotting commands: ']);
  
    time_hours=time_hour;
    %time_hours=vertcat(time_hour,23);
  
    figure;
    %figure(day_num);
    title(NIRSS_filename);
  
    subplot(3,1,1);
    if compute_res == 1 
      %surf(time_hour,nirss_altitude/1000,res');
      %shading interp;
      %imagesc(flipud(res'),[0 30]);
      iind=find(res == inf);
      res(iind) = 0;
      imagesc(time_hours,nirss_altitude/1000,res',[0 60]);
      hold on;
      %keyboard;
      %imagesc(time_hours,nirss_altitude/1000,res',[0 100]);
      %imagesc(time_hours,nirss_altitude(1,1:100)/1000,res(:,1:100)');
      colorbar;
      title('RES [mu]');
      ylabel('Altitude [kft]');
      xlabel('Time [HH]');
    else
      disp('    not plotting RES');
    end
  
    subplot(3,1,2);
    %nirss_altitudes = repmat(nirss_altitude,1440,1);
    %time_hours      = repmat(time_hour,1,480);
    %surf(time_hour,nirss_altitude/1000,nirss_hazard');
    %shading interp;
    %imagesc(flipud(nirss_hazard'));
    imagesc(time_hours,nirss_altitude/1000,(nirss_hazard'),[0 8]);
    %imagesc(time_hours,(nirss_altitude(1:100))/1000,(nirss_hazard(:,1:100)'));
    colorbar;
    title('NIRSS HAZ Category [0-8]');
    ylabel('Altitude [kft]');
    xlabel('Time [HH]');
    hold on;
    hours = [];
    minus = [];
    for iii=1:length(hour)
      houri = hour{iii};
      minui = minu{iii};
      hours = [hours houri];
      minus = [minus minui];
    end

    %keyboard;

    % loop over all PIREPs for this day
    % find index of dry PIREPs (ind_dry_PIREP) for later use to carry alond only radiom=dry PIREPs
    ind_dry_PIREP = [];
    ind_wet_PIREP = [];
    for bb = 1: aa
      index_pirep_time      = find(hours == pirep_hh(bb) & minus == pirep_mm(bb) );
      %disp(bb)
      %disp(index_pirep_time);
      %disp(wet_flag(index_pirep_time));
      if wet_flag(index_pirep_time) == 0
        ind_dry_PIREP = [ind_dry_PIREP bb];
      else
        ind_wet_PIREP = [ind_wet_PIREP bb]; 
      end
      %if bb == 24
      %   keyboard;
      %end
    end

    disp('ind_dry_PIREP')
    disp(ind_dry_PIREP)
    disp('ind_wet_PIREP')
    disp(ind_wet_PIREP)

    %keyboard;

    for bb = 1: aa

      test      = find(hours == pirep_hh(bb) & minus == pirep_mm(bb) );

      if size(test,2) > 1 
        x_ind(bb) = test(1);
      else
        x_ind(bb) = find(hours == pirep_hh(bb) & minus == pirep_mm(bb) );
      end     % end of if size(test,2) > 1

      time(bb)     = pirep_hh(bb)+pirep_mm(bb)/60;
      y_ind(bb)    = -pirep_alt(bb)/1000;

      % plot numeric representation of PIREP on plot
      text(x_ind(bb),y_ind(bb),num2str(pirep_sev(bb)),'Color','r');
      %text(time(bb),y_ind(bb),num2str(pirep_sev(bb)),'Color','r');

      % test for nirss severity closest to PIREP altitude
      %  round off vertical height in feet of each PIREP 
      vert_ind(bb) = round(pirep_alt(bb)/98);
      count_alt    = 0;
      stop_alt     = 0;
      most_steps   = 480 - vert_ind(bb);
  
      %keyboard;
 
      % loop downward over each vertical increment until non-zero hazard is found 
      for lll=1:vert_ind(bb)

        % if not stop, do following processing
        if stop_alt == 0

          % if pirep sev is null, set closest nirss to current vert index and stop to 'yes'
          if pirep_sev(bb) == -1
            nirss_hazard_PIREP_closest(bb) = nirss_hazard(x_ind(bb),vert_ind(bb));
            pirep_sev(bb)                  = 0;
            stop_alt                       = 1;

          % if PIREP not null .....
          else

            % and we look from PIREP height to sfc with no non 0 nirss haz, set closest to 0
            if count_alt+1 == vert_ind(bb)
              nirss_hazard_PIREP_closest(bb) = 0;

            % look down AND up from PIREP height for closest NIRSS hazard
            % and if nirss haz at alt of PIREP minus alt counter is not 0, this PIREPs nirss 
            % haz is set to the closest non 0 nirss haz and stop set to 'yes'
            elseif nirss_hazard(x_ind(bb),vert_ind(bb)-count_alt) ~= 0
              nirss_hazard_PIREP_closest(bb) = nirss_hazard(x_ind(bb),vert_ind(bb)-count_alt);
              stop_alt                       = 1;
            elseif nirss_hazard(x_ind(bb),vert_ind(bb)+count_alt) ~= 0 & count_alt < most_steps
              nirss_hazard_PIREP_closest(bb) = nirss_hazard(x_ind(bb),vert_ind(bb)+count_alt);
              stop_alt                       = 1;
            end
            %if bb == 4
            %  keyboard
            %  nirss_hazard_PIREP_closest(bb)
            %  count_alt
            %  vert_ind(bb)
            %end

            count_alt = count_alt+1;

          end
        end
      end

      %keyboard;

      clear count_alt, stop_alt;

      % save other fields which are indexed to PIREP times
      if compute_res == 1;
        res_PIREP_array              = res(x_ind(bb),:);
        ind_res_PIREP_array          = find(res_PIREP_array ~= Inf);
        res_PIREP_max(bb)            = max(res_PIREP_array(ind_res_PIREP_array));
      end

      nirss_hazard_PIREP_array       = nirss_hazard(x_ind(bb),:);
      ind_nirss_hazard_PIREP_array   = find(nirss_hazard_PIREP_array ~= 0);
      nirss_hazard_PIREP_mean(bb)    = mean(nirss_hazard_PIREP_array(ind_nirss_hazard_PIREP_array));

      if size(ind_nirss_hazard_PIREP_array,2) == 0;
       nirss_hazard_PIREP_max(bb)    = NaN;
      else
        nirss_hazard_PIREP_max(bb)   = max(nirss_hazard_PIREP_array(ind_nirss_hazard_PIREP_array));
      end    % end of if size(ind_nirss_hazard_PIREP_array,2) == 0

      % now find all model and CIP fields indexed to PIREP times and altitudes
      % first, round PIREP hour based on PIREP minutes
      if pirep_mm(bb) >= 35
        ind_pirep_hour = pirep_hh(bb) + 1;
      else
        ind_pirep_hour = pirep_hh(bb);
      end
      
      % since CIP hours go from 0 to 23, we add 1 hour to 'ind_pirep_hour' field to access correct array position 
      ind_pirep_hour = ind_pirep_hour + 1;

      % secondly, find closest model height to PIREP alt
      hgt_diff_array     = abs(pirep_alt(bb)*.3048 - model_hgt(:,ind_pirep_hour));
      %hgt_diff_array     = abs(pirep_alt(bb)*.3048 - model_hgt(:,ind_pirep_hour,day_num));
      hgt_diff_array_min = hgt_diff_array(1);
      for tt=1:length(hgt_diff_array)
        if hgt_diff_array(tt) <= hgt_diff_array_min
          hgt_diff_array_min = hgt_diff_array(tt);
          ind_hgt(bb)        = tt;
        end 
      end

      % lastly, make new field of cipice severity that corresponds to PIREP time and alt
      if ind_hgt(bb) < 37
        cipice_sev_PIREP(bb) = cipice_sev(ind_hgt(bb),ind_pirep_hour);
      end
      %keyboard;
      clear ind_pirep_hour;

    end      % end of for bb=1:aa
  
    subplot(3,1,3);
    if count_no_cip_data >= 23 | count_no_cip_file >= 20
      disp('    not plotting cipice_sev ... no data');
    else
      for ss=1:size(cipice_sev,2)
        for rr=1:size(cipice_sev,1)
          if cipice_sev(rr,ss) < 0.02
            cipice_sev_P(rr,ss) = 0;
          elseif cipice_sev(rr,ss) < 0.175 & cipice_sev(rr,ss) >= 0.02
            cipice_sev_P(rr,ss) = 1;
          elseif cipice_sev(rr,ss) < 0.275 & cipice_sev(rr,ss) >= 0.175
            cipice_sev_P(rr,ss) = 2;
          elseif cipice_sev(rr,ss) < 0.375 & cipice_sev(rr,ss) >= 0.275
            cipice_sev_P(rr,ss) = 3;
          elseif cipice_sev(rr,ss) < 0.538 & cipice_sev(rr,ss) >= 0.375
            cipice_sev_P(rr,ss) = 4;
          elseif cipice_sev(rr,ss) < 0.700 & cipice_sev(rr,ss) >= 0.538
            cipice_sev_P(rr,ss) = 5;
          elseif cipice_sev(rr,ss) < 0.800 & cipice_sev(rr,ss) >= 0.700
            cipice_sev_P(rr,ss) = 6;
          elseif cipice_sev(rr,ss) < 0.900 & cipice_sev(rr,ss) >= 0.800
            cipice_sev_P(rr,ss) = 7;
          elseif cipice_sev(rr,ss) <= 1.000 & cipice_sev(rr,ss) >= 0.900
            cipice_sev_P(rr,ss) = 8;
          end
        end
      end

      imagesc(flipud(cipice_sev_P(:,:)),[0 8]);
      %imagesc(flipud(cipice_sev(:,:)),[0 1]);
      %imagesc(flipud(cipice_sev(:,:,day_num)),[0 1]);
    end
    %%imagesc(time_hour,cip_height(1,1:14)/1000,flipud(cip_hazard(:,1:14)'/14.2857),[0 7]);
    colorbar;
    title('CIP HAZ Category [0-1]');
  
    %keyboard;

    %---------------------------------------------------------
    % output certain data fields for later multi-year analysis
    %----------------------------------------------------------
    % test output_data keyword
    if output_data == 1 
  
      %keyboard;

      % execute change dir linux commands
      cd_string = ['cd ' NIRSS_CIP_PIREP_output_path ];
      %cd_string = ['cd ' NIRSS_CIP_PIREP_output_path 'data'];
      eval(cd_string);
      % execute rm existing daily mat file linux commands
      %rm_matfile1_string = ['rm ' NIRSS_filename(1:8) '.mat'];
      %rm_matfile2_string = ['rm ' NIRSS_filename(1:8) '_sevmatchup.mat'];
      %if exist([NIRSS_filename(1:8) '.mat']) == 2
      %  unix(rm_matfile1_string);
      %end
      %if exist([NIRSS_filename(1:8) '_sevmatchup.mat']) == 2
      %  unix(rm_matfile2_string);
      %end

      % setup fields to output based on whether NIRSS radiom was wet or dry at PIREP time 
      if (size(ind_dry_PIREP,1)) > 0 
        pirep_sev_dry                  = pirep_sev(ind_dry_PIREP);
        nirss_hazard_PIREP_dry_closest = nirss_hazard_PIREP_closest(ind_dry_PIREP);
        nirss_hazard_PIREP_dry_mean    = nirss_hazard_PIREP_mean(ind_dry_PIREP);
        nirss_hazard_PIREP_dry_max     = nirss_hazard_PIREP_max(ind_dry_PIREP);
        cipice_sev_PIREP_dry           = cipice_sev_PIREP(ind_dry_PIREP);
        pirep_hh_dry                   = pirep_hh(ind_dry_PIREP);
        pirep_mm_dry                   = pirep_mm(ind_dry_PIREP);
        pirep_alt_dry                  = pirep_alt(ind_dry_PIREP);
        date_dry                       = date(ind_dry_PIREP);
      elseif (isempty(ind_dry_PIREP))
        pirep_sev_dry                  = [];
        nirss_hazard_PIREP_dry_closest = [];
        nirss_hazard_PIREP_dry_mean    = [];
        nirss_hazard_PIREP_dry_max     = [];
        cipice_sev_PIREP_dry           = [];
        pirep_hh_dry                   = [];
        pirep_mm_dry                   = [];
        pirep_alt_dry                  = [];
        date_dry                       = [];
      end

      %keyboard;

      if (size(ind_wet_PIREP,1)) > 0 
        pirep_sev_wet                  = pirep_sev(ind_wet_PIREP);
        nirss_hazard_PIREP_wet_closest = nirss_hazard_PIREP_closest(ind_wet_PIREP);
        nirss_hazard_PIREP_wet_mean    = nirss_hazard_PIREP_mean(ind_wet_PIREP);
        nirss_hazard_PIREP_wet_max     = nirss_hazard_PIREP_max(ind_wet_PIREP);
        cipice_sev_PIREP_wet           = cipice_sev_PIREP(ind_wet_PIREP);
        pirep_hh_wet                   = pirep_hh(ind_wet_PIREP);
        pirep_mm_wet                   = pirep_mm(ind_wet_PIREP);
        pirep_alt_wet                  = pirep_alt(ind_wet_PIREP);
        date_wet                       = date(ind_wet_PIREP);
      elseif (isempty(ind_wet_PIREP))
        pirep_sev_wet                  = [];
        nirss_hazard_PIREP_wet_closest = [];
        nirss_hazard_PIREP_wet_mean    = [];
        nirss_hazard_PIREP_wet_max     = [];
        cipice_sev_PIREP_wet           = [];
        pirep_hh_wet                   = [];
        pirep_mm_wet                   = [];
        pirep_alt_wet                  = [];
        date_wet                       = [];
      end

      %keyboard

      % save files with daily data fields of interest
      disp('  outputting mat files to data dir');

      if compute_res == 1;
        if count_no_cip_data >= 23 | count_no_cip_file >= 20
          if aa > 0
            %len = length(pirep_sev);
            %nirss_hazard_PIREP_closest = nirss_hazard_PIREP_closest(1:len);
            %nirss_hazard_PIREP_mean    = nirss_hazard_PIREP_mean(1:len);
            %nirss_hazard_PIREP_max     = nirss_hazard_PIREP_max(1:len);
            %%pirep_hh                   = pirep_hh(1:len);
            %%pirep_mm                   = pirep_mm(1:len);
            %pirep_alt                  = pirep_alt(1:len);
            save([NIRSS_filename(1:8) '_sevmatchup.11.mat'], '-mat', 'pirep_sev_dry','pirep_sev_wet','res_PIREP_max','nirss_hazard_PIREP_dry_closest','nirss_hazard_PIREP_wet_closest','nirss_hazard_PIREP_mean','nirss_hazard_PIREP_max','pirep_hh_dry','pirep_mm_dry','pirep_alt_dry','date_dry','pirep_hh_wet','pirep_mm_wet','pirep_alt_wet','date_wet','nirss_warn_vol_day' );
            %clear pirep_sev_dry,pirep_sev_wet,res_PIREP_max,nirss_hazard_PIREP_dry_closest, ...
            %      nirss_hazard_PIREP_wet_closest,nirss_hazard_PIREP_mean,...
            %      nirss_hazard_PIREP_max,pirep_hh_dry,pirep_mm_dry,pirep_alt_dry,date_dry, ...
            %      pirep_hh_wet,pirep_mm_wet,pirep_alt_wet,date_wet;
          else
            disp('    not outputting files ... no matching PIREPs');
          end
        else
          if aa > 0
            %len = length(pirep_sev);
            %nirss_hazard_PIREP_closest = nirss_hazard_PIREP_closest(1:len);
            %nirss_hazard_PIREP_mean    = nirss_hazard_PIREP_mean(1:len);
            %nirss_hazard_PIREP_max     = nirss_hazard_PIREP_max(1:len);
            %pirep_hh                   = pirep_hh(1:len);
            %pirep_mm                   = pirep_mm(1:len);
            %pirep_alt                  = pirep_alt(1:len);
            save([NIRSS_filename(1:8) '_sevmatchup.11.mat'], '-mat', 'pirep_sev_dry','pirep_sev_wet','res_PIREP_max','nirss_hazard_PIREP_dry_closest','nirss_hazard_PIREP_wet_closest','nirss_hazard_PIREP_mean','nirss_hazard_PIREP_max','cipice_sev_PIREP_dry', 'pirep_hh','pirep_mm','pirep_alt','date','cipice_sev_PIREP_wet', 'pirep_hh_wet','pirep_mm_wet','pirep_alt_wet','date_wet','cip_warn_vol_day','nirss_warn_vol_day' );
            %clear pirep_sev_dry,pirep_sev_wet,res_PIREP_max,nirss_hazard_PIREP_dry_closest, ...
            %      nirss_hazard_PIREP_wet_closest,nirss_hazard_PIREP_mean,...
            %      nirss_hazard_PIREP_max,pirep_hh_dry,pirep_mm_dry,pirep_alt_dry,date_dry, ...
            %      pirep_hh_wet,pirep_mm_wet,pirep_alt_wet,date_wet;
          else
            disp('    not outputting files ... no matching PIREPs');
          end
        end
      else
        if count_no_cip_data >= 23 | count_no_cip_file >= 20
          if aa > 0
            %len = length(pirep_sev);
            %nirss_hazard_PIREP_closest = nirss_hazard_PIREP_closest(1:len);
            %nirss_hazard_PIREP_mean    = nirss_hazard_PIREP_mean(1:len);
            %nirss_hazard_PIREP_max     = nirss_hazard_PIREP_max(1:len);
            %pirep_hh                   = pirep_hh(1:len);
            %pirep_mm                   = pirep_mm(1:len);
            %pirep_alt                  = pirep_alt(1:len);
            save([NIRSS_filename(1:8) '_sevmatchup.11.mat'], '-mat', 'pirep_sev_dry','pirep_sev_wet','nirss_hazard_PIREP_dry_closest','nirss_hazard_PIREP_wet_closest','nirss_hazard_PIREP_mean','nirss_hazard_PIREP_max','pirep_hh_dry','pirep_mm_dry','pirep_alt_dry','date_dry','pirep_hh_wet','pirep_mm_wet','pirep_alt_wet','date_wet','nirss_warn_vol_day' );
            %clear pirep_sev_dry,pirep_sev_wet,nirss_hazard_PIREP_dry_closest, ...
            %      nirss_hazard_PIREP_wet_closest,nirss_hazard_PIREP_mean,...
            %      nirss_hazard_PIREP_max,pirep_hh_dry,pirep_mm_dry,pirep_alt_dry,date_dry, ...
            %      pirep_hh_wet,pirep_mm_wet,pirep_alt_wet,date_wet;
          else
            disp('    not outputting files ... no matching PIREPs');
          end
        else
          if aa > 0
            %len = length(pirep_sev);
            %nirss_hazard_PIREP_closest = nirss_hazard_PIREP_closest(1:len);
            %nirss_hazard_PIREP_mean    = nirss_hazard_PIREP_mean(1:len);
            %nirss_hazard_PIREP_max     = nirss_hazard_PIREP_max(1:len);
            %pirep_hh                   = pirep_hh(1:len);
            %pirep_mm                   = pirep_mm(1:len);
            %pirep_alt                  = pirep_alt(1:len);
            save([NIRSS_filename(1:8) '_sevmatchup.11.mat'], '-mat', 'pirep_sev_dry','pirep_sev_wet','nirss_hazard_PIREP_dry_closest','nirss_hazard_PIREP_wet_closest','nirss_hazard_PIREP_mean','nirss_hazard_PIREP_max','cipice_sev_PIREP_dry', 'pirep_hh_dry','pirep_mm_dry','pirep_alt_dry','date_dry','cipice_sev_PIREP_wet', 'pirep_hh_wet','pirep_mm_wet','pirep_alt_wet','date_wet','cip_warn_vol_day','nirss_warn_vol_day' );
            %clear pirep_sev_dry,pirep_sev_wet,nirss_hazard_PIREP_dry_closest, ...
            %      nirss_hazard_PIREP_wet_closest,nirss_hazard_PIREP_mean,...
            %      nirss_hazard_PIREP_max,pirep_hh_dry,pirep_mm_dry,pirep_alt_dry,date_dry, ...
            %      pirep_hh_wet,pirep_mm_wet,pirep_alt_wet,date_wet;
          else
            disp('    not outputting files ... no matching PIREPs');
          end
        end
      end
    elseif output_data == 0
      disp('  not outputting data file:');
    end
  
    clear ind_wet_PIREP ind_dry_PIREP; 

    %keyboard;
  
    %---------------------------------------------------------
    % output daily plots 
    %----------------------------------------------------------
    cd_string = ['cd ' NIRSS_CIP_PIREP_output_path 'image'];
    eval(cd_string);
    saveas(1, [NIRSS_filename(1:8) '_NIRSS_CIP_PIREP.fig'], 'fig' );
  
  %   % clear key variables for next nirss_date
  %   %keyboard;
  %   clear time_hours,nirss_hazard,model_hgt,cipice_sev, pirep_sev, nirss_hazard_PIREP_mean, nirss_hazard_PIREP_max, pirep_hh, pirep_mm, pirep_alt;
  %  if compute_res == 1;
  %    clear res, res_PIREP_max;
  %  end

    %keyboard;

    %close all;

  end    % end of if(size(nirss_file_ind))

end      % end of for i=start_date:end_date
