%---------------------------------------
% Name: compare_S_K_ILW_to_PIREP.m
%
% Purpose: find and plot nearest S and K radar profiles
%          as well as radiometer ILW values
%          to PIREPs that are within
%          60 km of the remote sensing values
%          with uncorrelated PIREPs (not within 30 min of each other)
%
% Usage:
%
% Created: 9.1.2011 dserke
%         
%
%---------------------------------------

%---------------------------------------
% open PIREP file
%---------------------------------------
fid   = fopen('/d1/serke/projects/NIRSS_NASA/code/PIREP_retrieve/all_NIRSS@Platte_60km.txt');
%fid   = fopen('/d1/serke/projects/NIRSS_NASA/code/PIREP_retrieve/all_NIRSS@KCLE_25km.txt');
%fid   = fopen('/d1/serke/projects/NIRSS_NASA/code/PIREP_retrieve/all_NIRSS@Platte_25km.txt');

%---------------------------------------
% load PIREP fields
%---------------------------------------
count = 1;

while 1

   tline            = fgetl(fid);
   keyboard;
   %disp(tline)
   if ~ischar(tline), break, end
   unix_time(count) = str2num(tline(1:10));
   lat(count)       = str2num(tline(12:17));
   lon(count)       = str2num(tline(19:25));
   temp(count)      = str2num(tline(59:61));
   alt_lo(count)    = str2num(tline(71:73));
   alt_hi(count)    = str2num(tline(75:77));
   severity(count)  = str2num(tline(80));

   count            = count + 1;

end

fclose(fid);

%---------------------------------------
% find PIREPs that are at least 30 min 
% separated from other PIREPs 
%---------------------------------------
time_uncorl = ones(size(unix_time));
for i=2:length(unix_time)
  if time_uncorl(i-1) == 0
    time_diff(i) = unix_time(i) - unix_time(i-2);
  else
    time_diff(i) = unix_time(i) - unix_time(i-1);
  end
  if time_diff(i) < 1800
    time_uncorl(i) = 0;
  end
end

ind = find(time_uncorl == 1);

%plot(severity(ind))

%---------------------------------------
% convert unix time to yyyymmddhhmmss 
%---------------------------------------
for j=1:length(unix_time(ind))
  yyyymmddhhmmss_str{j}=date_time_utilities('yyyymmddhhmmss_from_unix_time',{unix_time(ind(j))});
end
%yyyymmddhhmmss = cell2mat(yyyymmddhhmmss_str{1});

%---------------------------------------
% prepare to find matching radiometer & radar files  
%---------------------------------------
cd /d1/fripp_d3/nirss/platteville/data/radiometer
% uncomment below if need to redo radiom file list
%unix('ls *lv2.csv > radiom_filelist.txt');
fid = fopen('radiom_filelist.txt');
radiom_file_list     = textscan(fid,'%s');
radiom_file_list_str = cell2mat(radiom_file_list{1});
fclose(fid);

for ll=1:length(radiom_file_list_str)
  radiom_file_date(ll) = str2num([radiom_file_list_str(ll,1:4) ...
                         radiom_file_list_str(ll,6:7) radiom_file_list_str(ll,9:10)]);
end

%---------------------------------------
% loop through each PIREP date & time
% and find matching radiom and radar data  
%---------------------------------------
%for ii=200:length(severity(ind))
for ii=1:130
%for ii=1:length(severity(ind))
  PIREP_date_str = [yyyymmddhhmmss_str{ii}(1:4),yyyymmddhhmmss_str{ii}(5:6),...
                    yyyymmddhhmmss_str{ii}(7:8)];
  PIREP_date     = str2num(PIREP_date_str);
  ind2           = find(radiom_file_date == PIREP_date);

  %-----------------------------
  % if/else to see if daily radiom file exists 
  %-----------------------------
  if length(ind2) == 0
    disp('Warning: radiometer data file does not exist');
    ILW(ii) = NaN;
    ILW_30(ii) = NaN;
    ILW_60(ii) = NaN;

  %-----------------------------
  % read in radiometer file 
  % from PIREP date 
  %-----------------------------
  else
    fid2           = fopen(radiom_file_list{1}{ind2});

    radiom_data    = textscan(fid2,'%f%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','Delimiter',',');
    fclose(fid2);
  
    %-----------------------------
    % find radiometer ILW at PIREP time 
    %-----------------------------
    PIREP_time_str = [yyyymmddhhmmss_str{ii}(9:10),yyyymmddhhmmss_str{ii}(11:12)];
    for mm=1:length(radiom_data{1})
      radiom_time_str{mm}= [radiom_data{2}{mm}(10:11), radiom_data{2}{mm}(13:14)];
      if radiom_time_str{mm} == PIREP_time_str
        ILW(ii) = radiom_data{10}(mm);
        ind_mm  = mm;
      end
    end

    %-----------------------------
    % if radiom file does not include PIREP time, 
    % then set ILW(ii) = Nan 
    %-----------------------------
    if str2num(radiom_time_str{1}) > str2num(PIREP_time_str)
      disp('Warning: radiom file START time does not extend to PIREP time.');
      ILW(ii) = NaN;
      ILW_30(ii) = NaN;
      ILW_60(ii) = NaN;
    elseif str2num(radiom_time_str{end}) < str2num(PIREP_time_str)
      disp('Warning: radiom file END time does not extend to PIREP time.');
      ILW(ii) = NaN;
      ILW_30(ii) = NaN;
      ILW_60(ii) = NaN;

    else
      %-----------------------------
      % find max radiometer ILW and time_offset 
      % up to +/- 30 min from PIREP time 
      % NOTE: 120 radiom records is 30 min
      %-----------------------------
      if str2num(PIREP_time_str) < 30 
        new_range = 4*(str2num(PIREP_time_str)) - 4;
        for pp=1:120
          ILW_plus_30(pp)   = radiom_data{10}(ind_mm+pp);
        end
        for qq=1:new_range
          ILW_minus_30(qq)  = radiom_data{10}(ind_mm-qq);
        end

      elseif  str2num(PIREP_time_str) > 2330
        new_range = 4*(str2num(PIREP_time_str) - 2330);
        for pp=1:new_range
          ILW_plus_30(pp)   = radiom_data{10}(ind_mm+pp);
        end
        for qq=1:120
          ILW_minus_30(qq)  = radiom_data{10}(ind_mm-qq);
        end
      else
        for pp=1:120
          ILW_plus_30(pp)  = radiom_data{10}(ind_mm+pp); 
          ILW_minus_30(pp) = radiom_data{10}(ind_mm-pp); 
        end    % end of for pp=1:120
      end      % end of if str2num(PIREP_time_str)

      ind_plus_30      = find(ILW_plus_30  == max(ILW_plus_30));
      ind_minus_30     = find(ILW_minus_30 == max(ILW_minus_30));
      ind_plus_30      = min(ind_plus_30);
      ind_minus_30     = min(ind_minus_30);

      if max(ILW_plus_30) > max(ILW_minus_30)
        ind_time_offset_30 = ind_plus_30;
      elseif max(ILW_plus_30) < max(ILW_minus_30)
        ind_time_offset_30 = -ind_minus_30;
      elseif max(ILW_plus_30) == max(ILW_minus_30)
        if ind_plus_30 < ind_minus_30
          ind_time_offset_30 = ind_plus_30;
        elseif ind_plus_30 > ind_minus_30 
          ind_time_offset_30 = -ind_minus_30;
        elseif ind_plus_30 == ind_minus_30 
          ind_time_offset_30 = ind_plus_30;
        end
      end

      time_offset_30(ii) = ceil(ind_time_offset_30/4);

      %-----------------------------
      % find max radiometer ILW and time_offset 
      % up to +/- 60 min from PIREP time 
      % NOTE: 240 radiom records is 60 min
      %-----------------------------
      if str2num(PIREP_time_str) < 101 
        if str2num(PIREP_time_str) == 100
          PIREP_time_str = '59';
        end 
        new_range_60 = 4*(str2num(PIREP_time_str)) - 4;
        for pp=1:240
          ILW_plus_60(pp)   = radiom_data{10}(ind_mm+pp);
        end
        for qq=1:new_range_60
          ILW_minus_60(qq)  = radiom_data{10}(ind_mm-qq);
        end

      elseif  str2num(PIREP_time_str) > 2300
        new_range_60= 4*(str2num(PIREP_time_str) - 2300);
        for pp=1:new_range_60
          ILW_plus_60(pp)   = radiom_data{10}(ind_mm+pp);
        end
        for qq=1:240
          ILW_minus_60(qq)  = radiom_data{10}(ind_mm-qq);
        end
      else
        for pp=1:240
          ILW_plus_60(pp)  = radiom_data{10}(ind_mm+pp); 
          ILW_minus_60(pp) = radiom_data{10}(ind_mm-pp); 
        end    % end of for pp=1:240
      end      % end of if str2num(PIREP_time_str)

      ind_plus_60      = find(ILW_plus_60  == max(ILW_plus_60));
      ind_minus_60     = find(ILW_minus_60 == max(ILW_minus_60));
      ind_plus_60      = min(ind_plus_60);
      ind_minus_60     = min(ind_minus_60);

      if max(ILW_plus_60) > max(ILW_minus_60)
        ind_time_offset_60 = ind_plus_60;
      elseif max(ILW_plus_60) < max(ILW_minus_60)
        ind_time_offset_60 = -ind_minus_60;
      elseif max(ILW_plus_60) == max(ILW_minus_60)
        if ind_plus_60 < ind_minus_60
          ind_time_offset_60 = ind_plus_60;
        elseif ind_plus_60 > ind_minus_60 
          ind_time_offset_60 = -ind_minus_60;
        elseif ind_plus_60 == ind_minus_60 
          ind_time_offset_60 = ind_plus_60;
        end
      end

      time_offset_60(ii) = ceil(ind_time_offset_60/4);

      clear mm ILW_plus_60 ILW_plus_30 ILW_minus_60 ILW_minus_30;

      %-----------------------------
      % read in K-band REFL profile
      % nearest to PIREP time 
      %-----------------------------
  
      %-----------------------------
      % read in S-band REFL profile
      % nearest to  PIREP time 
      %-----------------------------

      %-----------------------------
      % processing progress display
      %-----------------------------
      disp('------')
      disp(['PIREP #        = ', num2str(ii)])
      disp(['PIREP date     = ', PIREP_date_str])
      disp(['PIREP time     = ', PIREP_time_str])
      disp(['PIREP SEV      = ', num2str(severity(ind(ii))), '/8'])
      disp(['ILW            = ', num2str(ILW(ii))])
      disp(['time offset_30 = ', num2str(time_offset_30(ii)), ' min'])
      disp(['ILW@offset_30  = ', num2str(radiom_data{10}(ind_mm+ind_time_offset_30+1))])
      disp(['time offset_60 = ', num2str(time_offset_60(ii)), ' min'])
      disp(['ILW@offset_60  = ', num2str(radiom_data{10}(ind_mm+ind_time_offset_60+1))])
    
      %-----------------------------
      %  create ILW(ind) field 
      %-----------------------------
      ILW_30(ii) = radiom_data{10}(ind_mm+ind_time_offset_30+1);
      ILW_60(ii) = radiom_data{10}(ind_mm+ind_time_offset_60+1);

      %keyboard;

      %-----------------------------
      % save output fields to .mat file 
      %-----------------------------
      % save time_PIREP time_offset severity ILW REFL_prof_K REFL_prof_S /dir/filename.mat....

      %-----------------------------
      % clear fields for next pass thru loop 
      %-----------------------------
      clear radiom_data radiom_time_str ind2 ind_mm ind_time_offset;
   
    end % end if radiom file time does not extend to PIREP time

    clear radiom_data;

  end   % end of if radiometer date file exists

end     % end of for ii=1:length of PIREP list
    
      %-----------------------------
      % put ILW and PIREP SEV into bins 
      %-----------------------------
      % Key for NIRSS ILW (n) bins:
      %   n0 = 0.00 to 0.02 
      %   n1 = 0.03 to 0.09
      %   n2 = 0.10 to 0.14
      %   n3 = 0.15 to 0.19
      %   n4 = 0.20 to 0.24
      %   n5 = 0.25 to 0.29
      %   n6 = 0.30 to 0.34
      %   n7 = 0.30 to 0.49
      %   n8 = 0.5 and above
      bin_p0_n0 = 0;bin_p1_n0 = 0;bin_p2_n0 = 0;bin_p3_n0 = 0;bin_p4_n0 = 0;bin_p5_n0 = 0;bin_p6_n0 = 0;bin_p7_n0 = 0;bin_p8_n0 = 0;
      bin_p0_n1 = 0;bin_p1_n1 = 0;bin_p2_n1 = 0;bin_p3_n1 = 0;bin_p4_n1 = 0;bin_p5_n1 = 0;bin_p6_n1 = 0;bin_p7_n1 = 0;bin_p8_n1 = 0;
      bin_p0_n2 = 0;bin_p1_n2 = 0;bin_p2_n2 = 0;bin_p3_n2 = 0;bin_p4_n2 = 0;bin_p5_n2 = 0;bin_p6_n2 = 0;bin_p7_n2 = 0;bin_p8_n2 = 0;
      bin_p0_n3 = 0;bin_p1_n3 = 0;bin_p2_n3 = 0;bin_p3_n3 = 0;bin_p4_n3 = 0;bin_p5_n3 = 0;bin_p6_n3 = 0;bin_p7_n3 = 0;bin_p8_n3 = 0;
      bin_p0_n4 = 0;bin_p1_n4 = 0;bin_p2_n4 = 0;bin_p3_n4 = 0;bin_p4_n4 = 0;bin_p5_n4 = 0;bin_p6_n4 = 0;bin_p7_n4 = 0;bin_p8_n4 = 0;
      bin_p0_n5 = 0;bin_p1_n5 = 0;bin_p2_n5 = 0;bin_p3_n5 = 0;bin_p4_n5 = 0;bin_p5_n5 = 0;bin_p6_n5 = 0;bin_p7_n5 = 0;bin_p8_n5 = 0;
      bin_p0_n6 = 0;bin_p1_n6 = 0;bin_p2_n6 = 0;bin_p3_n6 = 0;bin_p4_n6 = 0;bin_p5_n6 = 0;bin_p6_n6 = 0;bin_p7_n6 = 0;bin_p8_n6 = 0;
      bin_p0_n7 = 0;bin_p1_n7 = 0;bin_p2_n7 = 0;bin_p3_n7 = 0;bin_p4_n7 = 0;bin_p5_n7 = 0;bin_p6_n7 = 0;bin_p7_n7 = 0;bin_p8_n7 = 0;
      bin_p0_n8 = 0;bin_p1_n8 = 0;bin_p2_n8 = 0;bin_p3_n8 = 0;bin_p4_n8 = 0;bin_p5_n8 = 0;bin_p6_n8 = 0;bin_p7_n8 = 0;bin_p8_n8 = 0; 

    %for jj = 1:length(severity(ind))   
    for jj = 1:length(ILW_30)   
      if severity(ind(jj)) < 0.5
        if ILW_30(jj) <= 0.02     
           bin_p0_n0 = bin_p0_n0 + 1;
        elseif ILW_30(jj) < 0.09 & ILW_30(jj) > 0.02  
           bin_p0_n1 = bin_p0_n1 + 1;
        elseif ILW_30(jj) < 0.14 & ILW_30(jj) >= 0.10     
           bin_p0_n2 = bin_p0_n2 + 1;
        elseif ILW_30(jj) < 0.19 & ILW_30(jj) >= 0.15  
           bin_p0_n3 = bin_p0_n3 + 1;
        elseif ILW_30(jj) < 0.24 & ILW_30(jj) >= 0.20  
           bin_p0_n4 = bin_p0_n4 + 1;
        elseif ILW_30(jj) < 0.29 & ILW_30(jj) >= 0.25 
           bin_p0_n5 = bin_p0_n5 + 1;
        elseif ILW_30(jj) < 0.34 & ILW_30(jj) >= 0.30        
           bin_p0_n6 = bin_p0_n6 + 1;
        elseif ILW_30(jj) < 0.49 & ILW_30(jj) >= 0.49    
           bin_p0_n7 = bin_p0_n7 + 1;
        elseif ILW_30(jj) >= 0.50    
           bin_p0_n8 = bin_p0_n8 + 1;
        end
      elseif severity(ind(jj)) < 1.5 & severity(ind(jj)) >= 0.5
        if ILW_30(jj) <= 0.02     
           bin_p1_n0 = bin_p1_n0 + 1;
        elseif ILW_30(jj) < 0.09 & ILW_30(jj) > 0.02  
           bin_p1_n1 = bin_p1_n1 + 1;
        elseif ILW_30(jj) < 0.14 & ILW_30(jj) >= 0.10     
           bin_p1_n2 = bin_p1_n2 + 1;
        elseif ILW_30(jj) < 0.19 & ILW_30(jj) >= 0.15  
           bin_p1_n3 = bin_p1_n3 + 1;
        elseif ILW_30(jj) < 0.24 & ILW_30(jj) >= 0.20  
           bin_p1_n4 = bin_p1_n4 + 1;
        elseif ILW_30(jj) < 0.29 & ILW_30(jj) >= 0.25 
           bin_p1_n5 = bin_p1_n5 + 1;
        elseif ILW_30(jj) < 0.34 & ILW_30(jj) >= 0.30        
           bin_p1_n6 = bin_p1_n6 + 1;
        elseif ILW_30(jj) < 0.49 & ILW_30(jj) >= 0.49    
           bin_p1_n7 = bin_p1_n7 + 1;
        elseif ILW_30(jj) >= 0.50    
           bin_p1_n8 = bin_p1_n8 + 1;
        end
      elseif severity(ind(jj)) < 2.5 & severity(ind(jj)) >= 1.5
        if ILW_30(jj) <= 0.02     
           bin_p2_n0 = bin_p2_n0 + 1;
        elseif ILW_30(jj) < 0.09 & ILW_30(jj) > 0.02  
           bin_p2_n1 = bin_p2_n1 + 1;
        elseif ILW_30(jj) < 0.14 & ILW_30(jj) >= 0.10     
           bin_p2_n2 = bin_p2_n2 + 1;
        elseif ILW_30(jj) < 0.19 & ILW_30(jj) >= 0.15  
           bin_p2_n3 = bin_p2_n3 + 1;
        elseif ILW_30(jj) < 0.24 & ILW_30(jj) >= 0.20  
           bin_p2_n4 = bin_p2_n4 + 1;
        elseif ILW_30(jj) < 0.29 & ILW_30(jj) >= 0.25 
           bin_p2_n5 = bin_p2_n5 + 1;
        elseif ILW_30(jj) < 0.34 & ILW_30(jj) >= 0.30        
           bin_p2_n6 = bin_p2_n6 + 1;
        elseif ILW_30(jj) < 0.49 & ILW_30(jj) >= 0.49    
           bin_p2_n7 = bin_p2_n7 + 1;
        elseif ILW_30(jj) >= 0.50    
           bin_p2_n8 = bin_p2_n8 + 1;
        end
      elseif severity(ind(jj)) < 3.5 & severity(ind(jj)) >= 2.5
        if ILW_30(jj) <= 0.02     
           bin_p3_n0 = bin_p3_n0 + 1;
        elseif ILW_30(jj) < 0.09 & ILW_30(jj) > 0.02  
           bin_p3_n1 = bin_p3_n1 + 1;
        elseif ILW_30(jj) < 0.14 & ILW_30(jj) >= 0.10     
           bin_p3_n2 = bin_p3_n2 + 1;
        elseif ILW_30(jj) < 0.19 & ILW_30(jj) >= 0.15  
           bin_p3_n3 = bin_p3_n3 + 1;
        elseif ILW_30(jj) < 0.24 & ILW_30(jj) >= 0.20  
           bin_p3_n4 = bin_p3_n4 + 1;
        elseif ILW_30(jj) < 0.29 & ILW_30(jj) >= 0.25 
           bin_p3_n5 = bin_p3_n5 + 1;
        elseif ILW_30(jj) < 0.34 & ILW_30(jj) >= 0.30        
           bin_p3_n6 = bin_p3_n6 + 1;
        elseif ILW_30(jj) < 0.49 & ILW_30(jj) >= 0.49    
           bin_p3_n7 = bin_p3_n7 + 1;
        elseif ILW_30(jj) >= 0.50    
           bin_p3_n8 = bin_p3_n8 + 1;
        end
      elseif severity(ind(jj)) < 4.5 & severity(ind(jj)) >= 3.5
        if ILW_30(jj) <= 0.02     
           bin_p4_n0 = bin_p4_n0 + 1;
        elseif ILW_30(jj) < 0.09 & ILW_30(jj) > 0.02  
           bin_p4_n1 = bin_p4_n1 + 1;
        elseif ILW_30(jj) < 0.14 & ILW_30(jj) >= 0.10     
           bin_p4_n2 = bin_p4_n2 + 1;
        elseif ILW_30(jj) < 0.19 & ILW_30(jj) >= 0.15  
           bin_p4_n3 = bin_p4_n3 + 1;
        elseif ILW_30(jj) < 0.24 & ILW_30(jj) >= 0.20  
           bin_p4_n4 = bin_p4_n4 + 1;
        elseif ILW_30(jj) < 0.29 & ILW_30(jj) >= 0.25 
           bin_p4_n5 = bin_p4_n5 + 1;
        elseif ILW_30(jj) < 0.34 & ILW_30(jj) >= 0.30        
           bin_p4_n6 = bin_p4_n6 + 1;
        elseif ILW_30(jj) < 0.49 & ILW_30(jj) >= 0.49    
           bin_p4_n7 = bin_p4_n7 + 1;
        elseif ILW_30(jj) >= 0.50    
           bin_p4_n8 = bin_p4_n8 + 1;
        end
      elseif severity(ind(jj)) < 5.5 & severity(ind(jj)) >= 4.5
        if ILW_30(jj) <= 0.02     
           bin_p5_n0 = bin_p5_n0 + 1;
        elseif ILW_30(jj) < 0.09 & ILW_30(jj) > 0.02  
           bin_p5_n1 = bin_p5_n1 + 1;
        elseif ILW_30(jj) < 0.14 & ILW_30(jj) >= 0.10     
           bin_p5_n2 = bin_p5_n2 + 1;
        elseif ILW_30(jj) < 0.19 & ILW_30(jj) >= 0.15  
           bin_p5_n3 = bin_p5_n3 + 1;
        elseif ILW_30(jj) < 0.24 & ILW_30(jj) >= 0.20  
           bin_p5_n4 = bin_p5_n4 + 1;
        elseif ILW_30(jj) < 0.29 & ILW_30(jj) >= 0.25 
           bin_p5_n5 = bin_p5_n5 + 1;
        elseif ILW_30(jj) < 0.34 & ILW_30(jj) >= 0.30        
           bin_p5_n6 = bin_p5_n6 + 1;
        elseif ILW_30(jj) < 0.49 & ILW_30(jj) >= 0.49    
           bin_p5_n7 = bin_p5_n7 + 1;
        elseif ILW_30(jj) >= 0.50    
           bin_p5_n8 = bin_p5_n8 + 1;
        end
      elseif severity(ind(jj)) < 6.5 & severity(ind(jj)) >= 5.5
        if ILW_30(jj) <= 0.02     
           bin_p6_n0 = bin_p6_n0 + 1;
        elseif ILW_30(jj) < 0.09 & ILW_30(jj) > 0.02  
           bin_p6_n1 = bin_p6_n1 + 1;
        elseif ILW_30(jj) < 0.14 & ILW_30(jj) >= 0.10     
           bin_p6_n2 = bin_p6_n2 + 1;
        elseif ILW_30(jj) < 0.19 & ILW_30(jj) >= 0.15  
           bin_p6_n3 = bin_p6_n3 + 1;
        elseif ILW_30(jj) < 0.24 & ILW_30(jj) >= 0.20  
           bin_p6_n4 = bin_p6_n4 + 1;
        elseif ILW_30(jj) < 0.29 & ILW_30(jj) >= 0.25 
           bin_p6_n5 = bin_p6_n5 + 1;
        elseif ILW_30(jj) < 0.34 & ILW_30(jj) >= 0.30        
           bin_p6_n6 = bin_p6_n6 + 1;
        elseif ILW_30(jj) < 0.49 & ILW_30(jj) >= 0.49    
           bin_p6_n7 = bin_p6_n7 + 1;
        elseif ILW_30(jj) >= 0.50    
           bin_p6_n8 = bin_p6_n8 + 1;
        end
      elseif severity(ind(jj)) < 7.5 & severity(ind(jj)) >= 6.5
        if ILW_30(jj) <= 0.02     
           bin_p7_n0 = bin_p7_n0 + 1;
        elseif ILW_30(jj) < 0.09 & ILW_30(jj) > 0.02  
           bin_p7_n1 = bin_p7_n1 + 1;
        elseif ILW_30(jj) < 0.14 & ILW_30(jj) >= 0.10     
           bin_p7_n2 = bin_p7_n2 + 1;
        elseif ILW_30(jj) < 0.19 & ILW_30(jj) >= 0.15  
           bin_p7_n3 = bin_p7_n3 + 1;
        elseif ILW_30(jj) < 0.24 & ILW_30(jj) >= 0.20  
           bin_p7_n4 = bin_p7_n4 + 1;
        elseif ILW_30(jj) < 0.29 & ILW_30(jj) >= 0.25 
           bin_p7_n5 = bin_p7_n5 + 1;
        elseif ILW_30(jj) < 0.34 & ILW_30(jj) >= 0.30        
           bin_p7_n6 = bin_p7_n6 + 1;
        elseif ILW_30(jj) < 0.49 & ILW_30(jj) >= 0.49    
           bin_p7_n7 = bin_p7_n7 + 1;
        elseif ILW_30(jj) >= 0.50    
           bin_p7_n8 = bin_p7_n8 + 1;
        end
      elseif severity(ind(jj)) >= 7.5
        if ILW_30(jj) <= 0.02     
           bin_p8_n0 = bin_p8_n0 + 1;
        elseif ILW_30(jj) < 0.09 & ILW_30(jj) > 0.02  
           bin_p8_n1 = bin_p8_n1 + 1;
        elseif ILW_30(jj) < 0.14 & ILW_30(jj) >= 0.10     
           bin_p8_n2 = bin_p8_n2 + 1;
        elseif ILW_30(jj) < 0.19 & ILW_30(jj) >= 0.15  
           bin_p8_n3 = bin_p8_n3 + 1;
        elseif ILW_30(jj) < 0.24 & ILW_30(jj) >= 0.20  
           bin_p8_n4 = bin_p8_n4 + 1;
        elseif ILW_30(jj) < 0.29 & ILW_30(jj) >= 0.25 
           bin_p8_n5 = bin_p8_n5 + 1;
        elseif ILW_30(jj) < 0.34 & ILW_30(jj) >= 0.30        
           bin_p8_n6 = bin_p8_n6 + 1;
        elseif ILW_30(jj) < 0.49 & ILW_30(jj) >= 0.49    
           bin_p8_n7 = bin_p8_n7 + 1;
        elseif ILW_30(jj) >= 0.50    
           bin_p8_n8 = bin_p8_n8 + 1;
        end
      end
    end

    bins_p_vs_n = [bin_p0_n8 bin_p1_n8 bin_p2_n8 bin_p3_n8 bin_p4_n8 bin_p5_n8 bin_p6_n8 bin_p7_n8 bin_p8_n8; bin_p0_n7 bin_p1_n7 bin_p2_n7 bin_p3_n7 bin_p4_n7 bin_p5_n7 bin_p6_n7 bin_p7_n7 bin_p8_n7; bin_p0_n6 bin_p1_n6 bin_p2_n6 bin_p3_n6 bin_p4_n6 bin_p5_n6 bin_p6_n6 bin_p7_n6 bin_p8_n6; bin_p0_n5 bin_p1_n5 bin_p2_n5 bin_p3_n5 bin_p4_n5 bin_p5_n5 bin_p6_n5 bin_p7_n5 bin_p8_n5; bin_p0_n4 bin_p1_n4 bin_p2_n4 bin_p3_n4 bin_p4_n4 bin_p5_n4 bin_p6_n4 bin_p7_n4 bin_p8_n4; bin_p0_n3 bin_p1_n3 bin_p2_n3 bin_p3_n3 bin_p4_n3 bin_p5_n3 bin_p6_n3 bin_p7_n3 bin_p8_n3; bin_p0_n2 bin_p1_n2 bin_p2_n2 bin_p3_n2 bin_p4_n2 bin_p5_n2 bin_p6_n2 bin_p7_n2 bin_p8_n2; bin_p0_n1 bin_p1_n1 bin_p2_n1 bin_p3_n1 bin_p4_n1 bin_p5_n1 bin_p6_n1 bin_p7_n1 bin_p8_n1; bin_p0_n0 bin_p1_n0 bin_p2_n0 bin_p3_n0 bin_p4_n0 bin_p5_n0 bin_p6_n0 bin_p7_n0 bin_p8_n0];


