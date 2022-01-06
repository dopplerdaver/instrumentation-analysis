%---------------------------------------
% Name: plot_3085radiom.m
%
% Purpose: 
%
% Usage:
%
% Created: 10.9.2012 dserke
%
%
%---------------------------------------


%---------------------
% set some important fields 
%---------------------
hgt = [0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.25, 2.50, 2.75, 3.00, 3.25, 3.50, 3.75, 4.00, 4.25, 4.50, 4.75, 5.00, 5.25, 5.50, 5.75, 6.00, 6.25, 6.50, 6.75, 7.00, 7.25, 7.50, 7.75, 8.00, 8.25, 8.50, 8.75, 9.00, 9.25, 9.50, 9.75,10.00]; 

lines_per_volume = 35;

%---------------------
% find number of lines in the file 
%---------------------
cd /d1/serke/projects/NIRSS_NASA/data/volumetric_prod_case_data/;
fid1 = fopen([case_date  '.3085']);

first_instance = 1;
for kk=1:20000
  tline = fgetl(fid1);
  if tline == -1 & first_instance == 1
    lines_per_file = kk-1;
    first_instance = 0;
  end
end
fclose(fid1);

%---------------------
% find number of complete volumes in the file 
%---------------------
num_volumes = floor(lines_per_file/lines_per_volume);

%---------------------
% skip header lines
%---------------------
fid1 = fopen([case_date  '.3085']);

% read first text line
tline = fgetl(fid1);

% check to see if header lines have been removed
if tline(1) == 'R'
  % header lines have not been removed
  % advance through 32 more header lines
  for mm=1:32
    tline = fgetl(fid1);
  end
else
  % header lines have been removed
  % close fid1 and reopen fid1 at the beginning
  fclose(fid1);
  fid1 = fopen([case_date  '.3085']);
end

%---------------------
% read in data 
%---------------------
%for i = 1:1
for i = 1:num_volumes-1

  %for j = 1:17
  for j = 1:lines_per_volume - 2

    %------------
    % read next line 
    %------------
    tline     = fgetl(fid1);

    ind       = find(tline == ',');
    class     = str2num(tline(ind(2)+1:ind(2)+3));
    processor = tline(ind(3)+1);

    if processor == 'Z'
      angle = 'Z'; 
      tilt  = 'Z';
    elseif processor == 'A'
      tilt  = tline(ind(3)+6:ind(3)+7);
      angle = tline(ind(3)+9) ;
      if tilt == 'EW' & angle == 'N'
        angle = 'E';
      elseif tilt == 'EW' & angle == 'S' 
        angle = 'W';
      end
    end

    %---------------------
    % get fields from line 
    %---------------------
    if (class == 31)
      [line_num(i),date_time(i),class,a,lat(i),lon(i),b,c,d,e,altm(i)] = ...
                         strread(tline,'%d%s%d%s%9.6f%10.5f%s%s%s%s%10.6f','delimiter',',');

    elseif (class == 201)
      [line_num(i),a,class,tempamb(i),rhamb(i),pamb(i),tir(i),rainflag(i)] = ...
                         strread(tline,'%d%s%d%8.5f%8.5f%8.5f%8.5f%d','delimiter',',');

    elseif (class == 401) & processor == 'Z'
      tline2  = tline(39:end);
      temp    = [];
      count   = 1;
      for k = 1:length(hgt)
        temp  = [temp str2num(tline2(count:count+6))];
        count = count + 6 + 2;
      end
      temp_z(i,:) = temp;
  
    elseif (class == 402) & processor == 'Z'
      tline2  = tline(39:end);
      vapd    = [];
      count   = 1;
      for k = 1:length(hgt)
        vapd  = [vapd str2num(tline2(count:count+6))];
        count = count + 6 + 2;
      end
      vapd_z(i,:) = vapd;

    elseif (class == 403)

    elseif (class == 404) & processor == 'Z'
      tline2  = tline(39:end);
      rh      = [];
      count   = 1;
      for k = 1:length(hgt)
        rh  = [rh str2num(tline2(count:count+6))];
        count = count + 6 + 2;
      end
      rh_z(i,:) = rh;

    elseif (class == 301) & angle == 'Z'
      [a,b,c,ivw_z(i),ilw_z(i),cbase_z(i)] = ...
                         strread(tline,'%s%s%s%7.4f%7.4f%7.4f','delimiter',',');

    elseif (class == 401) & processor == 'A' & angle == 'N'
      tline2  = tline(40:end);
      temp    = [];
      count   = 1;
      for k = 1:length(hgt)
        temp  = [temp str2num(tline2(count:count+6))];
        count = count + 6 + 2;
      end
      temp_n(i,:) = temp;

    elseif (class == 401) & processor == 'A' & angle == 'S'
      tline2  = tline(40:end);
      temp    = [];
      count   = 1;
      for k = 1:length(hgt)
        temp  = [temp str2num(tline2(count:count+6))];
        count = count + 6 + 2;
      end
      temp_s(i,:) = temp;

    elseif (class == 402) & processor == 'A' & angle == 'N'
      tline2  = tline(40:end);
      vapd    = [];
      count   = 1;
      for k = 1:length(hgt)
        vapd  = [vapd str2num(tline2(count:count+6))];
        count = count + 6 + 2;
      end
      vapd_n(i,:) = vapd;

    elseif (class == 402) & processor == 'A' & angle == 'S'
      tline2  = tline(40:end);
      vapd    = [];
      count   = 1;
      for k = 1:length(hgt)
        vapd  = [vapd str2num(tline2(count:count+6))];
        count = count + 6 + 2;
      end
      vapd_s(i,:) = vapd;

    elseif (class == 404) & processor == 'A' & angle == 'N'
      tline2  = tline(41:end);
      rh      = [];
      count   = 1;
      for k = 1:length(hgt)
        rh  = [rh str2num(tline2(count:count+6))];
        count = count + 6 + 2;
      end
      rh_n(i,:) = rh;

    elseif (class == 404) & processor == 'A' & angle == 'S'
      tline2  = tline(41:end);
      rh      = [];
      count   = 1;
      for k = 1:length(hgt)
        rh  = [rh str2num(tline2(count:count+6))];
        count = count + 6 + 2;
      end
      rh_s(i,:) = rh;

    elseif (class == 301) & angle == 'A' & tilt == 'NS'
      [a,b,c,ivw_n(i),ilw_n(i),cbase_n(i)] = ...
                         strread(tline,'%s%s%s%7.4f%7.4f%7.4f','delimiter',',');
      tline = fgetl(fid1);
      j     = j + 1;       
      [a,b,c,ivw_s(i),ilw_s(i),cbase_n(i)] = ...
                         strread(tline,'%s%s%s%7.4f%7.4f%7.4f','delimiter',',');

    elseif (class == 401) & processor == 'A' & angle == 'E'
      tline2  = tline(40:end);
      temp    = [];
      count   = 1;
      for k = 1:length(hgt)
        temp  = [temp str2num(tline2(count:count+6))];
        count = count + 6 + 2;
      end
      temp_e(i,:) = temp;

    elseif (class == 401) & processor == 'A' & angle == 'W'
      tline2  = tline(40:end);
      temp    = [];
      count   = 1;
      for k = 1:length(hgt)
        temp  = [temp str2num(tline2(count:count+6))];
        count = count + 6 + 2;
      end
      temp_w(i,:) = temp;

    elseif (class == 402) & processor == 'A' & angle == 'E'
      tline2  = tline(40:end);
      vapd    = [];
      count   = 1;
      for k = 1:length(hgt)
        vapd  = [vapd str2num(tline2(count:count+6))];
        count = count + 6 + 2;
      end
      vapd_e(i,:) = vapd;

    elseif (class == 402) & processor == 'A' & angle == 'W'
      tline2  = tline(40:end);
      vapd    = [];
      count   = 1;
      for k = 1:length(hgt)
        vapd  = [vapd str2num(tline2(count:count+6))];
        count = count + 6 + 2;
      end
      vapd_w(i,:) = vapd;

    elseif (class == 404) & processor == 'A' & angle == 'E'
      tline2  = tline(40:end);
      rh      = [];
      count   = 1;
      for k = 1:length(hgt)
        rh  = [rh str2num(tline2(count:count+6))];
        count = count + 6 + 2;
      end
      rh_e(i,:) = rh;

    elseif (class == 404) & processor == 'A' & angle == 'W'
      tline2  = tline(40:end);
      rh      = [];
      count   = 1;
      for k = 1:length(hgt)
        rh  = [rh str2num(tline2(count:count+6))];
        count = count + 6 + 2;
      end
      rh_w(i,:) = rh;

    elseif (class == 301) & angle == 'A' & tilt == 'EW'
      [a,b,c,ivw_e(i),ilw_e(i),cbase_e(i)] = ...
                         strread(tline,'%s%s%s%7.4f%7.4f%7.4f','delimiter',',');
      tline = fgetl(fid1);
      j     = j + 1;       
      [a,b,c,ivw_w(i),ilw_w(i),cbase_w(i)] = ...
                         strread(tline,'%s%s%s%7.4f%7.4f%7.4f','delimiter',',');


    end

    %keyboard

    clear class tline tline2 count;

  j = j + 1;

  end   % end of 'for j=1'

end   % end of 'for i=1'

 
fclose(fid1);

%------------------------------------
% setup hour array 
%------------------------------------
clear hour;
for ll=1:length(date_time)
  hour{ll} = date_time{ll}(10:11); 
end
hour = cell2mat(hour);

hr = [];
count = 1;
for ll=1:length(date_time)-1
  hr = [hr str2num(hour(count:count+1))];
  count = count + 2;
end
%hr = hr(1:395);

hgt = hgt*-1;

%------------------------------------
% plotting
%------------------------------------

%------------------------------------
% water vapor

%figure
%pcolor(((vapd_z')))
%shading interp

%------------------------------------
% relative humidity
 
%figure
%pcolor(((rh_z')))
%shading interp

%figure
%pcolor(((rh_s'+rh_n')/2))
%shading interp

%------------------------------------
% ilw 

%figure
%plot(ilw_s,'r-')
%hold on
%plot(ilw_z,'k-')
%plot(ilw_n,'b-')
%legend('Southwest','Zenith','Northeast')
%ylabel('ILW [gm-2]')
%xlabel('Radiom scan volume #')

%------------------------------------
% temp 

%figure
%imagesc(hr,hgt,((temp_s'+temp_n')/2)-(temp_z'),[-3 3])
%shading interp

%figure
%pcolor(((temp_z')))
%shading interp

cd /d1/serke/projects/NIRSS_NASA/code/mat;
