%------------------------------------
% loads and plots .csv format data from
% NASA's Narrowbeam Multichannel Microwave Radiometer
%
%  created 2.10.2011 - dserke
%------------------------------------

% first need to:
%   open .xls file in openoffice
%   save only data sheet to new file
%   rm header line
%   %s/"// until all quotes are gone

%------------------------------------
% load NNMMR csv file
%------------------------------------
a=readtext('/d1/serke/projects/NIRSS_NASA/data/NNMMR/NNMMRFeb1-3_data_all_nohdr.csv','[,]','%','"','numeric-empty2zero');
%a=readtext('NNMMRFeb1-3_data1_nohdr.csv','[,]','%','"','numeric-empty2zero');

line_num   = a(:,1);
date_time  = a(:,2);
az_angle   = a(:,4);
el_angle   = a(:,5);
LWP        = a(:,30);

count = 1;
for j=1:313
  for i=1:19
    LWP_tvsh(i,j) = LWP(count);
    count = count + 1;
  end
end

%------------------------------------
% make PIREP database
%------------------------------------

PIREP_time = [1 6 22 101 122 137 166 178];
PIREP_cat  = [0 0 4 4 5 5 3 3];

%------------------------------------
% present weather database 
%------------------------------------

wx_time  = [6 14 19 26.5 33.8 41 48.1 55.2 62.6 69.8 77.1 ...
            84.4 91.7 99 106.3 113.6 120.9 128.2 135.5 ...
            142.8 150.1 157.4 164.7 172 179.3 186.6 193.9];
wx_type  = ['s' ' ' 'z' 'z' 'z' 'z' 'z' 'z' 'z' 'z' 'r' 'r' ...
            ' ' ' ' 'r' 's' 'p' 's' 's' 's' 's' 's' 'S' ' ' ...
            's' 's' 's'];
%------------------------------------
% plotting code.... 
%------------------------------------
x=1:200;
y=-90:5:0;

imagesc(x,y,flipud(LWP_tvsh(:,1:200)),[0 0.8])
colorbar;
hold on;
for ll=1:length(PIREP_time)
  text(PIREP_time(ll), -85, num2str(PIREP_cat(ll)), 'Color','r');
end
for mm=1:length(wx_time)
  text(wx_time(mm), 0, num2str(wx_type(mm)), 'Color','r');
end
title('20110201 21:20 to 20110203 1:00 UTC NNMMR LWP @ 151 deg AZ')
ylabel('EL angle [deg]')
xlabel('Profile [#]')
