%---------------------------------------------------------
% integrate_sonde_lwcs.m
% 6.8.2015 dserke
%---------------------------------------------------------

%---------------------------------------------------------
% open sonde data file with no header and with fields separated by commas
%---------------------------------------------------------
fid  = fopen('/d1/serke/projects/case_studies/winter_2014_2015_cases/cases_with_sondes/2015_0122/data/20150122001.nohdr.txt');

%---------------------------------------------------------
% read sonde fields into 'data' with formatting
%---------------------------------------------------------
data = textscan(fid, '%6.0f%7.4f%8.4f%7.1f%7.2f%6.2f%6.2f%6.2f%6.3f%5.1f%4.0f%4.0f%5.2f%5.2f%7.6f%7.6f','Delimiter',',');

%---------------------------------------------------------
% close the input file
%---------------------------------------------------------
fclose(fid);

%---------------------------------------------------------
% convert hhmmss to total seconds into the day
%---------------------------------------------------------
time_in_seconds = [];

for j = 1:size(data{16},1)
  time_char       = num2str(data{1}(j));
  hour_to_sec     = str2num(time_char(1:2)) * 3600;
  min_to_sec      = str2num(time_char(3:4)) * 60;
  seconds         = str2num(time_char(5:6));
  time_in_sec     = hour_to_sec + min_to_sec + seconds;
  time_in_seconds = [time_in_seconds time_in_sec];
end

keyboard;

%---------------------------------------------------------
% calculate the integrated liquid water from sonde profile
% assume density of water is 1000kg/m3
% ILW units are mm
%---------------------------------------------------------
ILW = 0;
%for i = 2:13
for i = 2:size(data{16},1)
  ILW = ILW + data{16}(i)*data{9}(i)*(time_in_seconds(i) - time_in_seconds(i-1))*1000/1000/1000;
end
