%---------------------------------------------------------
% integrate_sonde_model_lwcs.m
% 6.8.2015 dserke
%---------------------------------------------------------

%%---------------------------------------------------------
%% open sonde data file with no header and with fields separated by commas
%%---------------------------------------------------------
%fid1  = fopen('/d1/serke/projects/SLWprobe_NASA/data/NASA_early2015/20150129_HRRR.csv');

%%---------------------------------------------------------
%% read sonde fields into 'data' with formatting
%%---------------------------------------------------------
%data = textscan(fid1, '%5.1f%10.6f%8.6f%8.6f%8.6f8.6%f%8.6f%8.6f%9.6f%10.6f%8.6f','Delimiter',',');

%%---------------------------------------------------------
%% close the sonde input file
%%---------------------------------------------------------
%fclose(fid1);

%---------------------------------------------------------
% open model data file with no header and with fields separated by commas
%---------------------------------------------------------
fid2  = fopen('/d1/serke/projects/SLWprobe_NASA/data/NASA_early2015/20150129_HRRR.csv');

%---------------------------------------------------------
% read model fields into 'model_data' with formatting
%---------------------------------------------------------
model_data = textscan(fid2, '%5.1f%10.6f%8.6f%7.5f%7.5f%7.5f%8.6f%8.6f%9.6f%12.6f%8.6f','Delimiter',',');

%---------------------------------------------------------
% close the input file
%---------------------------------------------------------
fclose(fid2);

%---------------------------------------------------------
% calculate the integrated liquid water from sonde profile
% assume density of water is 1000kg/m3
% ILW units are mm
%---------------------------------------------------------
model_ILW = 0;
%for i = 2:13
for i = 2:size(model_data{10},1)
  model_ILW = model_ILW + model_data{8}(i) * (model_data{10}(i)-model_data{10}(i-1)) * (1000/1000/1000)
end

