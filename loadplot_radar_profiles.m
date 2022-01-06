%----------------------------------------------
%
%
%
%
%
%
%
%
%
%----------------------------------------------

%----------------------------------------------
% load NIRSS Ka-band .txt radar profile from ... 
%----------------------------------------------
%k_data = load('/d1/fripp_d3/nirss/platteville/data/radar_output_profiles/20101207180500.txt');
%load('/d1/serke/projects/NIRSS_NASA/data/20110518_2054_kband_REFL_nohdr.txt');
k_data = load('/d1/serke/projects/NIRSS_NASA/data/20110518_1954_kband_REFL_nohdr.txt');
%k_data = load('/d1/fripp_d3/nirss/platteville/data/radar_output_profiles/20110518215400.txt');

%----------------------------------------------
% plot NIRSS Ka-band .txt radar profile  
%----------------------------------------------
figure;
plot(k_data(:,3) - 30,k_data(:,2),'b-');
grid on;


%----------------------------------------------
% pre-load CHILL RHI .nc radar profile from external drive
%----------------------------------------------
%s_file = '/media/nasa/CHILL@Platteville_2010/NC/2010_1207/cfrad.20101207_174813_000_v1_PPI.nc';
%s_file = '/media/nasa/CHILL@Platteville_2010/NC/20110518/cfrad.20110518_215403_000_v3_s02_e251.92_RHI.nc';
s_file = '/media/nasa/CHILL@Platteville_2010/NC/20110518/cfrad.20110518_195119_000_v6_s05_e209.03_RHI.nc';
s_data = ncload(s_file);

%----------------------------------------------
% at what range do you want to view from CHILL radar?
%   NOTE: CHILL is 30 km distance from Platteville
%----------------------------------------------
x = 200;
 
%----------------------------------------------
% handle RHI or PPI S-band scans 
%----------------------------------------------
%scan_type = 'PPI';
scan_type = s_file(end-5:end-3);

if scan_type == 'PPI'

  platte_azim = 196;
  ind         = find(ceil(azimuth) == platte_azim);

  alt     = [];
  x_dist  = [];
  for ii=length(ind)
    alt(ii)    = [alt range(x)*tan(elevation(ii)*pi/180)];
    x_dist(ii) = [x_dist range(ii)];
  end

elseif scan_type == 'RHI'

else
  disp('Error: Bad scan_type');
end

%----------------------------------------------
% create an array of altitudes AGL
%----------------------------------------------
alt = range(x)*tan(elevation*pi/180);
%alt = alt(102:end);

%----------------------------------------------
% plot NIRSS Ka-band (and S-band) radar profile  
%----------------------------------------------
figure;
plot(k_data(:,3) - 30,k_data(:,2),'b-');
%plot(X20110518_1954_kband_REFL_nohdr(:,3),X20110518_1954_kband_REFL_nohdr(:,2),'b-');
hold on;
plot(DZ(:,x)/1000,alt,'k-');
grid on;

%----------------------------------------------
% plot S-band moment profile data
%----------------------------------------------
figure;
subplot(4,1,1);
grid on;
hold on;
plot(RH(:,x-4)/20000,alt,'r-');
plot(RH(:,x-3)/20000,alt,'y-');
plot(RH(:,x-2)/20000,alt,'Color',[0.6 0 0]);
plot(RH(:,x-1)/20000,alt,'Color',[0 0.6 0]);
plot(RH(:,x+1)/20000,alt,'Color',[0 0 0.6]);
plot(RH(:,x+2)/20000,alt,'g-');
plot(RH(:,x+3)/20000,alt,'b-');
plot(RH(:,x+4)/20000,alt,'c-');
plot(RH(:,x)/20000,alt,'k-');
axis([-2 1 0 9100]);
subplot(4,1,2);
grid on;
hold on;
plot(DZ(:,x-4)/1000,alt,'r-');
plot(DZ(:,x-3)/1000,alt,'y-');
plot(DZ(:,x-2)/1000,alt,'Color',[0.6 0 0]);
plot(DZ(:,x-1)/1000,alt,'Color',[0 0.6 0]);
plot(DZ(:,x+1)/1000,alt,'Color',[0 0 0.6]);
plot(DZ(:,x+2)/1000,alt,'g-');
plot(DZ(:,x+3)/1000,alt,'b-');
plot(DZ(:,x)/1000,alt,'k-');
axis([-32 32 0 9100]);
%subplot(5,1,3);
%grid on;
%hold on;
%axis([-32 32 0 10000]);
%plot(LH(:,x-4)/1000,alt,'r-');
%plot(LH(:,x-3)/1000,alt,'y-');
%plot(LH(:,x-2)/1000,alt,'Color',[0.6 0 0]);
%plot(LH(:,x-1)/1000,alt,'Color',[0 0.6 0]);
%plot(LH(:,x+1)/1000,alt,'Color',[0 0 0.6]);
%plot(LH(:,x+2)/1000,alt,'g-');
%plot(LH(:,x+3)/1000,alt,'b-');
%plot(LH(:,x)/1000,alt,'k-');
%axis([-32 32 0 10000]);
subplot(4,1,3);
grid on;
hold on;
plot((DR(:,x-4)/10000)+1.51,alt,'r-');
plot((DR(:,x-3)/10000)+1.51,alt,'y-');
plot((DR(:,x-2)/10000)+1.51,alt,'Color',[0.6 0 0]);
plot((DR(:,x-1)/10000)+1.51,alt,'Color',[0 0.6 0]);
plot((DR(:,x+1)/10000)+1.51,alt,'Color',[0 0 0.6]);
plot((DR(:,x+2)/10000)+1.51,alt,'g-');
plot((DR(:,x+3)/10000)+1.51,alt,'b-');
plot((DR(:,x)/10000)+1.51,alt,'k-');
axis([-2.5 2.5 0 9100]);
subplot(4,1,4);
DR_0=DR(:,x)/10000;
DR_p1=DR(:,x+1)/10000;
DR_p2=DR(:,x+2)/10000;
DR_p3=DR(:,x+3)/10000;
DR_p4=DR(:,x+4)/10000;
DR_m1=DR(:,x-1)/10000;
DR_m2=DR(:,x-2)/10000;
DR_m3=DR(:,x-3)/10000;
DR_m4=DR(:,x-4)/10000;

for i =1:size(DR,1)
  DR_STD(i) = std([DR_0(i) DR_p1(i) DR_p2(i) DR_p3(i) DR_p4(i) DR_m1(i) DR_m2(i) DR_m3(i) DR_m4(i)]); 
end

plot(DR_STD,alt)
axis([0 1.3 0 9100]);
grid on;
