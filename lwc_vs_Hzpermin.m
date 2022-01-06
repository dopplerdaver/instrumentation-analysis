LWC      = [0.55 0.39 0.28];
Hzpermin = [0.38 0.29 0.22];

figure
plot(Hzpermin,LWC,'b*-');
hold on;
grid on;
plot(0.093,0.1,'k*');
plot(0.186,0.1,'k*');
plot(0.093,0.2,'k*');
plot(0.186,0.2,'k*');

xlabel('Hertz per min');
ylabel('Mean layer LWC [gm-3]');

