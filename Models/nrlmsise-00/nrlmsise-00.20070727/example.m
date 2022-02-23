IDAY = 172;
SEC = 1.5*3600;
ALT = 150;
GLAT = 60;
GLONG = 90; %360-70;
F107A = 150;
F107 = 150;
AP = 4*ones(7,1);  
MASS = 48;   % 
% altitudes for the profile
AltArray_km = [0:1:99 100:10:190 200:50:1000];

myNRLMSISE00 = nrlmsise00(IDAY,SEC,ALT,GLAT,GLONG,F107A,F107,AP,MASS);
calculated = myNRLMSISE00.calculateProfile(AltArray_km);

h = figure('color','white');
plot(calculated.O,calculated.Alt,'DisplayName','O'); hold on;
plot(calculated.N2,calculated.Alt,'DisplayName','N2');
plot(calculated.O2,calculated.Alt,'DisplayName','O2');
plot(calculated.He,calculated.Alt,'DisplayName','He');
plot(calculated.Ar,calculated.Alt,'DisplayName','Ar');
plot(calculated.H,calculated.Alt,'DisplayName','H');
plot(calculated.N,calculated.Alt,'DisplayName','N');
plot(calculated.OAnom,calculated.Alt,':','DisplayName','Anomalous Oxygen');
xlim([1e-10,1e20]);
set(gca,'xscale','log');
xlabel('Density [particles/cm^3]');
ylabel('Altitude [km]');
UD = calculated.Properties.UserData;
title({['\phi=' num2str(UD.GLAT) '\circN \lambda=' num2str(UD.GLONG) '\circE DOY '  ...
     num2str(UD.DAY) ' ' num2str(UD.SEC/3600) ' UTC LT = ' num2str(UD.LOCALSOLARTIME) ' Hr']; ...
    ['F10.7 = ' num2str(UD.F107) ' F10.7A = ' num2str(UD.F107A) ' sfu ' 'AP(1)=' num2str(AP(1)) ' 2nT']});
legend('location','northeast');