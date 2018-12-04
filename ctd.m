% reads in ctd data from mission (renamed ctd_mision#nr)
% rename .mat 
% gets downcasts
% NOT DOING THAT YET -> calculates temperature hourly average for each downcast
% exports to be used by the AZFP toolbox

clear
close all

%loads data 
% ctd_data=load('ctd_87.mat');
filename=load('ctd_id_87_delayed.mat');
StrName=fieldnames(filename);
StrName=StrName{1};

% separating the variables
time=filename.(StrName).time;
d=filename.(StrName).depth;
lat=filename.(StrName).latitude;
lon=filename.(StrName).longitude;
c=filename.(StrName).conductivity;
t=filename.(StrName).temperature;
p=filename.(StrName).pressure;
s=filename.(StrName).salinity;
rho=filename.(StrName).density;
id=filename.(StrName).profile_id;
id_time=filename.(StrName).profile_time;

% plot(time,d,'.k') 

% transforming unix time in matlab time
% converts unix time stamps (seconds since Jan 1, 1970) to
% Matlab serial date number (decimal days since Jan 1 0000).
unix_epoch=datenum(1970,1,1,0,0,0);
unix_time=time;
time=time/86400+unix_epoch;
id_time=id_time/86400+unix_epoch;

% figure,plot(time,d,'.k'),axis tight, datetick('x',0)
% close all

% figure, plot(d,'*'),set(gca,'ydir','reverse')

% find non nan indexes (they are the same for t,s,p,d...) and make time
% match the other sampling frequencies
ind=find(~isnan(t));
time1=time(ind);
id1=id(ind);

%remove nans from other variables just in case
p1=p(ind);
t1=t(ind);
s1=s(ind);
d1=d(ind);
rho1=rho(ind);
lat1=lat(ind);
lon1=lon(ind);

%% indices only valid for 20 July 2018
% ctd_time=time1(1153:20138);
% id=id1(1153:20138);
% ctd_depth=d1(1153:20138);
% ctd_lat=lat1(1153:20138);
% ctd_lon=lon1(1153:20138);

ctd_time=time1;
id=id1;
ctd_depth=d1;
ctd_lat=lat1;
ctd_lon=lon1;

%% Separates downcast from upcast by assuming that
ind_up=find(diff(ctd_depth)<0); % when diff between subsequent values is negative, it means the depth is decreasing, i.e. upcast
ind_down=find(diff(ctd_depth)>0); % when diff is positive (x2-x1), it means the depth is increasing, i.e. downcast

% downcast=d1(ind_down);
% upcast=d1(ind_up);

% calculates speed of sound, during downcasts only, using gsw toolbox
% function sound_speed = gsw_sound_speed(SA,CT,p)
% SA  =  Absolute Salinity                                        [ g/kg ]
% CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
% p   =  sea pressure                                             [ dbar ]
% sound_speed = gsw_sound_speed(s1(ind_down),t1(ind_down),p1(ind_down));

% plotting stuff
figure(1)
set(gcf,'color','w');
    plot(ctd_time(ind_down),ctd_depth(ind_down),'*r'),set(gca,'ydir','reverse'),datetick('x',21),axis tight
    hold on
    plot(ctd_time(ind_up),ctd_depth(ind_up),'ob')
    plot(ctd_time,ctd_depth,'.k'),set(gca,'ydir','reverse'),datetick('x',21),axis tight
    hold off

% figure(2)
% set(gcf,'color','w');
%     plot(ctd_time(ind_down),ctd_depth(ind_down),'*r'),set(gca,'ydir','reverse'),datetick('x',0)
%     hold on
%     plot(ctd_time(ind_up),ctd_depth(ind_up),'*b')
%     plot(ctd_time,ctd_depth,'-.k')
%     yyaxis right
%     plot(ctd_time(ind_down),sound_speed,'.g'),set(gca,'ydir','reverse'),datetick('x',0)
%     hold off

% which of these has greater variability (s,t,d) and how they influence the
% sound speed? temperature has the greatest variability and influence on
% speed of sound

% figure(3)
% set(gcf,'color','w');
% %     subplot(4,1,1)
% %     plot(rho1(ind_down),d1(ind_down),'.'),set(gca,'ydir','reverse')
% %     xlabel('Density (kg.m^{-3})'),ylabel('Depth (m)')
% %     grid on
% 
%     subplot(1,4,1)
%     plot(rho1(ind_down),d1(ind_down),'.'),set(gca,'ydir','reverse')
%     xlabel('Density (kg.m^{-3})'),ylabel('Depth (m)')
%     grid on
% 
%     subplot(1,4,2)
%     plot(s1(ind_down),d1(ind_down),'.'),set(gca,'ydir','reverse')
%     xlabel('Salinity (psu)'),ylabel('Depth (m)')
%     grid on
% 
%     subplot(1,4,3)
%     plot(t1(ind_down),d1(ind_down),'.'),set(gca,'ydir','reverse')
%     xlabel('Temperature (^{\circ}C)'),ylabel('Depth (m)')
%     grid on
% 
%     subplot(1,4,4)
%     plot(s1(ind_down),t1(ind_down),'.')%,set(gca,'ydir','reverse')
%     xlabel('Salinity (psu)'),ylabel('Temperature (^{\circ}C)')
%     title('TS')
%     grid on
    
%     figure
%     set(gcf,'color','w');grid on
%     plot(ctd_time,ctd_depth,'.'),set(gca,'ydir','reverse'),datetick('x',1),axis tight
%     ylabel('Depth (m)')
    
% ctd_time_up=ctd_time(ind_up);
% id_up=id(ind_up);
% ctd_depth_up=ctd_depth(ind_up);
% ctd_lat_up=ctd_lat(ind_up);
% ctd_lon_up=ctd_lon(ind_up);

ctd_time=ctd_time(ind_down);
id=id(ind_down);
ctd_depth=ctd_depth(ind_down);
ctd_lat=ctd_lat(ind_down);
ctd_lon=ctd_lon(ind_down);
    
    save '/Users/mpearson/Documents/AZFP MATLAB Toolbox_1/newest_versions/IN_USE/gosl/ctd4azfp.mat'...
        id ctd_time ctd_depth ctd_lat ctd_lon

% ups_downs=glider_profiler(cat(2,filename.(StrName).time,...
%     filename.(StrName).temperature,filename.(StrName).salinity,...
%     filename.(StrName).pressure,filename.(StrName).depth,'plot','yes'));

% ups_downs=glider_profiler(time,t,s,p,d);
