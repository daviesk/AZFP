% % this routine merges a ctd file with a Sv azfp file (prepared ealier)
% % it then assigns a glider depth to each ping range, 
% % creates dives, that will group all profiles belonging to one dive
% % so that they can be averaged at the same depth during one single profile
% % it then transforms each dive into a profile,
% % plots the whole thing
% 
% addpath(genpath('/Users/mpearson/Documents/AZFP MATLAB Toolbox_1/newest_versions/IN_USE/gosl'))

% don't need to add these separately because it is inside gosl
% addpath(genpath('/Users/mpearson/Documents/AZFP MATLAB Toolbox_1/newest_versions/IN_USE/gosl/m_map'))
% addpath(genpath('/Users/mpearson/Documents/AZFP MATLAB Toolbox_1/newest_versions/IN_USE/gosl/cmocean_v1'))
% % 

test nr 3

tic

% % --------------------
% runs AZFP toolbox and gets some of its Output data, i.e. volume
% backscatter and bin ranges
clear
ParametersAZFP;
[Output,Par]=ProcessAZFP(Parameters); % code can handle 3-4 days of data at a time
% if range averaging, there is an indexing bug on channels 2-4.  This removes
% the bad data (but doesn't fix the bug!)
ii=find(Output(2).N(:,1)<-1000000);
for kk=2:4
    Output(kk).N(ii,:)=[];
    Output(kk).TS(ii,:)=[];
    Output(kk).Sv(ii,:)=[];
end;

azfp_time=Output(1).Date;
echoSv1=Output(1).Sv;
echoSv2=Output(2).Sv;
% echoSv3=Output(3).Sv;
% echoSv4=Output(4).Sv;
echoTS1=Output(1).TS;
echoTS2=Output(2).TS;
range1=Output(1).Range;range1=range1(1,:);
range2=Output(2).Range;range2=range2(1,:);
% range3=Output(3).Range;range3=range3(1,:);
% range4=Output(4).Range;range4=range4(1,:);

% saves time,echo,range for all frequencies
% echoSv3 echoSv4 range3 range4
%save '/Users/mpearson/Documents/AZFP MATLAB Toolbox_1/newest_versions/IN_USE/gosl/echo11.mat' ...
%azfp_time echoSv1 echoSv2 echoTS1 echoTS2 range1 range2

save('echo11.mat','azfp_time','echoSv1','echoSv2','echoTS1','echoTS2','range1','range2')


% 
%cd ../../..
% 
% % --------------------
% works on CTD data through its own routine ctd.m
% ctd

clear all
close all

% --------------------
% loads only what we need
load('ctd4azfp.mat') % created using ctd.m
clear id
load('echo11.mat') % created using dave billenness routines

% value={'Choose frequency (1,2,3 or 4)'};
% value=inputdlg(value);
% value=str2double(value);
% % fprintf('value 1 = %d', value(1))

% add line to ask which frequency we want - how to transform the answer
% into a concatenated string?
echoSv=echoSv1; % or change it to whatever frequency we want 1(125),2(200),3(455),4(769)
bin_depth=range1;

clear echoSv1 echoSv2 range1 range2 echoTS1 echoTS2

% % quickly check inputs
% plot(ctd_time,ctd_depth,'k.'),datetick('x',0),set(gca,'Ydir','reverse'),axis tight
% 
% %figure, imagesc(echo1',[0 65535]) -> if using counts
% figure, imagesc(echoSv')

%% allocating space for the "new format"
glider_depth=nan(size(azfp_time)); % allocates some space for glider depth
lat=nan(size(azfp_time)); % allocates some space for lat & lon
lon=nan(size(azfp_time));
% prof_id=nan(size(azfp_time));

%% assuming AZFP_TIME IS IN matlab time, this bit finds the closest ctd time to azfp times and grabs a glider depth for each time stamp
    for ii=1:length(azfp_time)
%         index=find(ctd_time==azfp_time(ii) | ctd_time<(azfp_time(ii)-6.0069e-05) | ctd_time>(azfp_time(ii)+6.0069e-05)); % assuming 1 ping every 5 seconds , or 6.0069e-05 in Matlab time
        [~,index]=min(abs(ctd_time-azfp_time(ii)));% amongst all ctd times it finds the smallest difference between each azfp and ctd time, grabs its index and then grabs the depth,lat/lon associated to this index
            if ctd_depth(index)<15 % it gets rid off values shallower than 10.2m because glider sometimes gets funky when shallow (it will later screw dive counts if not removing these values)
            glider_depth(ii)=nan;
            lat(ii)=nan;
            lon(ii)=nan;
%             prof_id(ii)=nan;
            else
            glider_depth(ii)=ctd_depth(index);
            lat(ii)=ctd_lat(index);
            lon(ii)=ctd_lon(index);
%             prof_id(ii)=id(index);
            end
    end 

clear ii index
    
ind=find(isnan(glider_depth));
azfp_time(ind)=nan;
echoSv(ind,:)=nan;

glider_depth(ind)=[];
lat(ind)=[];
lon(ind)=[];
azfp_time(ind)=[];
echoSv(ind,:)=[];

% sanity check
figure
pcolor(echoSv'),shading flat
set(gca,'Ydir','reverse')
cmocean('balance') %balance, curl is greenish
colorbar

figure
% subplot(2,1,1)
set(gcf,'color','w');
plot(ctd_time,ctd_depth,'k.'),datetick('x','mmm dd HH PM'),set(gca,'Ydir','reverse'),axis tight %(12:717)
hold on
plot(azfp_time,glider_depth,'rx'),datetick('x','mmm dd HH PM'),set(gca,'Ydir','reverse'),axis tight
hold off

%% find indices location of max counts for bottom detection, start at index 10 so it skips the area near transducer
% to eliminate from plots all stuff below seafloor
[~,loc]=max(echoSv(:,10:end),[],2); % why 10?
ExtraY=9; % extra bins to plot beyond bottom % for echo5.mat -> ExtraY=15, otherwise 50
BottomLoc=loc+ExtraY;
[n,l]=max(BottomLoc); % this is to accomodate all different "depths" of profiles <o> in the original one n was BottomLoc(1)

rot_Sv=nan(size(echoSv,1),n);
rot_depth=nan(size(echoSv,1),n);
rot_time=nan(size(echoSv,1),n);

tmp_depth=glider_depth*ones(1,n)+ones(length(glider_depth),1)*bin_depth(1:n);  % dont understand this line.  BUG HERE!! bin_depth(1:n);
tmp_time=repmat(azfp_time,[1,n]);

% % sanity check
% subplot(2,1,1)
% pcolor(tmp_depth'),shading flat,set(gca,'Ydir','reverse'),colorbar
% subplot(2,1,2)
% pcolor(echoSv'),shading flat,set(gca,'Ydir','reverse'),colorbar

% starts "rotating" to have the floor signal at the bottom
for ii=1:size(echoSv,1) % first time, second time, so on...
   
    if(BottomLoc(l)-BottomLoc(ii)+1<1) % need to keep this bit otherwise it goes to crap...don't know exactly why
        break;
    end
   
    rot_Sv(ii,BottomLoc(l)-BottomLoc(ii)+1:BottomLoc(l))=echoSv(ii,1:BottomLoc(ii));
    rot_depth(ii,BottomLoc(l)-BottomLoc(ii)+1:BottomLoc(l))=tmp_depth(ii,1:BottomLoc(ii));
    rot_time(ii,BottomLoc(l)-BottomLoc(ii)+1:BottomLoc(l))=tmp_time(ii,1:BottomLoc(ii));
%     
end

clear ii

% x=isnan(rot_Sv);
% rot_Sv(x)=100;
% rot_depth(x)=nan;
%[x y]=ind2sub(size(Output.TS),max_idx); % turns the linear position into i,j index position
%[x y]=ind2sub(size(tmp_depth),ind_d);
% rot_Sv(isnan(rot_Sv))=0;
% rot_depth(isnan(rot_Sv))=0;
% rot_time(isnan(rot_Sv))=0;

figure,plot(rot_depth,'.'),set(gca,'Ydir','reverse'),grid on
figure,set(gcf,'color','w'),pcolor(rot_Sv'),shading flat,set(gca,'Ydir','reverse'),cmocean('balance'),colorbar
% contour(rot_Sv'),set(gca,'Ydir','reverse')
% figure,surf(rot_time,rot_depth,rot_Sv)

% figure
% % pcolor(tmp_time',tmp_depth',rot_Sv'),shading flat
% pcolor(rot_Sv'),shading flat
% % figure,imagesc(rot_Sv')
% % figure,pcolor(rot_depth'),shading flat
% set(gca,'Ydir','reverse')
% cmocean('balance') %balance, curl is greenish
% colorbar
% datetick('x','mmm dd HH PM')
% hold on 
% yyaxis right
% plot(azfp_time,glider_depth,'k.')
% set(gca,'Ydir','reverse')
% datetick('x','mmm dd HH PM')

% %sanity checks
% figure,plot(azfp_time,glider_depth,'.')
% set(gca,'Ydir','reverse')
% datetick('x','HH PM')
% 
% [X,Y]=meshgrid(azfp_time,glider_depth);
% pcolor(X,Y,rot_Sv');shading flat
% contour(rot_time',rot_depth',rot_Sv')
% 
% tmp=rot_depth(1,:);
% 
% [X,Y]=meshgrid(tmp_time,tmp_depth);

%% convert lat/long to decimal degrees ....  why?????
% Lat=((abs(lat)-floor(abs(lat)/100)*100)/60+floor(abs(lat)/100)).*lat./abs(lat);
% Lon=((abs(lon)-floor(abs(lon)/100)*100)/60+floor(abs(lon)/100)).*lon./abs(lon);

%% find the ends of the separate dives and assign each ping to their specific profile
% Kim's solution
% identifies the beginning of each downcast by time difference. Kim changed from 2*N to 25*N to account for variability in ping rate that occassionally occurs
pp=diff(datenum(azfp_time))*24*60*60; % number of seconds between each ping
i_end=find(pp>6); % identify the ends of profiles

% % % using logical indexing instead to find the end of each profile
% tmp=~isnan(glider_depth); % finds all indexes where there are depth values
% pp=diff(tmp); % looks for differences between 0 (a nan) and a 1 (a depth value exists) (x2-x1, i.e. 0-1)
% i_end=(find(pp==-1)); % +1 shows where a downcast starts (transition between a nan and a value will be +), and -1 where it ends (x2-x1)

i_end=[i_end; length(azfp_time)]; % this adds the last index to add the last profile in this file
% i_end(find(diff(i_end)==1))=[]; % this will tell us where each profile ends
profile_nr=[];% allocates a structure for each profile ~4 profiles per hour, it will give a profile number to each ping, i.e. group each ping that belongs in a specific profile

for i=1:length(i_end)
    profile_nr=[profile_nr; ones(i_end(i)-length(profile_nr),1)*i]; %it creates a sctruture with profile nr for each dive
end

clear i

%% combine data from each dive into a profile with .5 m bins
% from glider you need: echo.i_depth=depth, echo.lat, echo.lon each per ping number
% from ProcessAZFP you need:
% bin_depth=range
% Output.Sv and Output.N  -> Rename them to echo.Sv and echo.peak2peak
% Timestamp (per ping)
% echo.i_depth=glider_depth;
% bin_depth=range;
% echo.timestamp=azfp_time;

% clear glider_depth range azfp_time

nr_dives=max(profile_nr);
% dives_depth=floor((min(glider_depth)+min(bin_depth))*2)/2:.3:ceil((max(glider_depth)+max(bin_depth))*2)/2;
%dives_depth=(min(floor((min(rot_depth,[],2))))):.5:(max(ceil(max(rot_depth,[],2))));
dives_depth=0:0.5:200;

% dives_depth=floor((min(glider_depth)+min(bin_depth))*2)/2:.5:ceil((max(glider_depth)+max(bin_depth))*2)/2;
% this makes up a "fake idealized" equally spaced grid that will be use
% to assign echos within a given range; still to find out what the 2 does

% allocating some space...
% dives_lat2=nan(nr_dives,1); % why not just like this????
dives_lat=nan*ones(nr_dives,1);
dives_lon=nan*ones(nr_dives,1);
dives_time=nan*ones(nr_dives,1);
dives_Sv=nan*ones(nr_dives,length(dives_depth));
% dives_nr_averaged=nan*ones(nr_dives,length(dives_depth));

% gets a mean for the entire dive
for i=1:nr_dives %for each dive
    ind_prof=find(profile_nr==i); % find indices that belong to each profile
    dives_time(i)=nanmean(datenum(azfp_time(ind_prof))); % gets a mean time for each dive; can be done without datenum
    dives_lat(i)=nanmean(lat(ind_prof)); %gets the mean lat/lon for each dive
    dives_lon(i)=nanmean(lon(ind_prof));
    tmp_Sv=rot_Sv(ind_prof,:);
    tmp_Sv1=echoSv(ind_prof,:); % temporary; gets each echo value for each ping (and all ranges) belonging to each profile/dive (1,2,...)
    % this bit assigns a "true" depth for the bin, taking into consideration to the glider_depth
%     n=size(echoSv,2); %it gets the size of range dimension of azfp bins; it is the same as [~,n]=size(echoSv1)
%     original: temp_dep=echo.i_depth(ind_prof)*ones(1,length(ping{1,1}(nf_rang+1:param.rangebn-1)))+ones(length(ind_prof),1)*bin_depth;
%     tmp_depth=glider_depth(ind_prof)*ones(1,n)+ones(length(ind_prof),1)*bin_depth; %it makes the "real" echo depth matrix by adding the glider depth to the bin range from azfp to get the real depth of each bin
    tmp_depth=rot_depth(ind_prof,:);
    for j=1:length(dives_depth) % the "idealized" depth matrix that will include all ping values within a given range
%         ind_d=find(rot_depth(ind_prof,:)>dives_depth(j)-.25&rot_depth(:)<=dives_depth(j)+.25);
        ind_d=find(tmp_depth(:,:)>dives_depth(j)-.25&tmp_depth(:,:)<=dives_depth(j)+.25); %looks up on the real depth matrix for values within 0.5m from the specific idelized bin
%         dives_Sv(i,j)=10*log10(nanmean(10.^(tmp_Sv(ind_d)/10))); % WHY?????? % takes  the MEAN of ALL VALUES IN THAT DEPTH RANGE IN LINEAR SPACE
        dives_Sv(i,j)=nanmean(tmp_Sv(ind_d));
%         nr_averaged_cells(i,j)=length(find(~isnan(tmp_Sv(ind_d)))); % find number/quantity of averaged cells 
    end
end

% [m, n]=size(dives_Sv); % m is total nr. of dives, n is size of the idealized depth matrix
% dives.Sv=nan(m,600); % 600 is an arbitrary value to fill in with info to accomodate max idealized depths we could possibly expect
% dives.depth=nan(1,600);
% dives.no=nan(m,600);

% %just adding to an organized structure
% dives.time=dives_matlabtime;
% dives.lat=dives_lat;
% dives.lon=dives_lon;
% dives.Sv(1:m,1:n)=dives_Sv;
% dives.depth(1,1:n)=dives_depth;
% dives.no(1:m,1:n)=dives_no_averaged;

%% plots

% % Sv - volume backscatter
% figure
% % subplot(2,1,2)
% set(gcf,'color','w');
% % pcolor(dives_Sv'),shading flat,cmocean('balance')
% pcolor(dives_time,dives_depth,dives_Sv'),shading flat
% set(gca,'Ydir','reverse')
% colorbar
% % caxis([-85 -45])
% % datetick('x',21) %'mmm.dd,yyyy HH:MM:SS'
% datetick('x','mmm dd HH PM') 
% % xtickangle(45)
% ylabel('Water depth (m)')
% xlabel('Time (days)')
% axis tight
% cmocean('balance') % balance, curl is greenish

% figure,contourf(dives_time,dives_depth,dives_Sv'),set(gca,'Ydir','reverse')

figure
set(gcf,'color','w');
h1=imagesc(dives_time,dives_depth,dives_Sv')
set(h1,'alphadata',~isnan(dives_Sv'))
datetick('x','mmm dd HH PM','keepticks');
xtickangle(45)%,'Rotation',45.0) 
ylabel('Water depth (m)','FontSize',11)
xlabel('Time UTC (Atlantic time -3 UTC)','FontSize',11)
axis tight
cmocean('balance')
h2=colorbar; %('AxisLocation','in');
ylabel(h2,'Volume backscatter Sv (dB re m^{-1})','FontSize',11);%,'Rotation',270.0)
set(gca,'YLim',[15 90])
caxis([-90 -50])

% glider track
figure
% subplot(3,1,1)
set(gcf,'color','w');
m_proj('albers equal-area','lat',[47.25 48],'long',[-65 -63],'rect','off');
[C,h]=m_etopo2('contour',[-25 -50 -75 -100 -125],'edgecolor',[0.8 0.8 0.8]) 
clabel(C,h)
% m_elev('contour',[-80:1:-5],'edgecolor','b'); % add some bathymetry when possible
h3=m_line(ctd_lon,ctd_lat,'marker','o','color',[0 0 0],...
          'linest','none','markerfacecolor','w','clip','point'); 
h4=m_line(dives_lon,dives_lat,'marker','o','color',[0 .5 0],...
          'linest','none','markerfacecolor','w','clip','point');
xlim([-65 -63])
ylim([47.25 48])
m_gshhs_h('patch',[0 0 0]);
m_grid('linest','none','linewidth',1,'tickdir','in');
xlabel('Longitude','FontSize',11)
ylabel('Latitude','FontSize',11)
title ('Glider track GOSL Jul-Sep 2018','FontSize',11)

toc

save '/Users/mpearson/Documents/AZFP MATLAB Toolbox_1/newest_versions/IN_USE/gosl/dives_echo125/dives_echo125_11.mat' ...
dives_time dives_depth dives_Sv dives_lon dives_lat


%% add

% save stuff

% --------------------
% logical(diff(x)) % to find changes in the profile_ids
 