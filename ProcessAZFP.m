function [Output,Par] = ProcessAZFP(Parameters)

% KDavies 5 Dec 2018.  Line 145-146 for loop will 'break' if there are
% empty data files (.01A, .01B, .01C, .01D).  These must be removed from
% directory for this code to function.

% ProcessAZFP.m run as:
% [Output,Par] = ProcessAZFP(Parameters);
% 
% Inputs can be in any order or omitted and the defaults will be used:
% ProcDir = 0; % 1 will prompt for an entire directory to process
%   = 0 will prompt to load individual files in a directory
% datafilename = ''; % prompt for hourly AZFP file(s) to load
% xmlfilename = ''; % prompt for XML filename if no XML file exists in the directory
% Salinity = 35; % Salinity in 
% Bins2Avg = 10; % number of range bins to average
% Time2Avg = 60; % number of time values to average
% Pressure = 50; %in dbars (~ depth in meters)
% Plot = 0; % show an echogram plot for each channel
% 
% Outputs are: 
% Output: structured array with computed data with N, Sv and TS, averaged in range/time. Each
% freq stored in Output(1), Output(2) etc
% Parameters: the instrument parameters from the XML file
%
% Ver 1.3 September 2017 
% written by Dave Billenness
% ASL Environmental Sciences Inc.
% 1-6703 Rajpur Place, Victoria, B.C., V8M 1Z5, Canada
% T: +1 (250) 656-0177 ext. 126
% E: dbillenness@aslenv.com 
% w: http://www.aslenv.com/ 
% For any suggestions, comments, questions or collaboration, please contact me.

if(~nargin)
    error('Must pass in a list of Parameters => [Output,Par] = ProcessAZFP(Parameters)');
end

% set up input parser
p = inputParser;
% set up defaults
defaultProcDir = 0;
defaultxmlfilename = '';
defaultdatafilename = '';
defaultBins2Avg = 10;
defaultTime2Avg = 60;
defaultPressure = 50;
defaultSalinity = 35;
defaultPlot = 0;
defaultChannel = 1;
defaultValue2Plot = 2;
defaultNoiseFloor = 10000;
defaultOrientation = 1; 
defaultUseTiltCorr = 0;
addParameter(p,'ProcDir',defaultProcDir,@isnumeric);
addParameter(p,'xmlfilename',defaultxmlfilename,@ischar);
addParameter(p,'datafilename',defaultdatafilename,@ischar);
addParameter(p,'Bins2Avg',defaultBins2Avg,@isnumeric);
addParameter(p,'Time2Avg',defaultTime2Avg,@isnumeric);
addParameter(p,'Pressure',defaultPressure,@isnumeric);
addParameter(p,'Salinity',defaultSalinity,@isnumeric);
addParameter(p,'Plot',defaultPlot,@isnumeric);
addParameter(p,'Channel',defaultChannel,@isnumeric);
addParameter(p,'Value2Plot',defaultValue2Plot,@isnumeric);
addParameter(p,'NoiseFloor',defaultNoiseFloor,@isnumeric);
addParameter(p,'Orientation',defaultOrientation,@isnumeric);
addParameter(p,'UseTiltCorr',defaultUseTiltCorr,@isnumeric);
% parse input from Parameter file
parse(p,Parameters);
ProcDir = p.Results.ProcDir;
xmlfilename = p.Results.xmlfilename;
datafilename = p.Results.datafilename;
Bins2Avg = p.Results.Bins2Avg;
Time2Avg = p.Results.Time2Avg;
Pressure = p.Results.Pressure;
Salinity = p.Results.Salinity;
Plot = p.Results.Plot;
Channel = p.Results.Channel;
Value2Plot = p.Results.Value2Plot;
NoiseFloor = p.Results.NoiseFloor;
Orientation = p.Results.Orientation;    
UseTiltCorr = p.Results.UseTiltCorr; 

if(ProcDir)
    % select directory containing the hourly AZFP files to process
    dirname = uigetdir('', 'Select AZFP directory');
    cd(dirname);
    % get a list of all of the AZFP files
    filelist = dir('*.0*'); %the extra dot is to ignore .evi Echoview files, this does not work on MAC OS
    numfiles = length(filelist);
    % better? handling of multi phased data using time stamp 01A 01B files, 
    % ***assumes filelist.date hasn't been modified
    [~,I] = sort(datenum({filelist.date}));
    f = char({filelist.name});
    filelist = f(I,:);
else
    if(isempty(datafilename)) %if no datafilename input, then prompt
        [filelist, dirname] = uigetfile('*.*A;*.*B;*.*C;*.*D', 'Select AZFP hourly file(s)','MultiSelect', 'on');
        if(iscell(filelist))% multiple files selected
            numfiles = length(filelist);
        else %one file selected
            numfiles = size(filelist,1);
        end
        cd(dirname);
    else
        dirname = pwd;
        numfiles = 1;
        filelist = datafilename;
    end
end
pathname = pwd;
if(isempty(xmlfilename)) % if a single xml file is in the directory then load it, otherwise prompt
    xmlfile = dir('*.xml');
    if(length(xmlfile) == 1)
        xmlfilename = char(xmlfile.name);
    end
end
if(isempty(xmlfilename))
    % select an xml settings file
    [xmlfilename, pathname] = uigetfile('*.XML', 'Select instrument coefficients file');
end
if(contains(pathname,'\'))
    pathname = '';
end
Parameters = LoadAZFPxml(pathname,xmlfilename,Parameters);
Parameters.xmlfilename = xmlfilename;

Output(1).Date = [];
Output(1).Tx = [];
Output(1).Ty = [];
Output(1).T = [];
Output(1).filename = {''};
Output(1).HourlyAvgTemp = [];
Output(1).SoundSpeed = [];
for(ii=1:numfiles)
    if(ProcDir) % if proc an entire directory then the file list is in a structure
        fname = filelist(ii,:);
    else
        if(iscell(filelist))% multiple files selected
            fname = char(filelist(ii));
        else %one file selected
            fname = char(filelist(ii,:));
        end
    end
    [DataOut,Par]=LoadAZFP('Salinity',Salinity,'Bins2Avg',Bins2Avg,'Time2Avg',Time2Avg,'datafilename',fname,'Parameters',Parameters,'Pressure',Pressure);
    % check for an empty file and break
    if(isempty(DataOut))
        break;
    end
    for(jj=1:DataOut(1).NumChan)
        if(ii==1)
            Output(jj).N = [];
            Output(jj).Range = [];
            Output(jj).TiltCorrRange = [];
            Output(jj).Sv = [];
            Output(jj).TS = [];   
            Output(jj).seaAbs = [];
        else
            % if DataOut(jj).N has a different number of columns compared to
            % the previous file (stored in Output(jj).N) then catch this and
            % return an error
            if(size(DataOut(jj).N,2) ~= size(Output(jj).N,2))
                error('For a given freq - all files must have the same number of range bins');
            end
        end
        Output(jj).N(end+1:end+1+size(DataOut(jj).N,1)-1,:) = DataOut(jj).N;
        Output(jj).Range(end+1:end+1+size(DataOut(jj).Range,1)-1,:) = DataOut(jj).Range;
        Output(jj).TiltCorrRange(end+1:end+1+size(DataOut(jj).TiltCorrRange,1)-1,:) = DataOut(jj).TiltCorrRange;
        % check for bad tilt values
        if(jj==1 && UseTiltCorr)
            if(mean(DataOut(jj).TiltCorrRange./DataOut(jj).Range) < 0.9) % this would be about 20 deg Tx and Ty tilts
                fprintf('** Warning: tilt correction is set to ON but the tilts are large ... check\n');
            end
        end
        Output(jj).Sv(end+1:end+1+size(DataOut(jj).Sv,1)-1,:) = DataOut(jj).Sv;
        Output(jj).TS(end+1:end+1+size(DataOut(jj).TS,1)-1,:) = DataOut(jj).TS;
        Output(jj).Freq = DataOut(jj).Freq;
        Output(jj).seaAbs(end+1:end+1+size(DataOut(jj).seaAbs,1)-1,:) = DataOut(jj).seaAbs;
    end
    Output(1).Date(end+1:end+1+size(DataOut(1).Date,1)-1,:) = DataOut(1).Date;
    Output(1).Tx(end+1:end+1+size(DataOut(1).Tx,1)-1,:) = DataOut(1).Tx;
    Output(1).Ty(end+1:end+1+size(DataOut(1).Ty,1)-1,:) = DataOut(1).Ty;
    Output(1).T(end+1:end+1+size(DataOut(1).T,1)-1,:) = DataOut(1).T;
    Output(1).filename(ii) = {fname};
    Output(1).HourlyAvgTemp(end+1:end+1+size(DataOut(1).HourlyAvgTemp,1)-1,:) = DataOut(1).HourlyAvgTemp;
    Output(1).SoundSpeed(end+1:end+1+size(DataOut(1).SoundSpeed,1)-1,:) = DataOut(1).SoundSpeed;
end

% save the avg to the Output variable
Output(1).Bins2Avg = Bins2Avg;
Output(1).Time2Avg = Time2Avg;
Output(1).BurstInt = DataOut(1).BurstInt;
Output(1).PingPerProfile = DataOut(1).PingPerProfile;
Output(1).NumAcqPings = DataOut(1).NumAcqPings;
Output(1).DataType = DataOut(1).DataType;

% Plot results, all channels. Plot just Sv, default NoiseFloor and
% Orientation. To plot Counts or TS, and to change orientation see PlotAZFP.m
if Plot && ~isempty(Output(1).Date)
    % plot N, Sv, or Ts
    PlotAZFP(Output,Parameters);       
end

%**********************************************************
function Parameters = LoadAZFPxml(pathname,xmlfilename,Parameters)

if(~isempty(pathname))
    xDoc = xmlread([pathname '/' xmlfilename]);
else
    xDoc = xmlread(xmlfilename);
end

Parameters.NumFreq = str2double(xDoc.getElementsByTagName('NumFreq').item(0).getFirstChild.getData);
Parameters.SerialNumber = str2double(xDoc.getElementsByTagName('SerialNumber').item(0).getFirstChild.getData);
Parameters.BurstInterval = str2double(xDoc.getElementsByTagName('BurstInterval').item(0).getFirstChild.getData);
Parameters.PingsPerBurst = str2double(xDoc.getElementsByTagName('PingsPerBurst').item(0).getFirstChild.getData);
Parameters.AverageBurstPings = str2double(xDoc.getElementsByTagName('AverageBurstPings').item(0).getFirstChild.getData);

% temperature coeff
Parameters.ka = str2double(xDoc.getElementsByTagName('ka').item(0).getFirstChild.getData);
Parameters.kb = str2double(xDoc.getElementsByTagName('kb').item(0).getFirstChild.getData);
Parameters.kc = str2double(xDoc.getElementsByTagName('kc').item(0).getFirstChild.getData);
Parameters.A = str2double(xDoc.getElementsByTagName('A').item(0).getFirstChild.getData);
Parameters.B = str2double(xDoc.getElementsByTagName('B').item(0).getFirstChild.getData);
Parameters.C = str2double(xDoc.getElementsByTagName('C').item(0).getFirstChild.getData);

% tilts
Parameters.X_a = str2double(xDoc.getElementsByTagName('X_a').item(0).getFirstChild.getData);
Parameters.X_b = str2double(xDoc.getElementsByTagName('X_b').item(0).getFirstChild.getData);
Parameters.X_c = str2double(xDoc.getElementsByTagName('X_c').item(0).getFirstChild.getData);
Parameters.X_d = str2double(xDoc.getElementsByTagName('X_d').item(0).getFirstChild.getData);
Parameters.Y_a = str2double(xDoc.getElementsByTagName('Y_a').item(0).getFirstChild.getData);
Parameters.Y_b = str2double(xDoc.getElementsByTagName('Y_b').item(0).getFirstChild.getData);
Parameters.Y_c = str2double(xDoc.getElementsByTagName('Y_c').item(0).getFirstChild.getData);
Parameters.Y_d = str2double(xDoc.getElementsByTagName('Y_d').item(0).getFirstChild.getData);

% get parameters for each transducer freq
for(jj=1:Parameters.NumFreq)
    Parameters.DigRate(jj) = str2double(xDoc.getElementsByTagName('DigRate').item(jj-1).getFirstChild.getData);
    Parameters.LockOutIndex(jj) = str2double(xDoc.getElementsByTagName('LockOutIndex').item(jj-1).getFirstChild.getData);
    Parameters.Gain(jj) = str2double(xDoc.getElementsByTagName('Gain').item(jj-1).getFirstChild.getData);
    Parameters.PulseLen(jj) = str2double(xDoc.getElementsByTagName('PulseLen').item(jj-1).getFirstChild.getData);
    Parameters.DS(jj) = str2double(xDoc.getElementsByTagName('DS').item(jj-1).getFirstChild.getData);
    Parameters.EL(jj) = str2double(xDoc.getElementsByTagName('EL').item(jj-1).getFirstChild.getData);
    Parameters.TVR(jj) = str2double(xDoc.getElementsByTagName('TVR').item(jj-1).getFirstChild.getData);
    Parameters.VTX(jj) = str2double(xDoc.getElementsByTagName('VTX0').item(jj-1).getFirstChild.getData);
    Parameters.BP(jj) = str2double(xDoc.getElementsByTagName('BP').item(jj-1).getFirstChild.getData);
end

Parameters.SensorsFlag = str2double(xDoc.getElementsByTagName('SensorsFlag').item(0).getFirstChild.getData);
