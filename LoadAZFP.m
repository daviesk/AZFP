function [Output,Parameters] = LoadAZFP(varargin)

% LoadAZFP.m run as:
% [Output,Par] = LoadAZFP('datafilename','','xmlfilename','15051410.XML','Salinity',35,'Bins2Avg',5,'Time2Avg',100);
% 
% Inputs can be in any order or omitted and the defaults will be used:
% defaultxmlfilename = ''; % prompt for XML filename
% defaultdatafilename = ''; % AZFP hourly data file name
% defaultSalinity = 35; % Salinity in 
% defaultBins2Avg = 1; % number of range bins to average
% defaultTime2Avg = 1; % number of time values to average
% 
% Outputs are: 
% Output: computed data with N, Sv and TS, averaged in range/time. Each
% freq stored in Output(1), Output(2) etc
% Par: the instrument parameters from the XML file
%
% Ver 1.3 September 2017 
% written by Dave Billenness
% ASL Environmental Sciences Inc.
% 1-6703 Rajpur Place, Victoria, B.C., V8M 1Z5, Canada
% T: +1 (250) 656-0177 ext. 126
% E: dbillenness@aslenv.com 
% w: http://www.aslenv.com/ 
% For any suggestions, comments, questions or collaboration, please contact me.

p = inputParser;
defaultParameters = '';
defaultdatafilename = '';
defaultBins2Avg = 1;
defaultTime2Avg = 1;
defaultPressure = 50;
defaultSalinity = 35;
addParameter(p,'Parameters',defaultParameters,@isstruct);
addParameter(p,'datafilename',defaultdatafilename,@ischar);
addParameter(p,'Bins2Avg',defaultBins2Avg,@isnumeric);
addParameter(p,'Time2Avg',defaultTime2Avg,@isnumeric);
addParameter(p,'Pressure',defaultPressure,@isnumeric);
addParameter(p,'Salinity',defaultSalinity,@isnumeric);
parse(p,varargin{:});
Parameters = p.Results.Parameters;
datafilename = p.Results.datafilename;
Bins2Avg = p.Results.Bins2Avg;
Time2Avg = p.Results.Time2Avg;
Press = p.Results.Pressure;    % pressure (dBar)
Salinity = p.Results.Salinity;

% initialize variables
Data = [];
Output = [];
pathname = '';

% select an AZFP hourly file
if(isempty(datafilename))
    [filename, pathname] = uigetfile({'*.*A;*.*B;*.*C;*.*D'}, 'Select an AZFP hourly file');
else
    pathname = pwd;
    filename = datafilename;
end

% open binary file
fidAZFP = fopen([pathname '/' filename],'rb');
if(fidAZFP < 1)
    [filename, pathname] = uigetfile({'*.*A;*.*B;*.*C;*.*D'}, 'Select an AZFP hourly file');
    fidAZFP = fopen([pathname '/' filename],'rb');
end

% load the data
ii = 1;
while(~feof(fidAZFP))
    [Data,Success] = readAZFP(fidAZFP,ii,Data,Parameters);
    if(Success)
        % preallocate array if data averaging to #values in the hourly file x number of ranges
        if(ii==1 && (Bins2Avg > 1 || Time2Avg > 1))
            for(jj=1:Data(ii).NumChan)
                NumAvg = 1;
                if(Data(1).DataType(jj)) % if averaged data then divide by # pings
                    NumAvg = Data(1).NumAcqPings;
                end
                PavgArr(jj).data = zeros(3600/Data(1).BurstInt*Data(1).PingPerProfile/NumAvg,Data(1).NumBins(jj));
            end
        end
    
        if(ii==1)
            [PATHSTR,NAME,EXT] = fileparts(Parameters.xmlfilename);
            fprintf('File: %s - Loading Profile #%d %s using xml=%s Bins2Avg=%d Time2Avg=%d Salinity=%.2f Pressure=%.1f\n',filename,Data(ii).ProfileNumber,datestr(Data(ii).Date),[NAME EXT],Bins2Avg,Time2Avg,Salinity,Press);
        end
        % compute temperature from Data(ii).Ancillary(5)
        Data(ii).T = NaN;
        %if(Data(ii).SensorFlag)
            Data(ii).T = computeTemperature(Data(ii).Ancillary(5),Parameters);
        %end        
        
        % compute tilts from Data(ii).Ancillary(1) and
        % Data(ii).Ancillary(2), then cos(tiltmag)
        Data(ii).Tx = computeTilt(Data(ii).Ancillary(1),Parameters.X_a,Parameters.X_b,Parameters.X_c,Parameters.X_d);
        Data(ii).Ty = computeTilt(Data(ii).Ancillary(2),Parameters.Y_a,Parameters.Y_b,Parameters.Y_c,Parameters.Y_d);
        Data(ii).cosTiltMag = cos((sqrt(Data(ii).Tx^2+Data(ii).Ty^2))*pi/180);
        
        % compute power if we are averaging the data
        if(Bins2Avg > 1 || Time2Avg > 1)
            for(jj=1:Data(ii).NumChan)
                EL = Parameters.EL(jj) - 2.5/Parameters.DS(jj) + Data(ii).counts{jj}/(26214*Parameters.DS(jj));
                P = 10.^(EL/10);
                if(~isempty(P))
                    Data(ii).Pavg{jj} = P';
                    % create an array so easier to time average
                    PavgArr(jj).data(ii,:) = Data(ii).Pavg{jj}';                
                end
            end
        end
        ii = ii + 1;
    else
        % check if entire file is empty, or just an end block of data,
        % either return [] data or break out of loop
        if(isempty(Data))
            fclose(fidAZFP);
            return;
        else
            break;
        end
    end
end
% close the binary file
fclose(fidAZFP);

% compute hourly average temperature, then use this to compute SoundSpeed
% and the range bins, store in Data(1)
Data(1).HourlyAvgTemp = nanmean(arrayfun(@(x) Data(x).T, 1:length(Data)));
if(isnan(Data(1).HourlyAvgTemp))
    Data(1).HourlyAvgTemp = 6; %to do, move this to the setup parameter file
    fprintf('**** No AZFP temperature found - using a fixed temperature of %.1f degC to calc soundspeed and range\n', Data(1).HourlyAvgTemp);
end
% calculate sound speed using input Salinity
Data(1).SoundSpeed = computeSS(Data(1).HourlyAvgTemp,Press,Salinity);

% compute the hourly avg value of the cos(TiltMag)
Data(1).HourlyAvgcosTiltMag = mean(arrayfun(@(x) Data(x).cosTiltMag, 1:length(Data)));

% compute range and avg range (if depth averaging)
for(jj=1:Data(1).NumChan)
    % calc the number of averaged blocks and the range to the centre of the
    % sampling volume for bin m (from eqn. 11 on page 86 of the AZFP Operators Manual)
    m = 1:length(1:Bins2Avg:length(Data(1).counts{jj})-Bins2Avg+1);
    Data(1).Range{jj} = Data(1).SoundSpeed*Data(1).LockoutInd(jj)/(2*Data(1).DigRate(jj))+(Data(1).SoundSpeed/4)*(((2*m-1)*Data(1).RangeSamples(jj)*Bins2Avg-1)/Data(1).DigRate(jj)+Data(1).PulseLength(jj)/1e6);

    % calc absorption coeff for each freq
    Data(1).seaAbs(jj) = computeAbs(Data(1).HourlyAvgTemp,Press,Salinity,Data(1).Freq(jj));            
    
end

% check averaging time, max is all values in hourly file
if(Time2Avg > length(Data))
    Time2Avg = length(Data);
end

% now bin average, after all data from a single file has been loaded 
if(Bins2Avg > 1)
    for(jj=1:Data(1).NumChan)
        % calc #bins after averaging
        numBins = Data(1).counts{jj};
        numBins = length(arrayfun(@(i) mean(numBins(i:i+Bins2Avg-1)),1:Bins2Avg:length(numBins)-Bins2Avg+1));
        for nn = 0:numBins-1
            tmp1(:,nn+1) = mean(PavgArr(jj).data(:,nn*Bins2Avg+1:(nn+1)*Bins2Avg),2);
        end
        PavgArr(jj).data = tmp1;
        clear tmp1
    end
end

% now time average
if(Time2Avg > 1)
    NumTime = floor(length(Data)/Time2Avg);
    for(kk=1:NumTime)
        % Elements of array to average, Time2Avg = 30 then Elem=1-30, 31-60 etc
        Elem = [(kk-1)*Time2Avg+1:(kk-1)*Time2Avg+Time2Avg]';
        for(jj=1:Data(1).NumChan)
            % convert back to counts N
            ELavg = 10*log10(mean(PavgArr(jj).data(Elem,:),1));
            Output(jj).N(kk,:) = round(26214*Parameters.DS(jj)*(ELavg - Parameters.EL(jj) + 2.5/Parameters.DS(jj)));
            Output(jj).Range = Data(1).Range{jj};
            Output(jj).TiltCorrRange = Data(1).Range{jj}*Data(1).HourlyAvgcosTiltMag;
            % calc correction to Sv due to non square transmit pulse
            SvOffset = CalcSvOffset(Data(1).Freq(jj),Data(1).PulseLength(jj));
            Output(jj).Sv(kk,:) = Parameters.EL(jj)-2.5/Parameters.DS(jj)+Output(jj).N(kk,:)/(26214*Parameters.DS(jj))-Parameters.TVR(jj)-20*log10(Parameters.VTX(jj))+20*log10(Output(jj).Range)+2*Data(1).seaAbs(jj)*Output(jj).Range-10*log10(0.5*Data(1).SoundSpeed*Parameters.PulseLen(jj)/1e6*Parameters.BP(jj)) + SvOffset;
            Output(jj).TS(kk,:) = Parameters.EL(jj)-2.5/Parameters.DS(jj)+Output(jj).N(kk,:)/(26214*Parameters.DS(jj))-Parameters.TVR(jj)-20*log10(Parameters.VTX(jj))+40*log10(Output(jj).Range)+2*Data(1).seaAbs(jj)*Output(jj).Range;
            Output(jj).Freq = Data(1).Freq(jj);
            Output(jj).seaAbs = Data(1).seaAbs(jj);
        end
        Output(1).Date(kk,1) = mean(arrayfun(@(x) Data(x).Date, Elem))';
        Output(1).Tx(kk,1) = mean(arrayfun(@(x) Data(x).Tx, Elem))';
        Output(1).Ty(kk,1) = mean(arrayfun(@(x) Data(x).Ty, Elem))';
        Output(1).T(kk,1) = mean(arrayfun(@(x) Data(x).T, Elem))';
    end
    
    
else % no time averaging, but may still have range averaging
    if(Bins2Avg > 1)
        for(jj=1:Data(1).NumChan)
            % convert back to counts N
            ELavg = 10*log10(PavgArr(jj).data);
            Output(jj).N = round(26214*Parameters.DS(jj)*(ELavg - Parameters.EL(jj) + 2.5/Parameters.DS(jj)));
            Output(jj).Range = Data(1).Range{jj};
            Output(jj).TiltCorrRange = Data(1).Range{jj}*Data(1).HourlyAvgcosTiltMag;
        end
    else
        for(jj=1:Data(1).NumChan)
            Output(jj).N = cell2mat(arrayfun(@(x) Data(x).counts{jj}', [1:length(Data)]', 'UniformOutput', false));
            Output(jj).Range = Data(1).Range{jj};
            Output(jj).TiltCorrRange = Data(1).Range{jj}*Data(1).HourlyAvgcosTiltMag;
        end
    end
    for(jj=1:Data(1).NumChan)
        % calc correction to Sv due to non square transmit pulse
        SvOffset = CalcSvOffset(Data(1).Freq(jj),Data(1).PulseLength(jj));
        Output(jj).Sv = Parameters.EL(jj)-2.5/Parameters.DS(jj)+Output(jj).N./(26214*Parameters.DS(jj))-Parameters.TVR(jj)-20*log10(Parameters.VTX(jj))+20*log10(repmat(Output(jj).Range,[size(Output(jj).N,1) 1]))+2*Data(1).seaAbs(jj).*(repmat(Output(jj).Range,[size(Output(jj).N,1) 1]))-10*log10(0.5*Data(1).SoundSpeed*Parameters.PulseLen(jj)/1e6*Parameters.BP(jj)) + SvOffset;
        Output(jj).TS = Parameters.EL(jj)-2.5/Parameters.DS(jj)+Output(jj).N./(26214*Parameters.DS(jj))-Parameters.TVR(jj)-20*log10(Parameters.VTX(jj))+40*log10(repmat(Output(jj).Range,[size(Output(jj).N,1) 1]))+2*Data(1).seaAbs(jj).*(repmat(Output(jj).Range,[size(Output(jj).N,1) 1]));
        Output(jj).Freq = Data(1).Freq(jj);
        Output(jj).seaAbs = Data(1).seaAbs(jj);
    end
    Output(1).Date = arrayfun(@(x) Data(x).Date, 1:length(Data))';
    Output(1).Tx = arrayfun(@(x) Data(x).Tx, 1:length(Data))';
    Output(1).Ty = arrayfun(@(x) Data(x).Ty, 1:length(Data))';
    Output(1).T = arrayfun(@(x) Data(x).T, 1:length(Data))';
    
end
% trim N, Sv, TS to length of Output(1).Date
% size of Output(1).N if averaging was based on pre-alloc PavgArr array
% before the hourly file was read in (there may be missing time values)
if(size(Output(1).N,1) > size(Output(1).Date,1))
    Output(1).N(size(Output(1).Date,1)+1:end,:)= [];
    Output(1).TS(size(Output(1).Date,1)+1:end,:)= [];
    Output(1).Sv(size(Output(1).Date,1)+1:end,:)= [];
end

% test to extract a range of bins
ExtractBinRange = 0;
startbin = 4150;
endbin = 4350;
if(ExtractBinRange)
    Output(1).N = Output(1).N(:,startbin:endbin);
    Output(2).N = Output(2).N(:,startbin:endbin);
    Output(3).N = Output(3).N(:,startbin:endbin);
    Output(4).N = Output(4).N(:,startbin:endbin);
    Output(1).Sv = Output(1).Sv(:,startbin:endbin);
    Output(2).Sv = Output(2).Sv(:,startbin:endbin);
    Output(3).Sv = Output(3).Sv(:,startbin:endbin);
    Output(4).Sv = Output(4).Sv(:,startbin:endbin);
    Output(1).TS = Output(1).TS(:,startbin:endbin);
    Output(2).TS = Output(2).TS(:,startbin:endbin);
    Output(3).TS = Output(3).TS(:,startbin:endbin);
    Output(4).TS = Output(4).TS(:,startbin:endbin);
    Output(1).Range = Output(1).Range(startbin:endbin);
    Output(2).Range = Output(2).Range(startbin:endbin);
    Output(3).Range = Output(3).Range(startbin:endbin);
    Output(4).Range = Output(4).Range(startbin:endbin);
end


%fprintf('Date = %dx%d N = %dx%d\n',size(Output(1).Date),size(Output(1).N))
Output(1).HourlyAvgTemp = Data(1).HourlyAvgTemp;
Output(1).SoundSpeed = Data(1).SoundSpeed;
Output(1).NumChan = Data(1).NumChan;
Output(1).BurstInt = Data(1).BurstInt;
Output(1).PingPerProfile = Data(1).PingPerProfile;
Output(1).NumAcqPings = Data(1).NumAcqPings;
Output(1).DataType = Data(1).DataType;

end

% **************************************************************
% function readAZFP data block
function [Data,Success] = readAZFP(fidAZFP,ii,Data,Parameters)

% initialize variables
Success = 1;
FileType = 'FD02'; %specify profile data filetype

% read binary data file in big-endian format (ieee-be)
Flag = dec2hex(fread(fidAZFP,1,'uint16','ieee-be'));
if(~strcmpi(Flag,FileType))
    Success = 0;
    if(~feof(fidAZFP))
        fprintf('Error: Unknown file type: check that the correct XML file was loaded\n');
    end
    return;
end
Data(ii).ProfileFlag = Flag;
Data(ii).ProfileNumber = fread(fidAZFP,1,'uint16','ieee-be');
Data(ii).SerialNumber = fread(fidAZFP,1,'uint16','ieee-be');
Data(ii).PingStatus = fread(fidAZFP,1,'uint16','ieee-be');
Data(ii).BurstInt = fread(fidAZFP,1,'uint32','ieee-be');
date = fread(fidAZFP,7,'uint16','ieee-be'); % YY MM DD hh mm ss hh
Data(ii).Date = datenum(date(1),date(2),date(3),date(4),date(5),date(6)+date(7)/100);
Data(ii).DigRate = fread(fidAZFP,4,'uint16','ieee-be'); % digitization rate for each channel
Data(ii).LockoutInd = fread(fidAZFP,4,'uint16','ieee-be'); % lockout index for each channel
Data(ii).NumBins = fread(fidAZFP,4,'uint16','ieee-be'); % number of bins for each channel
Data(ii).RangeSamples = fread(fidAZFP,4,'uint16','ieee-be'); % range samples per bin for each channel
Data(ii).PingPerProfile = fread(fidAZFP,1,'uint16','ieee-be'); % number of pings per profile
Data(ii).AvgPings = fread(fidAZFP,1,'uint16','ieee-be'); % flag to indicate if pings avg in time
Data(ii).NumAcqPings = fread(fidAZFP,1,'uint16','ieee-be'); % # pings acquired in this burst
Data(ii).PingPeriod = fread(fidAZFP,1,'uint16','ieee-be'); % ping period in seconds
Data(ii).FirstLastPing = fread(fidAZFP,2,'uint16','ieee-be');
Data(ii).DataType = fread(fidAZFP,4,'uint8','ieee-be'); % datatype for each channel: 1=Avg Data (5bytes), 0=raw (2bytes)
Data(ii).DataError = fread(fidAZFP,1,'uint16','ieee-be');% error # is an error occurred
Data(ii).Phase = fread(fidAZFP,1,'uint8','ieee-be'); % phase # used to acquire this profile
Data(ii).Overrun = fread(fidAZFP,1,'uint8','ieee-be'); % 1 if an over run occurred
Data(ii).NumChan = fread(fidAZFP,1,'uint8','ieee-be'); % 1,2,3 or 4 (could acquire only 1 channel)
Data(ii).Gain = fread(fidAZFP,4,'uint8','ieee-be'); % gain chan 1-4
fread(fidAZFP,1,'uint8','ieee-be'); %spare chan
Data(ii).PulseLength = fread(fidAZFP,4,'uint16','ieee-be'); % pulselength chan 1-4 uS
Data(ii).BoardNum = fread(fidAZFP,4,'uint16','ieee-be'); % the board the data came from chan 1-4
Data(ii).Freq = fread(fidAZFP,4,'uint16','ieee-be'); % freq Hz for chan 1-4
Data(ii).SensorFlag = fread(fidAZFP,1,'uint16','ieee-be');% Flag to indicate if pressure sensor or temper sensor is avail
Data(ii).Ancillary = fread(fidAZFP,5,'uint16','ieee-be'); % Tilt-X, Y, Battery, Pressure, Temperature
Data(ii).AD = fread(fidAZFP,2,'uint16','ieee-be'); % AD chan 6 and 7
% read in the data, bytes depend on avg or raw, # channels 1 up to 4
for(jj=1:Data(ii).NumChan) % might have 4 freq but only 1 was acquired
    if(Data(ii).DataType(jj)) % averaged data = 32 bit summed up linear scale data followed by 8 bit overflow counts
        if(Data(ii).AvgPings)
            divisor = Data(ii).PingPerProfile * Data(ii).RangeSamples(jj);
        else
            divisor = Data(ii).RangeSamples(jj);
        end
        ls = fread(fidAZFP,Data(ii).NumBins(jj),'uint32','ieee-be'); %linearsum
        lso = fread(fidAZFP,Data(ii).NumBins(jj),'uchar','ieee-be'); %linearsumoverflow
        v = (ls + lso*4294967295)/divisor;
        v = (log10(v)-2.5)*(8*65535)*Parameters.DS(jj);
        v(isinf(v)) = 0;
        Data(ii).counts{jj} = v;
    else % raw data = 16 bit values Log values
        Data(ii).counts{jj} = fread(fidAZFP,Data(ii).NumBins(jj),'uint16','ieee-be');
    end    
end

end

% *************************
function T = computeTemperature(counts,Par)

Vin = 2.5 * (counts / 65535);
R = (Par.ka + Par.kb * Vin) / (Par.kc - Vin);
T = 1 / (Par.A + Par.B * (log(R)) + Par.C * (log(R)^3)) - 273;

end

function Tilt = computeTilt(N,a,b,c,d)

% X_a; X_b; X_c and X_d
Tilt = a + b*(N) + c*(N)^2 + d*(N)^3;

end

function seaC = computeSS(T,P,S)

% from Fundamentals of Physical Acoustics By David T. Blackstock
z = T/10;
seaC = 1449.05+z*(45.7+z*((-5.21)+0.23*z))+(1.333+z*((-0.126)+z*0.009))*(S-35.0)+(P/1000)*(16.3+0.18*(P/1000));

end

% calc Absorption coeff using Temperature, Pressure and Salinity and
% transducer frequency
function seaAbs = computeAbs(T,P,S,Freq)

% calculate relaxation frequencies
T_K = T + 273.0;
f1 = 1320.0*(T_K)*exp(-1700/T_K);
f2 = (1.55e7)*T_K*exp(-3052/T_K);

% coefficients for absorption equations
k = 1 + P/10.0;
a = (8.95e-8)*(1+T*((2.29e-2)-(5.08e-4)*T));
b = (S/35.0)*(4.88e-7)*(1+0.0134*T)*(1-0.00103*k+(3.7e-7)*(k*k));
c = (4.86e-13)*(1+T*((-0.042)+T*((8.53e-4)-T*6.23e-6)))*(1+k*(-(3.84e-4)+k*7.57e-8));
freqk = Freq*1000;
if(S == 0)
    seaAbs = c*(freqk.^2);
else
    seaAbs = (a*f1*(freqk.^2))./((f1*f1)+(freqk.^2))+(b*f2*(freqk.^2))./((f2*f2)+(freqk.^2))+c*(freqk.^2);
end


end

% A correction must be made to compensate for the effects of the
% finite response times of both the receiving and transmitting parts of the instrument. The
% magnitude of the correction will depend on the length of the transmitted pulse, and the response
% time (on both transmission and reception) of the instrument.
function SvOffset = CalcSvOffset(Frequency,PulseLength)

SvOffset = 0;

if(Frequency > 38) % 125,200,455,769 kHz
    if(PulseLength == 300)
        SvOffset = 1.1;
    elseif(PulseLength == 500)
        SvOffset = 0.8;
    elseif(PulseLength == 700)
        SvOffset = 0.5;
    elseif(PulseLength == 900)
        SvOffset = 0.3;
    elseif(PulseLength == 1000)
        SvOffset = 0.3;
    end
else % 38 kHz
    if(PulseLength == 500)
        SvOffset = 1.1;
    elseif(PulseLength == 1000)
        SvOffset = 0.7;
    end
end

end