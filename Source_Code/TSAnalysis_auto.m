%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Semi-Automated Surface Wave Two-Station Phase Velocity Dipersion Analysis
% by Huajian Yao, Jan 2007 at MIT
% Note:
% Please do not disclose this program to others! For reference, please use:
% Yao H., van der Hilst R.D., and de Hoop, M.V..Surface-wave array tomography 
% in SE Tibet from ambient seismic noise and two-station analysis : I -
% Phase velocity maps. 2006, Geophysical Journal  International, Vol. 166, 732-744, 
% doi: 10.1111/j.1365-246X.2006.03028.x.
% 
% Last modified date: Sep. 17, 2007
% Last modified data: Dec. 19, 2012 
%      for improving the calculation of group (energy) travel time using
%      Gaussian filtering
% Last modified data: Dec. 31, 2012 
%      fixing a bug when resampling data, improve the codes for instrument
%      response removal
% Last modified data: Sep. 15, 2014
%      adding the function to extract the dispesion curve automatically 
%      by: Chuanming Liu
% Last modified data: Mar. 17, 2015
%      applying the reference disperison curve between the station
%      by: Chuanming Liu

%%%------------------------------------------------------------------------
%%%                    Initialization Module  
%%%------------------------------------------------------------------------
% function varargout = mygui(varargin)  
function varargout = TSAnalysis_auto(varargin)
% TSANALYSIS_AUTO M-file for TSAnalysis_auto.fig
%      TSANALYSIS_AUTO, by itself, creates a new TSANALYSIS_AUTO or raises the existing
%      singleton*.
%
%      H = TSANALYSIS_AUTO returns the handle to a new TSANALYSIS_AUTO or the handle to
%      the existing singleton*.
%
%      TSANALYSIS_AUTO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TSANALYSIS_AUTO.M with the given input arguments.
%
%      TSANALYSIS_AUTO('Property','Value',...) creates a new TSANALYSIS_AUTO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TSAnalysis_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TSAnalysis_auto_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help TSAnalysis_auto

% Last Modified by GUIDE v2.5 07-Apr-2016 15:13:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TSAnalysis_auto_OpeningFcn, ...
                   'gui_OutputFcn',  @TSAnalysis_auto_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before TSAnalysis_auto is made visible.
function TSAnalysis_auto_OpeningFcn(hObject, eventdata, handles, varargin)

% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TSAnalysis_auto (see VARARGIN)

% Choose default command line output for TSAnalysis_auto
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% DefaultPath(handles);

% --- Outputs from this function are returned to the command line.
function varargout = TSAnalysis_auto_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes during object creation, after setting all properties.
function TSFigure_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TSFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% call data struct creation function
DataStructCreate;

% --- Creat structrual body
function DataStructCreate
global StationInfo SourceInfo WaveformInfo RecordInfo CrossCorrWaveInfo
global sta src wave rcd_Z filter cross group
global DataDirectory RespDirectory DispDirectory    
global MonthDays PrevStaLonLat TBeta


DataDirectory = pwd;
RespDirectory = pwd;
DispDirectory = pwd;
MonthDays = [31 29 31 30 31 30 31 31 30 31 30 31; 31 28 31 30 31 30 31 31 30 31 30 31];
PrevStaLonLat = [0 0;0 0];
TBeta = [1 10 20 40 80 100 150 200 300;15 12 10 8 6 5 4 3 2];

StationInfo = struct('Lon',0,...
	'Lat',0,...
	'GCDkm',0,...
	'GCDdeg',0,...
	'Azim',0,...
    'SampleT',0,...
	'SampleF',0,...
    'Name','');

RecordInfo = struct('YY',0,...
	'DD',0,...
	'HH',0,...
	'MM',0,...
	'SS',0,...
	'MS',0,...
	'DiffT',0,...
    'SampleT',0,...
	'SampleF',0,...
    'NumPt',0,...
    'Time',0);
    
WaveformInfo = struct('DatZ',zeros(1,100),...
    'AmpZ',0);

SourceInfo = struct('Lon',0, ...
	'Lat',0, ...
    'YY',0,...
    'Month',0,...
    'Day',0,...
    'DD',0,...
	'HH',0,...
	'MM',0,...
	'SS',0,...
	'MS',0);

RespInfo = struct('Amp',0,...
    'NumZeros',0,...
    'NumPoles',0,...
    'Zeros',0,...
    'Poles',0,...
    'poly_num',0,...
    'poly_den',0);

FilterInfo = struct('Mode',0,...
    'Domain',0,...
    'Window',0,...
    'CtrT',0,...
    'LowT',0,...
    'HighT',0,...
    'CtrF',0,...
    'LowF',0,...
    'HighF',0,...
    'SampleF',0,...
    'SampleT',0,...
    'Length',0,...
    'BandWidth',0,...
    'KaiserPara',1,...
    'GaussAlfa',2.5);

% Period range and interval for narrow band pass filtered cross correlation
CrossCorrWaveInfo = struct('StartT',0,...
    'EndT',0,...
    'StartF',0,...
    'EndF',0,...
    'DeltaT',0,...
    'DeltaF',0,...
    'NumCtrT',0,...
    'StartNum',0,...
    'EndNum',0,...
    'SENum',0,...
    'PointNum',0,...
    'WaveType',0,...
    'WinCode',0,...
    'HalfBand',0,...
    'wave1',zeros(1,100),...
    'wave2',zeros(1,100),...
    'group1',zeros(100,1),...
    'group2',zeros(100,1),...
    'VMin',0,...
    'VMax',0,...
    'DeltaV',0,...
    'PhasVImg',zeros(100,1));

GroupInfo = struct('ImgStartPt','0',...
    'ImgEndPt','0',...
    'ImgType','0',...
    'ArrPt1',zeros(1,100),...
    'ArrPt2',zeros(1,100),...
    'Velo',zeros(1,100));

sta = struct(StationInfo);
rcd_Z = struct(RecordInfo);
src = struct(SourceInfo);
resp = struct(RespInfo);
wave = struct(WaveformInfo);
cross = struct(CrossCorrWaveInfo);
filter = struct(FilterInfo);
group = struct(GroupInfo);

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function  TSFigure_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to TSFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ClickPoint
ClickPoint = get(handles.axes2, 'CurrentPoint');

% --- Executes during object deletion, before destroying properties.
function TSFigure_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to TSFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clear all

% --- Executes when TSFigure is resized.
function TSFigure_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to TSFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%------------------------------------------------------------------------
%%%                 TSAnalysis Setting Module  
%%%------------------------------------------------------------------------ 
% 1.Wave Tpye
% --- Executes on button press in Rayleighwave.
function Rayleighwave_Callback(hObject, eventdata, handles)
% hObject    handle to Rayleighwave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.Rayleighwave,'Value')==1
    set(handles.Lovewave,'Value',0);
end
% Hint: get(hObject,'Value') returns toggle state of Rayleighwave

% --- Executes on button press in Lovewave.
function Lovewave_Callback(hObject, eventdata, handles)
% hObject    handle to Lovewave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.Lovewave,'Value')==1
    set(handles.Rayleighwave,'Value',0);
end
% Hint: get(hObject,'Value') returns toggle state of Lovewave

% --- Executes on button press in Resp_switch.
function Resp_switch_Callback(hObject, eventdata, handles)
% hObject    handle to Resp_switch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.Resp_switch,'Value')==1
   set(handles.MsgEdit,'String','Choose to remove intrusment response, please press the ResFile List button!');
else
   set(handles.MsgEdit,'String','Choose to do not remove intrusment response!');
end


% Hint: get(hObject,'Value') returns toggle state of Resp_switch
% 3.Response Opition
% --- Executes on button press in RespFileInput.
% Open file containing all response data file names  
function RespFileInput_Callback(hObject, eventdata, handles)
% hObject    handle to RespFileInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
%{handles    structure with handles and user data (see GUIDATA)
global RespDirectory allresp

if get(handles.Resp_switch,'Value')==0
   words= 'Please press the Remove Resp botton first!';
   errordlg(words,'Warning'); 
elseif  get(handles.Resp_switch,'Value')==1
   allresp = struct('fname', '', 'staname', '');
   if RespDirectory == 0
     RespDirectory = pwd;
   end
   [pfile1, pname1, index] = uigetfile({'*.txt';'*.dat';'*.*'},'Open file containing all Response file names', RespDirectory);
   files = fullfile(pname1,pfile1);
   RespDirectory = pname1;
   if pname1==0
      set(handles.MsgEdit,'String','File containing all Response file read unsuccseefully !');
   else
      frf = fopen(files);
      if isempty(files)
         set(handles.MsgEdit,'String','Read Stations Response Information unsuccseefully !');
      else
         i = 1;
         while ~feof(frf)
    name = fscanf(frf, '%s', 1);
    if ~strcmp(name,'') 
        allresp(i).fname = name;
%         allresp(i).staname = allresp(i).fname(9:end-7);    
        tmpchar = allresp(i).fname(1:end);
        index = 1;
        Pt = zeros(1,8);
        for j = 1:length(tmpchar)
            if strcmp(tmpchar(j),'.')==1
                Pt(index) = j;
                index = index+1;
            end           
        end
        allresp(i).staname = tmpchar((Pt(2)+1):(Pt(3)-1));
        i = i+1;
    else
        break
    end
         end
         respnum = length(allresp);
         fclose(frf);
         set(handles.MsgEdit,'String','Read Stations Response Information Successfully!');
      end
   end   
end

%%%------------------------------------------------------------------------
%%%                 Load Station Waveform File Module  
%%%------------------------------------------------------------------------  
% --- Read data of two stations (main)
function RdTwoStaData(seisfile1, seisfile2)
global sta wave rcd_Z src source
% if Index_SAC_ASC == 1
%     % for SAC format
%     S = readsac(seisfile1);
%     [sta(1), rcd_Z(1), wave(1).DatZ, wave(1).AmpZ, source] = DataStructTrans(S);
%     S = readsac(seisfile2);
%     [sta(2), rcd_Z(2), wave(2).DatZ, wave(2).AmpZ, source] = DataStructTrans(S);
%    
% elseif Index_SAC_ASC == 2
%     % for SAC_ASC format
%     [sta(1), rcd_Z(1), wave(1).DatZ, wave(1).AmpZ, source] = Rd_Sac_ASC(seisfile1);
%     [sta(2), rcd_Z(2), wave(2).DatZ, wave(2).AmpZ, source] = Rd_Sac_ASC(seisfile2);
%     %[station, record, ReSampleWave, ampcoef, source] = Rd_Sac_ASC(seisfile)
% end

 % for SAC format
S = readsac(seisfile1);
[sta(1), rcd_Z(1), wave(1).DatZ, wave(1).AmpZ, source] = DataStructTrans(S);
S = readsac(seisfile2);
[sta(2), rcd_Z(2), wave(2).DatZ, wave(2).AmpZ, source] = DataStructTrans(S);

for i = 1:2
   % sta(i)
    if (sta(i).GCDkm == -12345) || (sta(i).Azim == -12345) || isnan(sta(i).GCDkm ) || isnan(sta(i).Azim)
        sta(i).GCDkm = deg2km(distance([src.Lat,src.Lon],[sta(i).Lat,sta(i).Lon]));
        sta(i).Azim = azimuth([src.Lat,src.Lon],[sta(i).Lat,sta(i).Lon]);
    end
end

for i = 1:2
    if mod(rcd_Z(i).YY, 4) == 0  
        YearDays = 366;
    else
        YearDays = 365;  
    end            
    DeltaTSrcSta = 0;
    if src.YY == rcd_Z(i).YY 
        DeltaTSrcSta = (src.DD - rcd_Z(i).DD)*24*3600;
    elseif src.YY == rcd_Z(i).YY + 1  
        DeltaTSrcSta = (src.DD + YearDays - rcd_Z(i).DD)*24*3600;
    else
        display('Year Error!');
        break
    end
    DeltaTSrcSta = DeltaTSrcSta + (src.HH - rcd_Z(i).HH)*3600 + (src.MM - rcd_Z(i).MM)*60; 
    DeltaTSrcSta = DeltaTSrcSta + (src.SS - rcd_Z(i).SS) + (src.MS - rcd_Z(i).MS)/1000;
    rcd_Z(i).DiffT = DeltaTSrcSta;
end

% --- Transform Date structure
function [station, record, ReSampleWave, ampcoef, source] = DataStructTrans(HdrData)
% HdrData contains information about the earthquake and receiver
global StationInfo SourceInfo RecordInfo
station = struct(StationInfo);
source = struct(SourceInfo);
record = struct(RecordInfo);
station.Lat = HdrData.STLA;
station.Lon = HdrData.STLO;
station.GCDkm = HdrData.DIST;
station.GCDdeg = HdrData.GCARC;
station.Azim = HdrData.AZ;
station.SampleT = HdrData.DELTA;
station.SampleF = 1/station.SampleT;
station.Name = HdrData.KSTNM;

record.YY = HdrData.NZYEAR;
record.DD = HdrData.NZJDAY;
record.HH = HdrData.NZHOUR;
record.MM = HdrData.NZMIN;
record.SS = HdrData.NZSEC;
record.MS = HdrData.NZMSEC;
record.DiffT = HdrData.O;
record.NumPt = HdrData.NPTS;
record.SampleT = station.SampleT;
record.SampleF = station.SampleF;

source.Lat = HdrData.EVLA;
source.Lon = HdrData.EVLO;

% offset = mean(SeisData(1:record.NumPt));
SeisData(1:record.NumPt) = detrend(HdrData.DATA1(1:record.NumPt));
record.Time = (record.NumPt - 1)*station.SampleT;

if station.SampleF > 1
    DecimateR = round(station.SampleF); 
    nn = floor(record.NumPt/DecimateR);
    if (nn*DecimateR+1) <= record.NumPt
        ReSampleWave = decimate(SeisData(1:nn*DecimateR+1)', DecimateR);
    else
        ReSampleWave = decimate([SeisData(1:nn*DecimateR)'; SeisData(nn*DecimateR)], DecimateR);
    end    

    % ReSampleWave = decimate(SeisData', DecimateR);
    station.SampleF = 1;
    station.SampleT = 1;
    record.NumPt = length(ReSampleWave);
    record.SampleT = 1;
    record.SampleF = 1;
else
    ReSampleWave = SeisData';
end
ReSampleWave = ReSampleWave/max(ReSampleWave);
ampcoef = HdrData.SCALE; % amplitude scaling factor
record.Time = (record.NumPt - 1)*station.SampleT;       

%%%------------------------------------------------------------------------
%%%                 Filter Design Module  
%%%------------------------------------------------------------------------ 
% --- Executes on button press in Makefilter.
function Makefilter_Callback_old(hObject, eventdata, handles)
% hObject    handle to Makefilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global filter sta TBeta    
filter.Domain = get(handles.Filter_domain,'Value');
filter.Window = get(handles.Win_type,'Value');

% if  filter.SampleF == 0
% %     errordlg('No data sampling frequency information! Please check it!');
    prompt = {'Please enter your data sampling frequency (Hz):'};
    title = ['Data Sampling Frequency'];
    line = 2;
    DataSPF = inputdlg(prompt,title,line);
    filter.SampleF = str2num(DataSPF{1});
    filter.SampleT = 1/filter.SampleF;
% end

switch filter.Domain
    case 1
		prompt = {'Enter Band_Pass_Filter Test Central Period (e.g. 40 s):',
            'Enter Band_Pass_Filter Band Width (e.g. 1 s):'};
		title = ['Set Filter Parameter'];
		line = 2; 
		FilterPara = inputdlg(prompt,title,line);
        filter.CtrT = str2num(FilterPara{1});
        filter.BandWidth = str2num(FilterPara{2});
        filter.CtrF = (2/filter.SampleF)/filter.CtrT;
        filter.LowF = (2/filter.SampleF)/(filter.CtrT + 0.5*filter.BandWidth);
        filter.HighF = (2/filter.SampleF)/(filter.CtrT - 0.5*filter.BandWidth);    
    case 2
		prompt = {'Enter Band_Pass_Filter Test Central Frequency (e.g. 0.02 HZ):',
            'Enter Band_Pass_Filter Band Width (e.g. 0.005 HZ):'};
		title = ['Set Filter Parameter']; 
		line = 2; 
		FilterPara = inputdlg(prompt,title,line);
        filter.CtrF = str2num(FilterPara{1});
        filter.BandWidth = str2num(FilterPara{2});
        filter.LowF = (2/filter.SampleF)*(filter.CtrF - 0.5*filter.BandWidth);
        filter.HighF = (2/filter.SampleF)*(filter.CtrF + 0.5*filter.BandWidth);
end

switch filter.Window
    case 1
        prompt = {'Please input beta value of Kaiser window (e.g.: 5-10) (Increasing beta widens the main lobe and decreases the amplitude of the sidelobes'};
        title = ['Kaiser Window Parameter'];
        line = 2;
        KaiserBeta = inputdlg(prompt,title,line);
        filter.KaiserPara = str2num(KaiserBeta{1});   
    case 2
        prompt = {'Please input alfa value of Gaussian window (e.g.: 2-6) (Increasing alfa widens the main lobe and decreases the amplitude of the sidelobes'};
        title = ['Gaussian Window Parameter'];
        line = 2;
        GaussAlfa = inputdlg(prompt,title,line);
        filter.GaussAlfa = str2num(GaussAlfa{1});                
end        

% % define filter.KaiserPara according to filter.CtrT (central period of filtering)        
% filter.KaiserPara = interp1(TBeta(1,:), TBeta(2,:), filter.CtrT, 'pchip');
%
filter.Length = max(1024*filter.SampleF, round(6*filter.CtrT*filter.SampleF));

if filter.Window == 1   % Kaiser Window
    filter.Data = fir1(filter.Length, [filter.LowF, filter.HighF], kaiser(filter.Length + 1,filter.KaiserPara));
elseif filter.Window == 2   % Gaussian Window
    filter.Data = fir1(filter.Length, [filter.LowF, filter.HighF], gausswin(filter.Length + 1,filter.GaussAlfa));
end

% wvtool(filter.Data); % Dome window
set(handles.MsgEdit,'String','Set Band-pass Filter Successfully!');

function Makefilter_Callback(hObject, eventdata, handles)
% hObject    handle to Makefilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global filter    
filter.Domain = get(handles.Filter_domain,'Value');
filter.Window = get(handles.Win_type,'Value');


switch filter.Domain
    case 1
		prompt = {'Enter Band_Pass_Filter Band Width (e.g. 1 s):'};
		title = ['Set Filter Parameter'];
		line = 2; 
		FilterPara = inputdlg(prompt,title,line);
        filter.BandWidth = str2num(FilterPara{1}); 
    case 2
		prompt = {'Enter Band_Pass_Filter Band Width (e.g. 0.005 HZ):'};
		title = ['Set Filter Parameter']; 
		line = 2; 
		FilterPara = inputdlg(prompt,title,line);
        filter.BandWidth = str2num(FilterPara{1});
end

switch filter.Window
    case 1
        prompt = {'Please input beta value of Kaiser window (e.g.: 5-10) (Increasing beta widens the main lobe and decreases the amplitude of the sidelobes'};
        title = ['Kaiser Window Parameter'];
        line = 2;
        KaiserBeta = inputdlg(prompt,title,line);
        filter.KaiserPara = str2num(KaiserBeta{1});   
    case 2
        prompt = {'Please input alfa value of Gaussian window (e.g.: 2-6) (Increasing alfa widens the main lobe and decreases the amplitude of the sidelobes'};
        title = ['Gaussian Window Parameter'];
        line = 2;
        GaussAlfa = inputdlg(prompt,title,line);
        filter.GaussAlfa = str2num(GaussAlfa{1});                
end        

set(handles.MsgEdit,'String','Set Band-pass Filter Parameters Successfully!');

% --- Executes on button press in Dome_button.
function Dome_button_Callback(hObject, eventdata, handles)
% hObject    handle to Dome_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global filter

prompt = {'Please enter your data sampling frequency (Hz):'};
title = ['Data Sampling Frequency'];
line = 2;
DataSPF = inputdlg(prompt,title,line);
filter.SampleF = str2num(DataSPF{1});
filter.SampleT = 1/filter.SampleF;

switch filter.Domain
    case 1
		prompt = {'Enter Band_Pass_Filter Test Central Period (e.g. 40 s):'};
		title = ['Set Filter Parameter'];
		line = 2; 
		FilterPara = inputdlg(prompt,title,line);
        filter.CtrT = str2num(FilterPara{1});
        filter.CtrF = (2/filter.SampleF)/filter.CtrT;
        filter.LowF = (2/filter.SampleF)/(filter.CtrT + 0.5*filter.BandWidth);
        filter.HighF = (2/filter.SampleF)/(filter.CtrT - 0.5*filter.BandWidth);    
    case 2
		prompt = {'Enter Band_Pass_Filter Test Central Frequency (e.g. 0.02 HZ):',
            'Enter Band_Pass_Filter Band Width (e.g. 0.005 HZ):'};
		title = ['Set Filter Parameter']; 
		line = 2; 
		FilterPara = inputdlg(prompt,title,line);
        filter.CtrF = str2num(FilterPara{1});
        filter.LowF = (2/filter.SampleF)*(filter.CtrF - 0.5*filter.BandWidth);
        filter.HighF = (2/filter.SampleF)*(filter.CtrF + 0.5*filter.BandWidth);
end



filter.Length = max(1024*filter.SampleF, round(6*filter.CtrT*filter.SampleF));
if filter.Window == 1   % Kaiser Window
    filter.Data = fir1(filter.Length, [filter.LowF, filter.HighF], kaiser(filter.Length + 1,filter.KaiserPara));
elseif filter.Window == 2   % Gaussian Window
    filter.Data = fir1(filter.Length, [filter.LowF, filter.HighF], gausswin(filter.Length + 1,filter.GaussAlfa));
end
wvtool(filter.Data); % Dome window

%%%------------------------------------------------------------------------
%%%                 Nessary Files Input Module  
%%%------------------------------------------------------------------------ 
% Default path for convenience
function DefaultPath(handles)
global PhaseDispDirectory 
global EventlistDirectory 
%% Default 
CF=pwd;
cd ..;
DefaultPath=pwd;
cd(CF);
% DefaultPath='F:\Project\Two_Station_Analysis\TAnalysis_Auotmatic';
%% Default output fold
PhaseDispDirectory =fullfile(DefaultPath,'OutPut');
%% Default input stationlist
global allsta eventlist
stalistfold=fullfile(DefaultPath,'Necessary');
stalistfname='Turkey_StaList.dat';
stslistfile=fullfile(stalistfold,stalistfname);

fstat = fopen(stslistfile,'r');
if fstat>0
   allsta = struct('name', {}, 'net', {}, 'lat', {}, 'lon', {});
%  read stations
   i=1;
   while ~feof(fstat)
	    allsta(i).name = fscanf(fstat,'%s',1); %station name
	    allsta(i).net = fscanf(fstat,'%s',1); %station network
	    allsta(i).lat = fscanf(fstat,'%f',1); %station longitude
	    allsta(i).lon = fscanf(fstat,'%f',1); %station latitude   
        temp = fgetl(fstat);
        i=i+1;
   end
   fclose(fstat);
end
%% Default input eventlist 
EventlistDirectory=fullfile(DefaultPath,'Example_Data','Turkey_BHZ_Data');
ev_dir=dir(fullfile(EventlistDirectory,'eventlist*'));
defaultname=ev_dir(1).name;
eventlist_file=fullfile(EventlistDirectory,defaultname);
fname = fopen(eventlist_file,'r');
loni=0;
while ~feof(fname)
    loni=loni+1;
    eventlist(loni).filename=fgetl(fname);
end
fclose(fname);

%% Default refence file
global ref_file
ref_fold=fullfile(DefaultPath,'Necessary');
ref_name='Rayleigh_phase_v_ref.dat';
ref_file=fullfile(ref_fold,ref_name);
 Load_reference_range(handles); 
%% Default Filter parameters
global filter
filter.BandWidth=1;
filter.KaiserPara=8;
%% Default Resp file path
global RespDirectory
RespDirectory=fullfile(DefaultPath,'Example_Data','Turkey_BHZ_Data','Turkey_RESP_PZ');
%% Default waveformfile path
global DataDirectory
DataDirectory=fullfile(DefaultPath,'Example_Data','Turkey_BHZ_Data','201203201802');


% --- Executes on button press in StationFileInput. (Station list path)
function StationFileInput_Callback(hObject, eventdata, handles)%
% hObject    handle to StationFileInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
clear allsta
global allsta 
[pfile1, pname1, index] = uigetfile({'*.*';'*.txt';'*.dat'},'Open file containing all stations information', pwd);
files = fullfile(pname1,pfile1);
if pname1==0
   set(handles.MsgEdit,'String','File containing all stations information read unsuccseefully !');
else
fstat = fopen(files,'r');
allsta = struct('name', {}, 'net', {}, 'lat', {}, 'lon', {});

% read stations
i=1;
while ~feof(fstat)
	   allsta(i).name = fscanf(fstat,'%s',1); 
	   allsta(i).net = fscanf(fstat,'%s',1);
	   allsta(i).lat = fscanf(fstat,'%f',1); 
	   allsta(i).lon = fscanf(fstat,'%f',1);   
       temp = fgetl(fstat);
      i=i+1;
end
fclose(fstat);

%plot topography, coast, and all stations
PlotAllStations(allsta, handles);
end

% --- Executes on button press in DisperFolder. (Output dispersion file
% path)
function DisperFolder_Callback(hObject, eventdata, handles)
% hObject    handle to DisperFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global PhaseDispDirectory
PhaseDispDirectory = uigetdir(pwd, 'Select folder to save phase velocity dispersion data:');
set(handles.MsgEdit,'String',['Phase velocity dispersion data folder is: ',PhaseDispDirectory]);

% --- Executes on button press in WaveformFileInput.(Single-event mode)
function WaveformFileInput_Callback(hObject, eventdata, handles)
% hObject    handle to WaveformFileInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DataDirectory src MonthDays evdata allsta 
cla(handles.axes1);
cla(handles.axes3);
cla(handles.axes4);
cla(handles.axes2);

if isempty(allsta)
    words='Please load the station list first!';
    errordlg(words,'Warning!');
else    

evdata = struct('fname', '', 'staname', '');

if DataDirectory == 0
   DataDirectory = pwd;
end

[pfile1, pname1, index] = uigetfile({'*.txt';'*.dat';'*.*'},'Open file containing all data file names', DataDirectory);
files = [pname1,pfile1];
DataDirectory = pname1;

if pname1==0
   set(handles.MsgEdit,'String','WaveformFile read unsuccseefully !');
else
							
fname = fopen(files);
temp = fscanf(fname, '%s', 1);
temp = fscanf(fname, '%s', 1);
src.YY = str2num(temp(1:4));
Month = str2num(temp(6:7));
Day = str2num(temp(9:10));
src.Month = Month;
src.Day = Day;
if mod(src.YY, 4) == 0 
    src.DD = sum(MonthDays(1,1:(Month - 1))) + Day;
else
    src.DD = sum(MonthDays(2,1:(Month - 1))) + Day;  
end
temp = fscanf(fname, '%s', 1);
src.HH = str2num(temp(1:2));
src.MM = str2num(temp(4:5));
src.SS = str2num(temp(7:8));
src.MS = 10*str2num(temp(10:11));
temp = fscanf(fname, '%s', 1);
mchar = size(temp,2);
src.Lat = str2num(temp(1:(mchar - 1)));
temp = fscanf(fname, '%s', 1);
mchar = size(temp,2);
src.Lon = str2num(temp(1:(mchar - 1)));
temp = fscanf(fname, '%s', 4);
temp = fscanf(fname, '%s', 1);
mchar = size(temp,2);
if strcmp(temp(mchar),',') == 1        
    srcmag = str2num(temp(1:(mchar - 1)));
else
    srcmag = str2num(temp(1:mchar));
end
temp = fgetl(fname);
   
i = 1;
while ~feof(fname)
    name = fscanf(fname, '%s', 1);
    if ~strcmp(name,'')
        evdata(i).fname = name;
        
        tmpchar = evdata(i).fname(25:32);
        index = 1;
        Pt = zeros(1,8);
        for j = 1:length(tmpchar)
            if strcmp(tmpchar(j),'.')==1
                Pt(index) = j;
                index = index+1;
            end           
        end
        evdata(i).staname = tmpchar((Pt(1)+1):(Pt(2)-1));        
        i =  i+1;
    else
        break
    end
end
filenum = length(evdata);
fclose(fname);
set(handles.EditSrcLon,'String',num2str(src.Lon));
set(handles.EditSrcLat,'String',num2str(src.Lat));
PlotAllStations(allsta, handles);
PlotEveStations(evdata, allsta, handles);
end

end

% --- Executes on button press in EventList. (event cycle mode)
function EventList_Callback(hObject, eventdata, handles)
% hObject    handle to EventList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear eventlist
global EventlistDirectory eventlist allsta


cla(handles.axes3);
cla(handles.axes4);
cla(handles.axes2);

if isempty(allsta)
    words='Please load the station list first!';
    errordlg(words,'Warning!');
else    




if EventlistDirectory == 0
   EventlistDirectory = pwd;
end
[pfile1, pname1, index] = uigetfile({'*.txt';'*.dat';'*.*'},'Open file containing all event fold names', EventlistDirectory);
files =[pname1,pfile1];
EventlistDirectory = pname1; 

%
if pname1==0
   set(handles.MsgEdit,'String','File containing all event fold names read unsuccseefully !');
else
% read the file containing event fold
   fname = fopen(files);
   loni=1;
   while ~feof(fname)
         eventlist(loni).filename=fgetl(fname);
         loni=loni+1;
   end
   fclose(fname);
   set(handles.MsgEdit,'String','File containing all event fold names read succseefully !');
end 
end

% --- load the waveformfile list (called by LaunchProcessing)
function signal = Silent_WaveformFileInput(hObject, eventdata, handles,eventfold_filename)
% hObject    handle to WaveformFileInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global  src MonthDays evdata allsta EventlistDirectory DataDirectory

evdata = struct('fname', '', 'staname', '');%
%%% find the file
event_path=fullfile(EventlistDirectory,eventfold_filename);
keyword=[eventfold_filename(1:4),'*File.txt']; % 2004*File.txt
file_list=dir(fullfile(event_path,keyword)); % should only one
% if length(file_list)~=1
%     disp(['the event fold ',eventfold_filename,' has no filelist or two txt files!!!']);
%     set(handles.StopDataProcessing,'Value',1);
% end
file_name=file_list.name;
DataDirectory=event_path; % every event fold path
waveformfile=fullfile(event_path,file_name);
if ~exist(waveformfile,'file')
    signal=0;
    fid_log=fopen('TSA_log.txt','a+');
    fprintf(fid_log,'\n Input waveformfile do not exist!!! the event is : %s',eventfold_filename);
    fclose(fid_log);
else
%%% read the file
fname = fopen(waveformfile);
% read weed format event information
temp = fscanf(fname, '%s', 1);
signal=strcmp(temp,'ISC/ISC,');
if signal==0
   fid_log=fopen('TSA_log.txt','a+');
   fprintf(fid_log,'\n Input waveformfile do not have the first line earthquake information,waveformifle is %s',waveformfile);
   fclose(fid_log);
else
   temp = fscanf(fname, '%s', 1);
   src.YY = str2num(temp(1:4));
   Month = str2num(temp(6:7));
   Day = str2num(temp(9:10));
   src.Month = Month;
   src.Day = Day;
   if mod(src.YY, 4) == 0  
      src.DD = sum(MonthDays(1,1:(Month - 1))) + Day;
   else
      src.DD = sum(MonthDays(2,1:(Month - 1))) + Day;  
   end
   temp = fscanf(fname, '%s', 1);
   src.HH = str2num(temp(1:2));
   src.MM = str2num(temp(4:5));
   src.SS = str2num(temp(7:8));
   src.MS = 10*str2num(temp(10:11));
   temp = fscanf(fname, '%s', 1);
   mchar = size(temp,2);
   src.Lat = str2num(temp(1:(mchar - 1)));
   temp = fscanf(fname, '%s', 1);
   mchar = size(temp,2);
   src.Lon = str2num(temp(1:(mchar - 1)));
   temp = fscanf(fname, '%s', 4);
   temp = fscanf(fname, '%s', 1);
   mchar = size(temp,2);
   if strcmp(temp(mchar),',') == 1        
      srcmag = str2num(temp(1:(mchar - 1)));
   else
      srcmag = str2num(temp(1:mchar));
   end
   temp = fgetl(fname);
   i = 1;
   while ~feof(fname)
    name = fscanf(fname, '%s', 1);
    if ~strcmp(name,'')
        evdata(i).fname = name;
        
        tmpchar = evdata(i).fname(25:32);
        index = 1;
        Pt = zeros(1,8);
        for j = 1:length(tmpchar)
            if strcmp(tmpchar(j),'.')==1
                Pt(index) = j;
                index = index+1;
            end           
        end
        evdata(i).staname = tmpchar((Pt(1)+1):(Pt(2)-1));      
        i =  i+1;
    else
        break
    end
   end
   
   filenum = length(evdata);
   fclose(fname);
   set(handles.EditSrcLon,'String',num2str(src.Lon));
   set(handles.EditSrcLat,'String',num2str(src.Lat));
   PlotAllStations(allsta, handles);
   PlotEveStations(evdata, allsta, handles);
   %%%%%%%%
   set(handles.MsgEdit,'String','WaveformFile read succseefully !');
end % signal
end % exist

%%%------------------------------------------------------------------------
%%%                 Infromation Renew Module   
%%%------------------------------------------------------------------------ 
 % updata message box information
function UpdataMsgBoxInfo(hObject, handles,index1,index2)
global src sta
set(handles.EditSrcLon,'String',num2str(src.Lon));
set(handles.EditSrcLat,'String',num2str(src.Lat));

% if sta(1).GCDkm < sta(2).GCDkm
%     set(handles.EditNameSta1,'String',num2str(sta(1).Name));
%     set(handles.EditNameSta2,'String',num2str(sta(2).Name));
%     set(handles.EditLonSta1,'String',num2str(sta(1).Lon));
%     set(handles.EditLonSta2,'String',num2str(sta(2).Lon));
%     set(handles.EditLatSta1,'String',num2str(sta(1).Lat));
%     set(handles.EditLatSta2,'String',num2str(sta(2).Lat));
%     set(handles.EditSrcStaDist1,'String',num2str(round(sta(1).GCDkm)));
%     set(handles.EditSrcStaDist2,'String',num2str(round(sta(2).GCDkm)));
%     set(handles.EditSrcStaAzim1,'String',num2str(sta(1).Azim));
%     set(handles.EditSrcStaAzim2,'String',num2str(sta(2).Azim));
% else
%     set(handles.EditNameSta1,'String',num2str(sta(2).Name));
%     set(handles.EditNameSta2,'String',num2str(sta(1).Name));
%     set(handles.EditLonSta1,'String',num2str(sta(2).Lon));
%     set(handles.EditLonSta2,'String',num2str(sta(1).Lon));
%     set(handles.EditLatSta1,'String',num2str(sta(2).Lat));
%     set(handles.EditLatSta2,'String',num2str(sta(1).Lat));
%     set(handles.EditSrcStaDist1,'String',num2str(round(sta(2).GCDkm)));
%     set(handles.EditSrcStaDist2,'String',num2str(round(sta(1).GCDkm)));
%     set(handles.EditSrcStaAzim1,'String',num2str(sta(2).Azim));
%     set(handles.EditSrcStaAzim2,'String',num2str(sta(1).Azim));
% end
    set(handles.EditIndex1,'String',num2str(index1));
    set(handles.EditIndex2,'String',num2str(index2));
    set(handles.EditNameSta1,'String',num2str(sta(1).Name));
    set(handles.EditNameSta2,'String',num2str(sta(2).Name));
    set(handles.EditLonSta1,'String',num2str(sta(1).Lon));
    set(handles.EditLonSta2,'String',num2str(sta(2).Lon));
    set(handles.EditLatSta1,'String',num2str(sta(1).Lat));
    set(handles.EditLatSta2,'String',num2str(sta(2).Lat));
    set(handles.EditSrcStaDist1,'String',num2str(round(sta(1).GCDkm)));
    set(handles.EditSrcStaDist2,'String',num2str(round(sta(2).GCDkm)));
    set(handles.EditSrcStaAzim1,'String',num2str(sta(1).Azim));
    set(handles.EditSrcStaAzim2,'String',num2str(sta(2).Azim));
%%%------------------------------------------------------------------------
%%%                 Process Control Module  
%%%------------------------------------------------------------------------ 
%    
function StartDataProcessing_Callback(hObject, eventdata, handles)
% hObject    handle to StartDataProcessing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% StationFileInput
 LaunchProcessing(hObject, eventdata, handles);
 
% --- Executes on button press in ProcessAllData.
function ProcessAllData_Callback(hObject, eventdata, handles)
% hObject    handle to ProcessAllData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ProcessAllData

if get(handles.ProcessAllData, 'Value')
    set(handles.ProcessSelectedData, 'Value', 0);
else
    set(handles.ProcessSelectedData, 'Value', 1);
end

% --- Executes on button press in ProcessSelectedData.
function ProcessSelectedData_Callback(hObject, eventdata, handles)
% hObject    handle to ProcessSelectedData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ProcessSelectedData

if get(handles.ProcessSelectedData, 'Value')
    set(handles.ProcessAllData, 'Value', 0);
else
    set(handles.ProcessAllData, 'Value', 1);
end

% --- Executes on button press in Exitbutton.
function Exitbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Exitbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Signal
set(handles.MsgEdit,'String','Eixt the process.');
Signal=0;


% --- Executes on button press in Pausebutton.
function Pausebutton_Callback(hObject, eventdata, handles)
% hObject    handle to Pausebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiwait();

% --- Executes on button press in Goonbutton.
function Goonbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Goonbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%%%%%%%
uiresume();
%%%------------------------------------------------------------------------
%%%                 Reference Setting Module  
%%%------------------------------------------------------------------------ 
% --- Executes on button press in RefC_Disper. (For needing reference mode)
function RefC_Disper_Callback(hObject, eventdata, handles)
% hObject    handle to RefC_Disper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global  ref_file
% input phase velocity reference dispersion and ranges
[pfile1, pname1, index] = uigetfile({'*.dat';'*.txt';'*.*'},'Open reference phase velocity file', pwd);
ref_file = fullfile(pname1,pfile1);
if pname1==0
   set(handles.MsgEdit,'String','Reference phase velocity dispersion read unsuccessfully!');
   set(handles.RefC_Disper,'Value',0);
else
   if get(handles.RefC_Disper,'Value')==1
      Load_reference_range(handles); 
  end
end

% --- Load the refernce and plot
function Load_reference_range(handles)

global g_refphasedisp ref_file

refcdisp=load(ref_file);
g_refphasedisp = refcdisp;

[nT, nc] = size(g_refphasedisp);
if nc == 3
   g_refphasedisp(:,2) = refcdisp(:,2) - refcdisp(:,3);
   g_refphasedisp(:,3) = refcdisp(:,2) + refcdisp(:,3);
   set(handles.MsgEdit,'String','Reference phase velocity dispersion read successfully!');
elseif nc == 2
   dc = str2num(get(handles.deltaPhaseV, 'String'));
   g_refphasedisp(:,2) = refcdisp(:,2) - dc;
   g_refphasedisp(:,3) = refcdisp(:,2) + dc;
   set(handles.MsgEdit,'String','Reference phase velocity dispersion read successfully!');
end

h2 = handles.axes2;
set(gcf,'CurrentAxes',h2);
hold(h2,'off');
plot(g_refphasedisp(:,1), g_refphasedisp(:,2),'b--');
hold on; plot(g_refphasedisp(:,1), g_refphasedisp(:,3),'b--');
plot(g_refphasedisp(:,1), refcdisp(:,2), 'b', 'LineWidth',2);
set(gca,'YTick',2.5:0.2:5.5);
set(gca,'XTick',20:10:150);
xlim([20 150]);
ylim([2.5 5.5]);
xlabel('Period (sec)', 'FontSize',11,'FontWeight','bold','fontname','Times New Roman','VerticalAlignment','middle');
ylabel('Ref Phase Velocity (km/s)', 'FontSize',11,'FontWeight','bold','fontname','Times New Roman');

function deltaPhaseV_Callback(hObject, eventdata, handles)
% hObject    handle to deltaPhaseV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global g_refphasedisp
if (isempty(g_refphasedisp))
  set(handles.MsgEdit,'String','The change of the refernece range will be adopted when the refernce file does not contain the error range.');
else
  Load_reference_range(handles);
end

% --- Executes on button press in RefStaPairFolder.(Optional)
function RefStaPairFolder_Callback(hObject, eventdata, handles)
% hObject    handle to RefStaPairFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% note1. the station pair reference file must begging with 'Ref'
% note2. the station pair files must cover all of the station-pair which
%        mean Num(files)=Num(station)(Num(station)-1)/2 
% note3. the period range reference dierperison must be equal to the GUI
%        input Period Range:From handles.StartPeriod to handles.EndPeriod
global allsta RefstaData
RefstaData=struct('name',{},'data',{});
RefstaPFold = uigetdir(pwd, 'Open the folder saving station-pair reference phase velocity dispersion data:');
filelist=dir(fullfile(RefstaPFold,'Ref*'));
stapairNum=length(allsta)*(length(allsta)-1)/2;

if length(filelist)==stapairNum
   set(handles.MsgEdit,'String',['Phase velocity dispersion data folder is: ',RefstaPFold,'Please Wait For Loading!']);
   for loni=1:length(filelist)
       refDataIn=load(fullfile(RefstaPFold,filelist(loni).name));
       RefstaData(loni).name=strrep(strrep(filelist(loni).name,'Ref_',''),'.dat','');
       [nr, nc] = size(refDataIn);
       refData=refDataIn;
       if nc == 3
          refData(:,2) = refDataIn(:,2) - refDataIn(:,3);
          refData(:,3) = refDataIn(:,2) + refDataIn(:,3);
          RefstaData(loni).data=refData;
       elseif nc == 2
          dc = str2num(get(handles.deltStaPairC, 'String'));
          refData(:,2) = refDataIn(:,2) - dc;
          refData(:,3) = refDataIn(:,2) + dc;
          RefstaData(loni).data=refData;
       end
       set(handles.MsgEdit,'String','Reference of station-pair phase velocity dispersion read successfully!');
       clear refData refDataIn
   end
else
   set(handles.MsgEdit,'String','Reference of station-pair phase velocity dispersion read Unsuccessfully!');
end


%%%------------------------------------------------------------------------
%%%                 TSAnalysis Mode Swich Module  
%%%------------------------------------------------------------------------ 
% Automatic Events Cycle Mode
% --- Executes on button press in Visible_cycle_mode. (Visible Mode)
function Visible_cycle_mode_Callback(hObject, eventdata, handles)
% hObject    handle to Visible_cycle_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of Visible_cycle_mode
 if  get(handles.Visible_cycle_mode,'Value')
     set(handles.Silent_cycle_mode, 'Value',0);
     set(handles.Semi_single_mode,'Value',0);
     set(handles.Auto_single_mode,'Value',0);
     set(handles.MsgEdit,'String','Enter Automatic Event Cycle Visble Mode!');
 end
 
% --- Executes on button press in Silent_cycle_mode.(Silent Mode)
function Silent_cycle_mode_Callback(hObject, eventdata, handles)
% hObject    handle to Silent_cycle_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 if get(handles.Silent_cycle_mode,'Value')
    set(handles.Silent_cycle_mode, 'Value',0);
    set(handles.Semi_single_mode,'Value',0);
    set(handles.Auto_single_mode,'Value',0);
    set(handles.MsgEdit,'String','Enter Automatic Event Cycle Silent Mode!');
 end


% Single Event Mode
% --- Executes on button press in Auto_single_mode.
function Auto_single_mode_Callback(hObject, eventdata, handles)
% hObject    handle to Auto_single_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 if get(handles.Auto_single_mode,'Value')
    set(handles.Silent_cycle_mode, 'Value',0);
    set(handles.Visible_cycle_mode,'Value',0);
    set(handles.Semi_single_mode,'Value',0);
    set(handles.MsgEdit,'String','Enter Single Event Seim-automatic Mode!');
 end
 
% --- Executes on button press in Semi_single_mode.
function Semi_single_mode_Callback(hObject, eventdata, handles)
% hObject    handle to Semi_single_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 if get(handles.Semi_single_mode,'Value')
    set(handles.Silent_cycle_mode, 'Value',0);
    set(handles.Visible_cycle_mode,'Value',0);
    set(handles.Auto_single_mode,'Value',0);
    set(handles.MsgEdit,'String','Enter Single Event Automatic Mode!');
 end

%%%------------------------------------------------------------------------
%%%                 Remove Instrument Response Module   
%%%------------------------------------------------------------------------  

% --- Remove instrument response: Main function
function RmInstruResponse_new(handles)
global sta rcd_Z wave RespDirectory
global allresp

% Setting
% find the corresponding files
respnum = length(allresp);
RespFile = cell(1,2);
IndexRespFile = zeros(1,2);
fileindex=zeros(1,2);
for k = 1:2
    for i = 1:respnum
        if strcmp(allresp(i).staname,sta(k).Name)
            RespFile{k} = fullfile(RespDirectory, allresp(i).fname);
            IndexRespFile(k) = 1;
            % get the type of resp file
            tmpchar = allresp(i).fname(1:end);
            index = 1;
            Pt = zeros(1,8);
            for j = 1:length(tmpchar)
                if strcmp(tmpchar(j),'.')==1
                   Pt(index) = j;
                   index = index+1;
                end           
            end
            if strcmpi(tmpchar(1:Pt(1)-1),'RESP')
               fileindex(k)=1;
            elseif strcmpi(tmpchar(1:Pt(1)-1),'PZ')
               fileindex(k)=2;
            end
                       
            break;
        end
    end
end    

if IndexRespFile(1) == 0  
    error('Error: response file for station 1 is not existing!');
end
if IndexRespFile(2) == 0 
    error('Error: response file for station 2 is not existing!');
end

% filter range
TLow = 0.8*str2double(get(handles.StartPeriod,'String'));  
TMin = 0.5*str2double(get(handles.StartPeriod,'String'));  
THigh = 1.2*str2double(get(handles.EndPeriod,'String')); 
TMax = 1.5*str2double(get(handles.EndPeriod,'String')); 
HighF = 1/TLow;  % high pass freq  0.0625Hz  f3
HighFMax = min(1/TMin, sta(k).SampleF/2); 
LowF = 1/THigh;   % low pass freq   0.0056  f2
LowFMin = max(1/TMax, 0); % lowest freq for response removal  0.0044  f1
freqrange=[LowFMin,LowF,HighF,HighFMax];

water_level_deconvolution = 0.001;% 0.001,for remove instrument response

% remove the instrument response
% fileindex = 1 : read RESP.* response file
%           = 2 : polezeros response file
% freqrange = [f_low_cut f_low_pass f_high_pass f_high_cut] % filter freq band
%              f_low_cut < f_low_pass < f_high_pass < f_high_cut, e.g.,
%              [0.005 0.01 5 10] Hz

for k = 1:2
    fs=rcd_Z(k).SampleF;
    seisdata=wave(k).DatZ;    
    wave(k).DatZ =rmResp_bpfilter(seisdata, fs, freqrange,RespFile{k}, fileindex(k), water_level_deconvolution);    
end

% --- function to process 
function bfwave = rmResp_bpfilter(seisdata, fs, freqrange, respfile, fileindex, water_level_deconvolution)
% seisdata: data vector
% fs: sampling frequency
% freqrange = [f_low_cut f_low_pass f_high_pass f_high_cut] % filter freq band
%              f_low_cut < f_low_pass < f_high_pass < f_high_cut, e.g.,
%              [0.005 0.01 5 10] Hz
% respfile: response file name
% fileindex = 1 : read RESP.* response file
%           = 2 : polezeros response file

staNumPt = length(seisdata); % number of points in data
bfwave=zeros(staNumPt,1);

%% demean and detrend the data
seisdata = detrend(seisdata - mean(seisdata)); 

%% remove instrument response   
if fileindex == 1 | fileindex == 2
    
    % define instrument response information
    resp = struct('Amp',0,...
        'NumZeros',0,...
        'NumPoles',0,...
        'Zeros',0,...
        'Poles',0,...
        'poly_num',0,...
        'poly_den',0);
    
    %%
    HighF = freqrange(3);  % high pass freq
    LowF = freqrange(2);   % low pass freq
    %% HighF must be less than fs/2, otherwise set to be 0.99*(fs/2)
    if HighF > (fs/2)
        HighF = 0.99*(fs/2);
        display(['HighF error! HighF is reset to ', num2str(HighF)]);
    end
    
    %% read instrument response files
    if fileindex == 1 % read RESP.* response file
        [resp.Amp, resp.Numzeros, resp.Numpoles, resp.Zeros, resp.Poles] = Rd_InstruRespFile(respfile);
    elseif fileindex == 2 % read PoleZero* response file
        [resp.Amp, resp.Numzeros, resp.Numpoles, resp.Zeros, resp.Poles] = Rd_InstruPoleZeroFile(respfile);
    end
    
    [resp.poly_num, resp.poly_den] = zp2tf(resp.Zeros', resp.Poles', resp.Amp);
    
    LowFMin = max(freqrange(1), 0); % lowest freq for response removal
    HighFMax = min(freqrange(4), fs/2); % highest freq for reponse removal
    
    SegLength = 2^20; % length of data parts for fft to avoid memory issue
    fftnum = ceil(staNumPt/SegLength); % number of segments for the data
    
    %%
    
    for k=0:(fftnum-1)
        %%
        if k ~=(fftnum-1)
            datapart = detrend(seisdata((1 + k*SegLength):(k + 1)*SegLength));
            fftlength=SegLength;
            fftdata=fft(datapart,fftlength);
            clear datapart
        else
            datapart = detrend(seisdata((1 + k*SegLength):staNumPt));
            fftlength = 2^(nextpow2(staNumPt - k*SegLength));
            fftdata=fft(datapart, fftlength);
            clear datapart
        end
        fftdata = reshape(fftdata, 1, fftlength);
        %%
        f(1:(fftlength/2+1)) = fs*(0:(fftlength/2))/fftlength;
        delta_f = fs/fftlength;
        
        w(1:(fftlength/2+1)) = 2*pi*f(1:(fftlength/2+1));
        h = freqs(resp.poly_num, resp.poly_den, w); % obtain instrument response
        
        MinFPoint = max(2, ceil(LowFMin/delta_f));
        MaxFPoint = min(fftlength/2, floor(HighFMax/delta_f));
        
        % remove instrument response: the first half frequency spectrum
        nn =  MinFPoint:MaxFPoint;
        % Y = XH -> X = Y/H -> X = Y*conj(H)/abs(H)^2
        for i_nn=nn(1):nn(end)
            h2_abs(i_nn)=max(abs(h(i_nn)).^2,water_level_deconvolution);
        end
        clear i_nn;
        fftdata(nn) = fftdata(nn).*conj(h(nn))./h2_abs(nn);  % water-level deconvolution
        fftdata(1:MinFPoint) = 0;
        fftdata(MaxFPoint:(fftlength/2+1)) = 0;
        fftdata((fftlength/2+2):fftlength)=conj(fftdata((fftlength/2):-1:2)); % treat another half spectrum
        
        %% band pass filtering
        LowPtN = round(LowF/delta_f);
        HighPtN = round(HighF/delta_f);
        nptdfs = round((LowF - LowFMin)/delta_f);
        if nptdfs >= 4
            nn = (LowPtN - nptdfs):(LowPtN-1);
            taperwin = hann(2*nptdfs-1)';
            fftdata(1:(LowPtN - nptdfs -1))=0;
            fftdata(nn) = taperwin(1:nptdfs).*fftdata(nn);
        end
        
        nptdfs = round((HighFMax - HighF)/delta_f);
        nn = (HighPtN + 1):(HighPtN + nptdfs);
        if nptdfs >= 4
            taperwin = hann(2*nptdfs-1)';
            fftdata(nn)= taperwin(nptdfs:end).*fftdata(nn);
            fftdata((HighPtN + nptdfs + 1):(fftlength/2+1)) = 0;
        end
        
        fftdata((fftlength/2+2):fftlength)=conj(fftdata((fftlength/2):-1:2));
        
        
        %% time domain filtered data
        band_data_T=real(ifft(fftdata));
        
        if k~=(fftnum-1)
            bfwave((1 + k*SegLength):(k + 1)*SegLength)=band_data_T;
        else
            bfwave((1 + k*SegLength):staNumPt)=band_data_T(1:(staNumPt - k*SegLength));
        end
        
    end
else
    bfwave = seisdata;
end

% --- Read  read RESP.* response file
function [Normfactor, Numzeros, Numpoles, Respzero, Resppole] = Rd_InstruRespFile(RespFile)
%%
fname = fopen(RespFile,'r');

%read Normfactor_A0 @ line 19
isok=0;
while (isok == 0)
    temp1 = fgetl(fname);
    if (length(findstr('Response out units lookup:',temp1)) > 0)
        temp2 = fscanf(fname,'%s',4);
        Normfactor_A0 = fscanf(fname,'%f',1);
        temp3 = fgetl(fname);
        isok=1;
    end
end

%read line 20
temp1 = fgetl(fname);
%read line 21
temp1 = fscanf(fname,'%s',4);
Numzeros = fscanf(fname,'%f',1);
temp1 = fgetl(fname);
%read line 22
temp1 = fscanf(fname,'%s',4);
Numpoles = fscanf(fname,'%f',1);
temp1 = fgetl(fname);
%read line 23, 24: zeroes header
temp1 = fgetl(fname);
temp1 = fgetl(fname);
%read zeros
Respzero = zeros(1, Numzeros);
for i = 1:Numzeros
   temp1 = fscanf(fname,'%s',1);
   temp = fscanf(fname,'%d',1);
   realpart = fscanf(fname,'%e',1);
   imagpart = fscanf(fname,'%e',1);
   Respzero(i) = complex(realpart, imagpart);
   temp1 = fgetl(fname);
end
%read 2 lines: poles header
temp1 = fgetl(fname);
temp1 = fgetl(fname);
%read poles
Resppole = zeros(1, Numpoles);
for i = 1:Numpoles
   temp1 = fscanf(fname,'%s',1);
   temp = fscanf(fname,'%d',1);
   realpart = fscanf(fname,'%e',1);
   imagpart = fscanf(fname,'%e',1);
   Resppole(i) = complex(realpart, imagpart);
   temp1 = fgetl(fname);
end

%read Normfactor_Sensitivity @ line end-3
isok=0;
while (isok == 0)
    temp1 = fgetl(fname);
    if (length(findstr('Stage sequence number:',temp1)) > 0)
        temp2 = fscanf(fname,'%s',2);
        if (length(findstr('Sensitivity:',temp2)) > 0)
            Normfactor_Sensitivity = fscanf(fname,'%f',1);
            temp3 = fgetl(fname);
            isok=1;
        else
            temp3 = fgetl(fname);
        end
    end
end

fclose(fname);

Normfactor=Normfactor_A0*Normfactor_Sensitivity;

% --- read PoleZero* response file
function [Normfactor, Numzeros, Numpoles, Respzero, Resppole] = Rd_InstruPoleZeroFile(RespFile)

fname = fopen(RespFile,'r');

temp = fscanf(fname, '%s', 1);

while strcmp(temp(1), '*')    
   temp4 = fgetl(fname);
   temp = fscanf(fname, '%s', 1);
end

if strcmp(temp, 'ZEROS')
    Numzeros = fscanf(fname, '%d', 1);
    temp1 = fgetl(fname);
    Respzero = zeros(1, Numzeros);
    for i = 1:Numzeros
       realpart = fscanf(fname,'%e',1);
       imagpart = fscanf(fname,'%e',1);
       RespzeroTemp(i) = complex(realpart, imagpart);
       temp1 = fgetl(fname);    
    end    
end
Numzeros = Numzeros - 1; % do not convert to displacement 
Respzero = RespzeroTemp(2:end);

temp = fscanf(fname, '%s', 1);
if strcmp(temp, 'POLES')
    Numpoles = fscanf(fname, '%d', 1);
    temp1 = fgetl(fname);
    Resppole = zeros(1, Numpoles);
    for i = 1:Numpoles
       realpart = fscanf(fname,'%e',1);
       imagpart = fscanf(fname,'%e',1);
       Resppole(i) = complex(realpart, imagpart);
       temp1 = fgetl(fname);    
    end
end

temp = fscanf(fname, '%s', 1);
Normfactor_CONSTANT = fscanf(fname, '%e', 1);
fclose(fname);

Normfactor=Normfactor_CONSTANT;


% --- Remove instrument response
function RmInstruResponse_old(handles)

global sta rcd_Z wave RespDirectory
global allresp

h2 = handles.axes2;

respnum = length(allresp);
RespFile = cell(1,2);
IndexRespFile = zeros(1,2);
display('The two station response files: ');
for k = 1:2
    for i = 1:respnum
        if strcmp(allresp(i).staname,sta(k).Name)
            RespFile{k} = fullfile(RespDirectory, allresp(i).fname);
            IndexRespFile(k) = 1;
            display(['  ' num2str(k) ':  ' allresp(i).fname]);
            break;
        end
    end
end            

if IndexRespFile(1) == 0  
    display('Error: response file for station 1 is not existing!');
end
if IndexRespFile(2) == 0 
    display('Error: response file for station 2 is not existing!');
end
    
%% remove instrument response and bandpass filtering
for k = 1:2
    
    [resp(k).Amp, resp(k).Numzeroes, resp(k).Numpoles, resp(k).Zeros, resp(k).Poles] = Rd_InstruRespFile(RespFile{k});
    [resp(k).poly_num, resp(k).poly_den] = zp2tf(resp(k).Zeros', resp(k).Poles', resp(k).Amp);
    
    if mod(rcd_Z(k).NumPt,2) == 1
        rcd_Z(k).NumPt = rcd_Z(k).NumPt + 1;
        wave(k).DatZ(rcd_Z(k).NumPt) = 0;
    end
    fftlength = rcd_Z(k).NumPt;
    fftdata = fft(wave(k).DatZ, fftlength);
    
    f(1:(fftlength/2+1)) = sta(k).SampleF*(0:(fftlength/2))/fftlength;


    % The band pass filtering: [TMin TLow THigh TMax] <-> [0.0 1.0 1.0 0.0]
    TLow = 0.8*str2double(get(handles.StartPeriod,'String'));  
    TMin = 0.5*str2double(get(handles.StartPeriod,'String'));  

    THigh = 1.2*str2double(get(handles.EndPeriod,'String')); 
    TMax = 1.5*str2double(get(handles.EndPeriod,'String'));  

    HighF = 1/TLow;  % high pass freq  0.0625Hz  f3
    HighFMax = min(1/TMin, sta(k).SampleF/2); 

    LowF = 1/THigh;   % low pass freq   0.0056  f2
    LowFMin = max(1/TMax, 0); % lowest freq for response removal  0.0044  f1
    
%     figure
%     subplot(2,1,1), plot(f(1:(fftlength/2+1)), abs(fftdata(1:(fftlength/2+1))));
%     subplot(2,1,2),plot(f(1:(fftlength/2+1)), angle(fftdata(1:(fftlength/2+1))));
    
    delta_f = sta(k).SampleF/fftlength;
    w(1:(fftlength/2+1)) = 2*pi*f(1:(fftlength/2+1));
    h = freqs(resp(k).poly_num, resp(k).poly_den, w); 

    %% remove instrument response  f1-f4
    MinFPoint = max(2, ceil(LowFMin/delta_f));
    MaxFPoint = min(fftlength/2, floor(HighFMax/delta_f)); 
    nn =  MinFPoint:MaxFPoint;
    % Y = XH -> X = Y/H -> X = Y*conj(H)/abs(H)^2
    h = h/max(abs(h));
    fftdata = reshape(fftdata, 1, fftlength);
    fftdata(nn) = fftdata(nn).*conj(h(nn))./(abs(h(nn)).^2 + 0.01);  
    fftdata(1:MinFPoint) = 0;
    fftdata(MaxFPoint:(fftlength/2+1)) = 0;
    fftdata((fftlength/2+2):fftlength)=conj(fftdata((fftlength/2):-1:2)); 
    
    %% band pass filtering f1-f2  f3-f4
    LowPtN = round(LowF/delta_f);%f2
    HighPtN = round(HighF/delta_f); % f3
    nptdfs = round((LowF - LowFMin)/delta_f);
    if nptdfs >= 4
        nn = (LowPtN - nptdfs):(LowPtN-1);
        taperwin = hann(2*nptdfs-1)';
        fftdata(1:(LowPtN - nptdfs -1))=0; 
        % figure(99); hold off; subplot(2,1,1); hold off; plot(abs(fftdata(nn)),'r');
        fftdata(nn) = taperwin(1:nptdfs).*fftdata(nn);
        % hold on; plot(abs(fftdata(nn)),'b--');
    end

    nptdfs = round((HighFMax - HighF)/delta_f);%  f3-f4
    nn = (HighPtN + 1):(HighPtN + nptdfs);
    if nptdfs >= 4
        taperwin = hann(2*nptdfs-1)';
        % subplot(2,1,2); hold off; plot(abs(fftdata(nn)),'r');
        fftdata(nn)= taperwin(nptdfs:end).*fftdata(nn);   
        fftdata((HighPtN + nptdfs + 1):(fftlength/2+1)) = 0;
        % hold on; plot(abs(fftdata(nn)),'b--');
    end

    fftdata((fftlength/2+2):fftlength)=conj(fftdata((fftlength/2):-1:2));        
    fftdata = reshape(fftdata, fftlength, 1);

    wave(k).DatZ = real(ifft(fftdata,fftlength));
    
%     hold(h2, 'off');
%     if k == 1
%         plot(h2, -rcd_Z(k).DiffT+(0:(rcd_Z(k).NumPt-1))*sta(k).SampleT, wave(k).DatZ+2, 'r');
%     elseif k==2
%         hold(h2, 'on');
%         plot(h2, -rcd_Z(k).DiffT+(0:(rcd_Z(k).NumPt-1))*sta(k).SampleT, wave(k).DatZ, 'r');
%     end
 
end



%%%------------------------------------------------------------------------
%%%                  Plotting Module   
%%%------------------------------------------------------------------------
% --- Appendix plot C-T spectrum
function plotC_Tspectrumh2(handle,clim)
global cross
TPoint = cross.StartT:cross.DeltaT:cross.EndT;
VPoint = cross.VMin:cross.DeltaV:cross.VMax;
CImgPt = size(VPoint,2);

X=TPoint(1:cross.NumCtrT);
Y=VPoint(1:CImgPt);
C=cross.PhasVImg(1:CImgPt,1:cross.NumCtrT);

set(gcf,'CurrentAxes',handle);
hold(handle,'off');
imagesc(X,Y,C,clim); 
set(gca,'YDir','normal');
xlabel('Period (s)', 'FontSize', 10, 'FontWeight', 'bold','fontname','Times New Roman');
ylabel('Phase Velocity (km/s)', 'FontSize', 10, 'FontWeight', 'bold','fontname','Times New Roman');
set(gca, 'FontSize', 10, 'FontWeight', 'bold');
xlim(gca,[cross.StartT,cross.EndT]);
set(gca, 'XTick',cross.StartT:10:cross.EndT,'YTick',2:0.25:6,'XGrid','on','YGrid','on');

%%%------------------------------------------------------------------------
%%%                 Functional Module  
%%%------------------------------------------------------------------------ 
% --- Plot topography, coast, and all stations
function PlotAllStations_old(allsta, handles)

h1 = handles.axes1;
hold(h1,'off');
set(gcf,'CurrentAxes',h1)
load topomap
x=linspace(0,360,1080);
y=linspace(-90,90,2160);
[cmap clim] = demcmap(topomap);
hold on
imagesc(x,y,topomap,clim);%colorbar('vert');
colormap(cmap);axis image; grid on; axis on;

load coast % the coast line is lat, long
kk = find(long < 0);
long(kk) = 360 - abs(long(kk));
ii = find( abs(long) < 0.5);
long(ii) = NaN;
lat(ii) = NaN;
hold(h1,'on');
plot(h1, long,lat,'k', 'LineWidth',2);
set(gca,'ydir','normal');
set(gca,'FontSize',16,'FontWeight','bold');
axis equal

stanum = length(allsta);
for i = 1:stanum
    Lon(i) = allsta(i).lon;
    Lat(i) = allsta(i).lat;
    plot(h1, allsta(i).lon, allsta(i).lat, 'k^', 'MarkerSize',6, 'MarkerFaceColor','k');
end

ExtraLon = 0.2*(max(Lon)-min(Lon));
ExtraLat = 0.2*(max(Lat)-min(Lat));
ExtraLon = min(ExtraLon, 1);
ExtraLat = min(ExtraLat, 1);
LonMin=min(Lon)-ExtraLon;
LonMax=max(Lon)+ExtraLon;
LatMin=min(Lat)-ExtraLat;
LatMax=max(Lat)+ExtraLat;
xlim(h1,[LonMin LonMax]);
ylim(h1,[LatMin LatMax]);
set(gca, 'XTickMode','auto','YTickMode','auto');
% set(gca, 'XTick',LonMin:round((LonMax-LonMin)/5):LonMax);
% set(gca, 'YTick',LatMin:round((LatMax-LatMin)/5):LatMax);
% title('Station Location', 'FontSize',16,'FontWeight','bold');
stanum=i-1;

set(handles.MsgEdit,'String','Read Stations Information Successfully!');

function PlotAllStations(allsta, handles)

h1 = handles.axes1;
hold(h1,'off');
set(gcf,'CurrentAxes',h1)
load topomap
x=linspace(0,360,1080);
y=linspace(-90,90,2160);
[cmap clim] = demcmap(topomap);
hold on
imagesc(x,y,topomap,clim);%colorbar('vert');
colormap(cmap);axis image; grid on; axis on;

load coast % the coast line is lat, long
kk = find(long < 0);
long(kk) = 360 - abs(long(kk));
ii = find( abs(long) < 0.5);
long(ii) = NaN;
lat(ii) = NaN;
hold(h1,'on');
plot(h1, long,lat,'k', 'LineWidth',2);
set(gca,'ydir','normal');
set(gca,'FontSize',16,'FontWeight','bold');
axis equal

stanum = length(allsta);
for i = 1:stanum
    Lon(i) = allsta(i).lon;
    Lat(i) = allsta(i).lat;
    hold(h1, 'on');
    plot(h1, allsta(i).lon, allsta(i).lat, 'k^', 'MarkerSize',6, 'MarkerFaceColor','k');
end

ExtraLon = 0.1*(max(Lon)-min(Lon));
ExtraLat = 0.1*(max(Lat)-min(Lat));
ExtraLon = min(ExtraLon, 1);
ExtraLat = min(ExtraLat, 1);
xlim(h1,[min(Lon)-ExtraLon max(Lon)+ExtraLon]);
ylim(h1,[min(Lat)-ExtraLat max(Lat)+ExtraLat]);
set(gca, 'XTickMode','auto','YTickMode','auto');
% title('Station Location', 'FontSize',16,'FontWeight','bold');
stanum=i-1;

set(handles.MsgEdit,'String','Read Stations Information Successfully!');

% --- Plot the stations with waveform file loaded
function PlotEveStations(evdata, allsta, handles)

h1 = handles.axes1;
hold(h1,'on');
set(gcf,'CurrentAxes',h1);

filenum = length(evdata);
stanum = length(allsta);
for i = 1:filenum
    % get station lat and lon from station name
    for kk = 1:stanum        
        if strcmp(evdata(i).staname,allsta(kk).name)
            stalat = allsta(kk).lat;
            stalon = allsta(kk).lon;
            hold(h1, 'on');
            plot(h1, stalon, stalat, 'r^', 'MarkerSize',8, 'MarkerFaceColor','r');
            hold(h1, 'on');
            text(stalon, stalat, evdata(i).staname, 'FontSize',8);
            break
        end
    end
end

set(handles.MsgEdit,'String','Read Waveform Data File List Successfully!');

% --- Judge whether the earthquake is almost on the interstation great
% circle path with alfa < selminangle1 (e.g. 3 deg), beta < selminangle2 (e.g. 6 deg), 
% please refer to Yao et al.,2006, GJI paper for the detail illustration of angles alfa and beta
% And minimum event-station distance is MinEveStaDist (e.g.,1000 km) and
% minimum differential distance from event to two stations is MinDiffDist
% (e.g., 50 km)
function isTSpath = IsAlongTSGCPath(eventlat, eventlon, rstalat1, rstalon1, rstalat2, rstalon2, hObject, handles)

srdist1 = deg2km(distance([eventlat, eventlon], [rstalat1, rstalon1]));
srdist2 = deg2km(distance([eventlat, eventlon], [rstalat2, rstalon2]));
MinEveStaDist = 1000;  % minimum event to station distance
MinDiffDist = 50;      % minimum differential distance from event to two stations

if min(srdist1, srdist2) >= MinEveStaDist && abs(srdist1 - srdist2) >= MinDiffDist
    selminangle1 = str2num(get(handles.IsTSPathAlfa,'String'));
    selminangle2 = str2num(get(handles.IsTSPathBeta,'String'));

    az1 = azimuth([eventlat, eventlon], [rstalat1, rstalon1]);
    az2 = azimuth([eventlat, eventlon], [rstalat2, rstalon2]);

    if srdist2 > srdist1
        az3 = azimuth([rstalat1, rstalon1], [rstalat2, rstalon2]);
        az4 = azimuth([rstalat1, rstalon1], [-eventlat, eventlon+180]);
    else
        az3 = azimuth([rstalat2, rstalon2], [rstalat1, rstalon1]);
        az4 = azimuth([rstalat2, rstalon2], [-eventlat, eventlon+180]);
    end

    angle1 = min(abs(az1-az2), abs(abs(az1-az2)-360));
    angle2 = min(abs(az3-az4), abs(abs(az3-az4)-360));

    % on a great cirle(deviation angles alpha & beta less than selminangle degrees)
    if angle1 <= selminangle1 && angle2 <= selminangle2    
        isTSpath = 1;
    else
        isTSpath = 0;
    end
else
    isTSpath = 0;
end    

% --- Executes in function StartDataProcessing_Callback
function MatchSignal=ReferenceMatch(staName1,staName2)
% this function is to match the station-pair name to RefstaData and produce
% the global variable g_refphasedisp
global g_refphasedisp RefstaData
% for lonj = 1:max(size(staName1,2), size(staName2,2))
%     if double(staName1(lonj)) > double(staName2(lonj))
%        StaKey = strcat(staName1,'-',staName1);
%        break
%     elseif double(staName1(lonj)) < double(staName2(lonj))
%        StaKey = strcat(staName2,'-', staName1);
%        break
%     else
%        continue
%     end
% end
% stapairNum=length(RefstaData);
% MatchSignal=0;
% for loni=1:stapairNum
%     if strcmp(StaKey,RefstaData(loni).name)
%        g_refphasedisp=RefstaData(loni).data;
%        MatchSignal=1;
% %      disp(['Find station-pair ref',staName1,staName2])
%     end   
% end
StaKey = strcat(staName1,'-',staName2);
MatchSignal=0;
stapairNum=length(RefstaData);
for loni=1:stapairNum
    if strcmp(StaKey,RefstaData(loni).name)
       g_refphasedisp=RefstaData(loni).data;
       MatchSignal=1;
       break;
    end   
end

if MatchSignal==0
   StaKey = strcat(staName2,'-',staName1);
   for loni=1:stapairNum
       if strcmp(StaKey,RefstaData(loni).name)
          g_refphasedisp=RefstaData(loni).data;
          MatchSignal=1;
          break;
       end   
   end
end
if MatchSignal==0
   disp([' Can not Find station-pair ref file:',StaKey])
   g_refphasedisp=[]; 
end

% --- Write the log file
function loglog(words)
global TSAlog
fid_log=fopen(TSAlog,'a+');
fprintf(fid_log,'%s \n',words);
fclose(fid_log);

% --- Check the necessary files
function pass=CheckNecessaries(handles)
global allsta allresp g_refphasedisp evdata eventlist filter
% Check the necessary
pass=1;
if isempty(allsta)
   words= 'Please load the station list file!';
   errordlg(words,'Warning'); 
   pass=0;
end

if isempty(evdata) && (isempty(eventlist))
   words='Please load the event list file for cycle mode or waveform list file for single event mode! ';
   errordlg(words,'Warning'); 
   pass=0;
end   

if get(handles.Resp_switch,'Value')==1 && (isempty(allresp))
   words= 'If you choose to remove instrument response, please load the resp list file!';
   errordlg(words,'Warning');
   pass=0;
end

if isempty(filter.KaiserPara)
   words= 'Please pess ''Make Filter '' botton to set the parameters of fliter!';
   errordlg(words,'Warning');
   pass=0;
end

if get(handles.Semi_single_mode,'Value')==0 && (isempty(g_refphasedisp))
   words= 'Please pess ''Ref C Disper '' botton to load the reference dispersion curve which is needed in the automatic mode!';
   errordlg(words,'Warning');
   pass=0;
end     

function Default_path(handles)
global filter
filter.KaiserPara=9;

%%%------------------------------------------------------------------------
%%%                 TSAnalysis Main Program   
%%%------------------------------------------------------------------------  
% --- Launch
% --- Executes on button press in StartDataProcessing.
function  LaunchProcessing(hObject, eventdata, handles)
% hObject    handle to StartDataProcessing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% StationFileInput
global allsta 
% RespFileInput
global  allresp
% WaveformFileInput
global DataDirectory src evdata TSAlog
global wave  
% EventList_Callback
global eventlist EventlistDirectory PhaseDispDirectory

global Signal


Default_path(handles);
% sac file btyes limits, if the file is cutted the lim must be small
btyes_lim=3000;
h2 = handles.axes2;
tic;

Signal=1;
Date_now=date;
TSAlog=['TSAnalysisLog_',Date_now,'.txt'];

% Choose the TSanalysis Mode
cycleSignal=2;
if get(handles.Silent_cycle_mode,'Value') || get(handles.Visible_cycle_mode,'Value')
   cycleSignal=1;
elseif get(handles.Semi_single_mode,'Value') || get(handles.Auto_single_mode,'Value')
   cycleSignal=0;
elseif cycleSignal==2
   words= 'Please choose a mode!';
   errordlg(words,'Warning'); 
end
 
% cycle events mode
if cycleSignal~=2
   pass=CheckNecessaries(handles);
   if pass==0
      cycleSignal=2;
   end
end

if cycleSignal==1
%%% preparation: 1.stationinput 2.resp input 3.wavefile input
   event_Num=length(eventlist);
   words=[datestr(now,0),' CycleMode'];
   disp(words);
   loglog(words);
   words=['Output path:',PhaseDispDirectory];
   loglog(words);
   for index=1:event_Num  
       if Signal==0
          set(handles.MsgEdit,'String','Data Processing Stopped!');
          break
       end           
       words=['NO.',num2str(index),'/',num2str(event_Num),': earthquake event date:',eventlist(index).filename];
       disp(words);
       loglog(words);
% wavefile input
       fold_path=fullfile(EventlistDirectory,eventlist(index).filename);
       if exist(fold_path,'dir')==0
          words=['the event ',fold_path,'is Null, go to the next event!'];
          loglog(words);
          continue
       end       
       WF_signal=Silent_WaveformFileInput(hObject, eventdata, handles,eventlist(index).filename);
       if WF_signal==0
          continue
       end   
       stanum = length(allsta);
       respnum = length(allresp);
       filenum = length(evdata);
% Initialization
       StartIndex1 = 1;
       EndIndex1 = filenum - 1;
       EndIndex2 = filenum;          
% prossess begin  
       for i = StartIndex1:EndIndex1   
           if Signal==0
              set(handles.MsgEdit,'String','Data Processing Stopped!');
              break
          end          
          StartIndex2 = i + 1;           
          for j = StartIndex2:EndIndex2
              if Signal==0
                 break
              end 
% get station lat and lon from station name
              for kk = 1:stanum
                  if strcmp(evdata(i).staname,allsta(kk).name)
                     stalat1 = allsta(kk).lat;
                     stalon1 = allsta(kk).lon;
                  end
                  if strcmp(evdata(j).staname,allsta(kk).name)
                     stalat2 = allsta(kk).lat;
                     stalon2 = allsta(kk).lon;
                  end            
              end
% judge whether the event and two stations are satisfying two-station criterion 
% to avoid the same (or nearly the same) two stations%dist = distance(lat1,lon1,lat2,lon2)
              if deg2km(distance(stalat1,stalon1,stalat2,stalon2)) > 1 
                 isTSpath = IsAlongTSGCPath(src.Lat, src.Lon, stalat1, stalon1, stalat2, stalon2, hObject, handles);
              else
                 isTSpath = 0;
%                fid_log=fopen(TSAlog,'a+');
%                fprintf(fid_log,'\n not a qualified two-station path:\n %s \n %s',evdata(i).fname,evdata(j).fname);
%                fclose(fid_log);
              end
              if isTSpath == 1
                 seisfile1 = fullfile(DataDirectory, evdata(i).fname);
                 seisfile2 = fullfile(DataDirectory, evdata(j).fname);
                 file_sta1=dir(seisfile1);
                 file_sta2=dir(seisfile2);
                 if isempty(file_sta1) || isempty(file_sta2)
                    fid_log=fopen(TSAlog,'a+');
                    fprintf(fid_log,'The one of station files is Null,\n %s \n %s',seisfile1,seisfile2);
                    fclose(fid_log)
                 elseif (file_sta1.bytes< btyes_lim) || (file_sta2.bytes<btyes_lim)
                    fid_log=fopen(TSAlog,'a+');
                    fprintf(fid_log,'The one of station files is too short:\n %s bytes is %f \n %s bytes is %f ',seisfile1,file_sta1.bytes,seisfile2,file_sta2.bytes);
                    fclose(fid_log);
                 else
                    % read data for two station 
                    RdTwoStaData(seisfile1, seisfile2);   
% For the adding specific reference dispersion curve.                    
%                     if get(handles.SpecificRefIndex,'Value')
%                        MatchSignal=ReferenceMatch(evdata(i).staname,evdata(j).staname);
%                        if MatchSignal==0
%                           set(handles.MsgEdit,'String','For AutoPhaseDisperProcessing lack the corresponding reference phase velocity!!! ');
%                           words=['For station pair : ',evdata(i).staname,'-',evdata(j).staname,' Lack corresponding reference dispersion file of the stations! skip this pair!'];
%                           loglog(TSAlog,words);
%                           error(words);
%                        end
%                     end
                    SumAmp1 = sum(abs(wave(1).DatZ(1:end)));
                    SumAmp2 = sum(abs(wave(2).DatZ(1:end)));
                    if min(SumAmp1,SumAmp2) > 0
                       set(gcf,'CurrentAxes',h2);            
                       if get(handles.Visible_cycle_mode,'Value')==1
                           UpdataMsgBoxInfo(hObject, handles,i,j);
                       end
                       % main function for two-station data processing
                       [ProcessIndex,warning]=TSDataProcessing(hObject, handles);
                       if ProcessIndex~=1
                           word=['File: ',evdata(i).fname,' ',evdata(j).fname,' Warning:',warning];
                           loglog(word);
                           disp(word);
                       end
                    end % TSD main fuction
                 end  % length
              end % isTSpath      
          end %%(j = StartIndex2:EndIndex2)
      end %% i = StartIndex1:EndIndex1 
  end %index=1:event_Num
%%% -----------------------------------------------------------------------

elseif cycleSignal==0
  stanum = length(allsta);
  respnum = length(allresp);
  filenum = length(evdata);
%%% Initialization
  if get(handles.ProcessAllData, 'Value')  
     StartIndex1 = 1;
     EndIndex1 = filenum - 1;
     EndIndex2 = filenum;
     set(handles.DataStartIndex1, 'String', num2str(StartIndex1));
     set(handles.DataStartIndex2, 'String', 'i+1');
     set(handles.DataEndIndex1, 'String', num2str(EndIndex1));
     set(handles.DataEndIndex2, 'String', num2str(EndIndex2));
  else
     StartIndex1 = str2num(get(handles.DataStartIndex1, 'String'));
     StartIndex2 = str2num(get(handles.DataStartIndex2, 'String'));
     EndIndex1 = str2num(get(handles.DataEndIndex1, 'String'));
     EndIndex2 = str2num(get(handles.DataEndIndex2, 'String'));
  end
% new plot
  PlotAllStations(allsta, handles);
  PlotEveStations(evdata, allsta, handles);
% prossess begin  
  words=[datestr(now,0),' Single Event Mode'];
  disp(words);
  words=['Output path:',PhaseDispDirectory];
  disp(words);
  for i = StartIndex1:EndIndex1
      set(handles.MsgEdit,'String',['First Statio Index = ' num2str(i)]);    
      if Signal==0
         set(handles.MsgEdit,'String','Data Processing Stopped!');
         break
      end    
      if get(handles.ProcessAllData, 'Value')
         StartIndex2 = i + 1;
      elseif strcmp('i+1',get(handles.DataStartIndex2, 'String'))
         StartIndex2 = i + 1;
      end    
      for j = StartIndex2:EndIndex2
          if Signal==0
             break
          end  
% get station lat and lon from station name
         for kk = 1:stanum
             if strcmp(evdata(i).staname,allsta(kk).name)
                stalat1 = allsta(kk).lat;
                stalon1 = allsta(kk).lon;
             end
             if strcmp(evdata(j).staname,allsta(kk).name)
                stalat2 = allsta(kk).lat;
                stalon2 = allsta(kk).lon;
             end            
         end
% judge whether the event and two stations are satisfying two-station criterion
% to avoid the same (or nearly the same) two stations %dist = distance(lat1,lon1,lat2,lon2)
         if deg2km(distance(stalat1,stalon1,stalat2,stalon2)) > 1 
            isTSpath = IsAlongTSGCPath(src.Lat, src.Lon, stalat1, stalon1, stalat2, stalon2, hObject, handles);
         else
            isTSpath = 0;
            % set(handles.MsgEdit,'String',[evdata(i).staname '-' evdata(j).staname ': not a qualified two-station path']);
            display([evdata(i).staname '-' evdata(j).staname ': not a qualified two-station path']);
         end
         if isTSpath == 1
            seisfile1 = fullfile(DataDirectory, evdata(i).fname);
            seisfile2 = fullfile(DataDirectory, evdata(j).fname);
            display('The two station data files: ');
            display(['  1:  ' evdata(i).fname]);
            display(['  2:  ' evdata(j).fname]);  
            % read data for two station 
            RdTwoStaData(seisfile1, seisfile2); 
% For the adding specific reference dispersion curve. 
%             if get(handles.SpecificRefIndex,'Value')
%                MatchSignal=ReferenceMatch(evdata(i).staname,evdata(j).staname);
%                if MatchSignal==0
%                   set(handles.MsgEdit,'String','For AutoPhaseDisperProcessing lack the reference phase velocity!!! ');
%                   words=['For station pair : ',evdata(i).staname,'-',evdata(j).staname,' Lack corresponding reference dispersion file of the stations! skip this pair!'];
%                   loglog(TSAlog,words);
%                   error(words);
%                end
%             end            
            SumAmp1 = sum(abs(wave(1).DatZ(1:end)));
            SumAmp2 = sum(abs(wave(2).DatZ(1:end)));
            if min(SumAmp1,SumAmp2) > 0
               set(gcf,'CurrentAxes',h2);          
               UpdataMsgBoxInfo(hObject, handles,i,j);
               % main function for two-station data processing
               [ProcessIndex,warning]=TSDataProcessing(hObject, handles);
               if ProcessIndex~=1
                  word=['File: ',evdata(i).fname,' ',evdata(j).fname,' Warning:',warning];
                  disp(word);
               end
               set(handles.MsgEdit,'String','Data Processing is continuing!');             
            end 
         end      
      end %(j = StartIndex2:EndIndex2)
   end % i = StartIndex1:EndIndex1
end % silent 
set(handles.MsgEdit,'String','Data Processing Finished!');
toc;


% --- Main function for two-station data processing
function [ProcessIndex,warning]=TSDataProcessing(hObject, handles)
global cross  g_refphasedisp
%% 1.Remove instrument response
if get(handles.Lovewave,'Value')==1
   cross.WaveType = 2;
else
   cross.WaveType = 1;  
end

if get(handles.Resp_switch,'Value')==1
   RmInstruResponse_new(handles);
%    words='Setp 1.Instrument Response Removed!';
%    set(handles.MsgEdit,'String',words);
end

%% 2.Get group image and group arrival,MFT to get the group arrival
ProcessIndex = GroupImageArrival(hObject, handles);
if ProcessIndex~=1
   warning='Surface wave window not long enough for the processing!';
else 
   warning='Normal Staus.';
end

if get(handles.Silent_cycle_mode,'Value')~=1 && get(handles.Visible_cycle_mode,'Value')~=1
   set(handles.MsgEdit,'String','Setp 2.MFT,Group arrival determined.');
   pause(0.5);
end

%% 3.Moving window the waves with respect to the group arrival and do cross-correlation at each narrow period band   
if ProcessIndex == 1  
   CrossCorrelation(hObject, handles);
   if get(handles.Silent_cycle_mode,'Value')==0
      set(handles.MsgEdit,'String','Setp 3.Moving-window and cross-correlation with narrowband filter.');
   end
   IsDispGood = 1;  
%% 4.For pick up phase velocity
% IsDispGood=1.repetition(semiautomatic);2.save dispersion
%            0.bad non dipsersion save(automatic)
   while IsDispGood == 1               
%    plot cross dot C-T spectrum %modification,move from crosscorrelation
     if get(handles.Silent_cycle_mode,'Value')~=1
        h2 = handles.axes2;
        plotC_Tspectrumh2(h2,[-1, 1]);
     end      
     
     if get(handles.Semi_single_mode, 'Value')~=1 
        if isempty(g_refphasedisp) 
           set(handles.MsgEdit,'String','For AutoPhaseDisperProcessing lack the reference phase velocity!!! ');
           error('For AutoPhaseDisperProcessing lack the reference phase velocity!!! ');
        else
           IsDispGood = AutoPhaseDisperProcessing(hObject, handles);
           if get(handles.Auto_single_mode,'Value')==1
              pause(1);
           elseif get(handles.Visible_cycle_mode,'Value')==1
              pause(0.5);
           end
        end     
     else
        IsDispGood = PhaseVDisper(hObject, handles); 
     end  
   end  

end 



% --- Step 1. Get the group velocity arrival time
%     group.ArrPt1
function ProcessIndex = GroupImageArrival(hObject, handles)
global sta wave rcd_Z src filter cross DeltaTInitial
global TravPtV G_VPoint
global group
% global evestadist1 evestadist2 % add for the near station distance

% you are expecting to obtain the phase velocity dispersion curve from 20 s to 100 s with 1 s interval. 
cross.StartT = str2num(get(handles.StartPeriod,'String'));
cross.EndT = str2num(get(handles.EndPeriod,'String'));
EndT = cross.EndT;
cross.DeltaT = str2num(get(handles.DeltaPeriod,'String'));

filter.Domain = get(handles.Filter_domain,'Value');
filter.Window = get(handles.Win_type,'Value');

h2 = handles.axes2; 
h3 = handles.axes3; 
h4 = handles.axes4;

GVWinR = [2.5 5.2]; % group velocity window (km/s) for Rayleigh Wave
GVWinL = [2.8 5.8]; % group velocity window (km/s) for Love Wave
RefTravTMin = [0 0];
RefTravTMax = [0 0];
%%%   count the time difference from the the time recorded
if cross.WaveType == 1 
    for i = 1:2
        RefTravTMin(i) = sta(i).GCDkm/GVWinR(2) + rcd_Z(i).DiffT;
        RefTravTMax(i) = sta(i).GCDkm/GVWinR(1) + rcd_Z(i).DiffT;  
    end
else
    for i = 1:2
        RefTravTMin(i) = sta(i).GCDkm/GVWinL(2) + rcd_Z(i).DiffT;
        RefTravTMax(i) = sta(i).GCDkm/GVWinL(1) + rcd_Z(i).DiffT;
    end
end

% Please  make sure the surface wave travel time in the range of the record!

if (min(RefTravTMin) >= 0) && (min(rcd_Z(1).Time - RefTravTMax(1), rcd_Z(2).Time - RefTravTMax(2)) > 0)
    ProcessIndex = 1;   % continue processing
    
    StaDistance = abs(sta(1).GCDkm - sta(2).GCDkm);
    cross.EndT = min(EndT, 2*StaDistance/2.5);
    cross.EndT = round(cross.EndT - mod(cross.EndT, cross.DeltaT));

    if cross.EndT <= cross.StartT
        ProcessIndex = 0;  % stop processing
    else
        TPoint = cross.StartT:cross.DeltaT:cross.EndT;
        cross.NumCtrT = length(TPoint);
    end
else
    ProcessIndex = 0;  % stop processing
end
%% --- Begin to process 
if ProcessIndex == 1
%% --- 1. Cut the surface waveform    
   % set filter sample frequency to be the data sampling frequency
   filter.SampleF = rcd_Z(1).SampleF;% 1Hz
   filter.SampleT = 1/filter.SampleF;
   
   TaperTime = 50*min(sta(1).GCDkm, sta(2).GCDkm)/2000; % unit: second   
   TaperNum = round(TaperTime*rcd_Z(1).SampleF);
   [CrossDeltaT, WaveIndex] = max(RefTravTMax - RefTravTMin);
   
   cross.StartNum = ceil(RefTravTMin./[filter.SampleT filter.SampleT]);
   cross.EndNum = floor(RefTravTMax./[filter.SampleT filter.SampleT]);
   cross.SENum = cross.EndNum - cross.StartNum + 1;
   cross.PointNum = max(cross.SENum);
   
   % for converting  travel time to group velocity
   for loni=1:2
     time = (cross.StartNum(loni):( cross.StartNum(loni)+cross.PointNum-1 ) )*filter.SampleT-rcd_Z(loni).DiffT;
     TravPtV(loni).v_inv=sta(loni).GCDkm./time;
   end
   
  % extract the effective surface waves
   MaxFilterLength = round(max(512*filter.SampleF, 5*cross.EndT*filter.SampleF));
   MaxHalfFilterNum =  floor(MaxFilterLength/2);
   cross.wave1 = zeros(1, cross.PointNum + MaxHalfFilterNum);
   cross.wave2 = zeros(1, cross.PointNum + MaxHalfFilterNum);
   
   % set taper window
   window1 = ones(cross.SENum(1),1);
   window2 = ones(cross.SENum(2),1);
   window1(1:TaperNum(1)) = sin(0.5*pi*(1:TaperNum)/TaperNum);
   window1((cross.SENum(1)-TaperNum+1):cross.SENum(1)) = window1(TaperNum:-1:1);
   window2(1:TaperNum) = sin(0.5*pi*(1:TaperNum)/TaperNum);
   window2((cross.SENum(2)-TaperNum+1):cross.SENum(2)) = window2(TaperNum:-1:1);
   
   % obtain waveforms for group arrival analysis and cross-correlation
   if sta(1).GCDkm < sta(2).GCDkm
       DeltaTInitial = rcd_Z(1).DiffT - rcd_Z(2).DiffT + (cross.StartNum(2)-1)*rcd_Z(2).SampleT - (cross.StartNum(1)-1)*rcd_Z(1).SampleT;
       cross.wave1(1:cross.SENum(1)) = wave(1).DatZ(cross.StartNum(1):cross.EndNum(1)).*window1;
       cross.wave2(1:cross.SENum(2)) = wave(2).DatZ(cross.StartNum(2):cross.EndNum(2)).*window2;
       
       % plot the wavefrom image % 
       if get(handles.Silent_cycle_mode,'Value')~=1
          hold(h2,'off')
          plot(h2,(1:cross.PointNum)*rcd_Z(1).SampleT, 2+cross.wave1(1:cross.PointNum)/max(cross.wave1(1:cross.PointNum)),'k');
          hold(h2,'on');
          plot(h2,(1:cross.PointNum)*rcd_Z(2).SampleT + DeltaTInitial, cross.wave2(1:cross.PointNum)/max(cross.wave2(1:cross.PointNum)),'k');
       end  
       evestadist1 =  sta(1).GCDkm;
       evestadist2 =  sta(2).GCDkm;
   else
       DeltaTInitial = rcd_Z(2).DiffT - rcd_Z(1).DiffT + (cross.StartNum(1)-1)*rcd_Z(1).SampleT - (cross.StartNum(2)-1)*rcd_Z(2).SampleT;
       cross.wave1(1:cross.SENum(2)) = wave(2).DatZ(cross.StartNum(2):cross.EndNum(2)).*window2;
       cross.wave2(1:cross.SENum(1)) = wave(1).DatZ(cross.StartNum(1):cross.EndNum(1)).*window1;
       % adjustment
       group_v=TravPtV(2).v_inv;
       TravPtV(2).v_inv=TravPtV(1).v_inv;
       TravPtV(1).v_inv=group_v;
                   
       if get(handles.Silent_cycle_mode,'Value')~=1
          hold(h2,'off');
          plot(h2,(1:cross.PointNum)*rcd_Z(2).SampleT, 2+cross.wave1(1:cross.PointNum)/max(cross.wave1(1:cross.PointNum)),'k');
          hold(h2,'on');
          plot(h2,(1:cross.PointNum)*rcd_Z(1).SampleT + DeltaTInitial, cross.wave2(1:cross.PointNum)/max(cross.wave2(1:cross.PointNum)),'k');
       end 
       evestadist2 =  sta(1).GCDkm;
       evestadist1 =  sta(2).GCDkm;      
   end
%% --- 2. Obtain wave group image
   %  set period vector
   TPoint = cross.StartT:cross.DeltaT:cross.EndT;
   Ypoint = 1:cross.PointNum;

   % bandpass filtering the waveforms and obtain wave group image    
   crossgroup1 = EnvelopeImageCalculation(cross.wave1(1:cross.PointNum), filter.SampleF, TPoint, evestadist1);
   AmpS_T1 = max(crossgroup1,[],2);
   crossgroup2 = EnvelopeImageCalculation(cross.wave2(1:cross.PointNum), filter.SampleF, TPoint, evestadist2);
   AmpS_T2 = max(crossgroup2,[],2);
   
   for numt = 1:cross.NumCtrT
       cross.group1(1:cross.PointNum, numt) = crossgroup1(numt,:)'/AmpS_T1(numt);
       cross.group2(1:cross.PointNum, numt) = crossgroup2(numt,:)'/AmpS_T2(numt);
   end
   % clear crossgroup1 AmpS_T1 crossgroup2 AmpS_T2 
   % clear  window1 window2  
   
   group.ArrPt1 = zeros(1, cross.NumCtrT);
   group.ArrPt2 = zeros(1, cross.NumCtrT);
   % The final result of this function, using in the next part.
   [MaxAmp,group.ArrPt1(1:cross.NumCtrT)] = max(cross.group1(1:cross.PointNum,1:cross.NumCtrT));
   [MaxAmp,group.ArrPt2(1:cross.NumCtrT)] = max(cross.group2(1:cross.PointNum,1:cross.NumCtrT)); 
   
%% --- 3. Covert arrival time to group velocity for figure
   DeltaV=0.002;
   for loni=1:2 % 1 for the near station,2 for far
     G_VPoint(loni).v_inv = TravPtV(loni).v_inv(1):-1*DeltaV:TravPtV(loni).v_inv(end);
     G_VPoint(loni).Num = length(G_VPoint(loni).v_inv);
     G_VPoint(loni).v = G_VPoint(loni).v_inv(end:-1:1);
   end

   for numt = 1:cross.NumCtrT     
        cross.GroupVImg1(1:G_VPoint(1).Num, numt) = interp1( TravPtV(1).v_inv ,crossgroup1(numt,:)'/AmpS_T1(numt),G_VPoint(1).v_inv );
        cross.GroupVImg2(1:G_VPoint(2).Num, numt) = interp1( TravPtV(2).v_inv ,crossgroup2(numt,:)'/AmpS_T2(numt),G_VPoint(2).v_inv );
   end 
   % obtain the group arrival by search the max amplitude in the velocity
    group.ArrPoint1 = zeros(1, cross.NumCtrT);
    group.ArrPoint2 = zeros(1, cross.NumCtrT);
    group.ArrV1 = zeros(1,cross.NumCtrT);
    group.ArrV2 = zeros(1,cross.NumCtrT);
    [MaxAmp,group.ArrPoint1(1:cross.NumCtrT)] = max(cross.GroupVImg1(1:G_VPoint(1).Num,1:cross.NumCtrT));
    [MaxAmp,group.ArrPoint2(1:cross.NumCtrT)] = max(cross.GroupVImg2(1:G_VPoint(2).Num,1:cross.NumCtrT)); 
    group.ArrV1(1:cross.NumCtrT)= G_VPoint(1).v_inv( group.ArrPoint1(1:cross.NumCtrT) );
    group.ArrV2(1:cross.NumCtrT)= G_VPoint(2).v_inv( group.ArrPoint2(1:cross.NumCtrT) ); 
  
    clear crossgroup1 AmpS_T1 crossgroup2 AmpS_T2    
    clear window1 window2 
    clear Travtime group_v  
    clear time
%% --- 4. plot part.
if get(handles.Silent_cycle_mode,'Value')~=1
   % the nearer station
   set(gcf,'CurrentAxes',h3);
   hold(h3,'off');
   % imagesc(TPoint,Ypoint,cross.group1(1:cross.PointNum,1:cross.NumCtrT),[0,1]); 
   minamp = min(min(cross.GroupVImg1(1:G_VPoint(1).Num,1:cross.NumCtrT ) ));
   imagesc(TPoint, G_VPoint(1).v_inv, cross.GroupVImg1(1:G_VPoint(1).Num,1:cross.NumCtrT ),[minamp,1]);
   set(gca,'ydir','normal');
   ylabel('Group Velocity (km/s)','FontName','Times New Roman','Fontweight','bold');
   hold on
   plot(TPoint,group.ArrV1(1:cross.NumCtrT), 'g');
   
   % the farer station
   set(gcf,'CurrentAxes',h4);
   hold(h4,'off');
   % imagesc(TPoint(1:cross.NumCtrT),Ypoint,cross.group2(1:cross.PointNum,1:cross.NumCtrT),[0,1]);    
   minamp = min(min(cross.GroupVImg2(1:G_VPoint(2).Num,1:cross.NumCtrT ) ));
   imagesc(TPoint, G_VPoint(2).v_inv, cross.GroupVImg2(1:G_VPoint(2).Num,1:cross.NumCtrT ),[minamp,1]);
   set(gca,'ydir','normal')
   ylabel('Group Velocity (km/s)','FontName','Times New Roman','Fontweight','bold');
   xlabel('Period (s)','FontName','Times New Roman','Fontweight','bold','VerticalAlignment','middle');
   hold on
   plot(TPoint,group.ArrV2(1:cross.NumCtrT), 'g');
end   
   

end    


% --- Calculate envelope image, i.e., to obtain envelope at each T 
function EnvelopeImage = EnvelopeImageCalculation(WinWave, fs, TPoint, StaDist)
%%% new code for group velocity analysis using frequency domain Gaussian filter
%%% crossgroup1 = EnvelopeImageCalculation(cross.wave1(1:cross.PointNum), filter.SampleF, Tpoint, evestadist1);

alfa = [0 100 250 500 1000 2000 4000 20000; 5  8  12  20  25  35  50 75];

guassalfa = interp1(alfa(1,:), alfa(2, :), StaDist);

NumCtrT = length(TPoint);
PtNum = length(WinWave);

nfft = 2^nextpow2(max(PtNum,1024*fs));
xxfft = fft(WinWave, nfft);
fxx = (0:(nfft/2))/nfft*fs; 
IIf = 1:(nfft/2+1);
JJf = (nfft/2+2):nfft;

EnvelopeImage = zeros(NumCtrT, PtNum);
for i = 1:NumCtrT
    CtrT = TPoint(i);
    fc = 1/CtrT;         
    Hf = exp(-guassalfa*(fxx - fc).^2/fc^2);
    yyfft = zeros(1,nfft);
    yyfft(IIf) = xxfft(IIf).*Hf;
    yyfft(JJf) = conj(yyfft((nfft/2):-1:2));
    yy = real(ifft(yyfft, nfft));
    filtwave = abs(hilbert(yy(1:nfft)));
    EnvelopeImage(i, 1:PtNum) = filtwave(1:PtNum);
end


% --- Step 2. Recalulate group velicity and get cross-correlation spectrum.
%  MFT, get group image and group arrival
%  Moving window the waves with respect to the group arrival and do
%  cross-correlation at each narrow period band，
function CrossCorrelation(hObject, handles)
% hObject    handle to CrossCorrelation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global sta filter cross group DeltaTInitial ClickPoint
global TravPtV G_VPoint
h3 = handles.axes3; 
h4 = handles.axes4;
TPoint = cross.StartT:cross.DeltaT:cross.EndT;
StaDistance = abs(sta(1).GCDkm - sta(2).GCDkm);

if cross.WaveType == 1  
   cross.VMin = 2.5;
   cross.VMax = 5.2;
   cross.DeltaV = 0.001;
elseif cross.WaveType == 2 
   cross.VMin = 2.8;
   cross.VMax = 5.8;
   cross.DeltaV = 0.001;
end
% cross.HalfBand = 125; % window half band with respect to the group arrival
cross.StartNum = 1;
cross.EndNum = cross.PointNum;

MaxTravT = StaDistance/cross.VMin;
MinTravT = StaDistance/cross.VMax;
%DeltaTInitial=sta(2).GCDkm/GVWinR(2)-sta(1).GCDkm/GVWinR(2)
MaxShiftNum = floor((MaxTravT - DeltaTInitial)/filter.SampleT+1);
MinShiftNum = floor((MinTravT - DeltaTInitial)/filter.SampleT-1);
ShiftNum = MaxShiftNum - MinShiftNum + 1;
CrossCorrImg = zeros(ShiftNum, cross.NumCtrT);

ShiftPtV = zeros(1,ShiftNum);
ShiftPtT = zeros(1,ShiftNum);

for i = 1:ShiftNum
    ShiftPtT(i) = (MinShiftNum + i - 1)*filter.SampleT + DeltaTInitial;
    ShiftPtV(i) = StaDistance/ShiftPtT(i);
end



for numt = 1:cross.NumCtrT 
    filter.CtrT = cross.StartT + (numt - 1)*cross.DeltaT;
    filter.CtrF = (2/filter.SampleF)/filter.CtrT;
    filter.Length = round(max(512*filter.SampleF, round(5*filter.CtrT*filter.SampleF)));
    HalfFilterNum =  floor(filter.Length/2);
%% --- Narrow-band filtering by moving window   
    switch filter.Domain
        case 1
            filter.LowF = (2/filter.SampleF)/(filter.CtrT + 0.5*filter.BandWidth);
            filter.HighF = (2/filter.SampleF)/(filter.CtrT - 0.5*filter.BandWidth);
        case 2
            filter.LowF = (2/filter.SampleF)*(filter.CtrF - 0.5*filter.BandWidth);
            filter.HighF = (2/filter.SampleF)*(filter.CtrF + 0.5*filter.BandWidth);
    end
    
    if filter.Window == 1   % Kaiser Window
        filterData = fir1(filter.Length, [filter.LowF, filter.HighF], kaiser(filter.Length + 1,filter.KaiserPara));
  
    elseif filter.Window == 2   % Gaussian Window
        filterData = fir1(filter.Length, [filter.LowF, filter.HighF], gausswin(filter.Length + 1,filter.GaussAlfa));
    end
    
    % Moving Window the waves 
    MovingWin1 = ones(1, cross.PointNum + HalfFilterNum);
    MovingWin2 = ones(1, cross.PointNum + HalfFilterNum);

% moving window type - 1 : Gaussian=
%     alfa = 100;
%     MovingWin1(PtNum) = exp(-0.5*((PtNum - group.ArrPt1(numt))/sta(1).SampleT).^2/(alfa^2));
%     MovingWin2(PtNum) = exp(-0.5*((PtNum - group.ArrPt2(numt))/sta(2).SampleT).^2/(alfa^2));

% moving window type - 2 : boxcar + Gaussian taper   
    cross.HalfBand = round(2*filter.CtrT*filter.SampleF);
    WinLowerPt1 = round(group.ArrPt1(numt) - cross.HalfBand*filter.SampleF) - cross.StartNum + 1;
    WinLowerPt2 = round(group.ArrPt2(numt) - cross.HalfBand*filter.SampleF) - cross.StartNum + 1;
    WinUpperPt1 = round(group.ArrPt1(numt) + cross.HalfBand*filter.SampleF) - cross.StartNum + 1;
    WinUpperPt2 = round(group.ArrPt2(numt) + cross.HalfBand*filter.SampleF) - cross.StartNum + 1;
       
    alfa = 50*filter.SampleF; % in points
    if WinLowerPt1 > 1
        MovingWin1(1:WinLowerPt1) = exp(-(WinLowerPt1:-1:1).^2/(alfa^2));
    end
    if WinLowerPt2 > 1
        MovingWin2(1:WinLowerPt2) = exp(-(WinLowerPt2:-1:1).^2/(alfa^2));
    end
        
    if WinUpperPt1 <= cross.PointNum
        MovingWin1(WinUpperPt1:cross.PointNum) = exp(-((WinUpperPt1:cross.PointNum) - WinUpperPt1).^2/(alfa^2));
        MovingWin1((cross.PointNum+1):(cross.PointNum + HalfFilterNum)) = 0;
    else
        MovingWin1(WinUpperPt1:(cross.PointNum + HalfFilterNum)) = 0;
    end
   
    if WinUpperPt2 <= cross.PointNum
        MovingWin2(WinUpperPt2:cross.PointNum) = exp(-((WinUpperPt2:cross.PointNum) - WinUpperPt2).^2/(alfa^2));
        MovingWin2((cross.PointNum+1):(cross.PointNum + HalfFilterNum)) = 0;
    else
        MovingWin2(WinUpperPt2:(cross.PointNum + HalfFilterNum)) = 0;
    end   

    WinWave1 = zeros(1, cross.PointNum + HalfFilterNum);
    WinWave2 = zeros(1, cross.PointNum + HalfFilterNum);

    WinWave1 = cross.wave1(1:(cross.PointNum + HalfFilterNum)).*MovingWin1;
    WinWave2 = cross.wave2(1:(cross.PointNum + HalfFilterNum)).*MovingWin2;
 
    FilteredWave1 = fftfilt(filterData, WinWave1(1:(cross.PointNum + HalfFilterNum)));
    FilteredWave1 = FilteredWave1((cross.PointNum + HalfFilterNum):-1:1);
    FilteredWave1 = fftfilt(filterData, FilteredWave1(1:(cross.PointNum + HalfFilterNum)));
    FilteredWave1 = FilteredWave1((cross.PointNum + HalfFilterNum):-1:1);

    FilteredWave1(1:cross.PointNum) = FilteredWave1(1:cross.PointNum)/max(FilteredWave1(1:cross.PointNum));
	wavehilbert(1:cross.PointNum) = hilbert(FilteredWave1(1:cross.PointNum));
	cross.group1(1:cross.PointNum, numt) = abs(wavehilbert(1:cross.PointNum));
    cross.group1(1:cross.PointNum, numt) = cross.group1(1:cross.PointNum, numt)/max(cross.group1(1:cross.PointNum, numt));
    
    FilteredWave2 = fftfilt(filterData,WinWave2(1:(cross.PointNum + HalfFilterNum)));
    FilteredWave2 = FilteredWave2((cross.PointNum + HalfFilterNum):-1:1);
    FilteredWave2 = fftfilt(filterData, FilteredWave2(1:(cross.PointNum + HalfFilterNum)));
    FilteredWave2 = FilteredWave2((cross.PointNum + HalfFilterNum):-1:1);

    FilteredWave2(1:cross.PointNum) = FilteredWave2(1:cross.PointNum)/max(FilteredWave2(1:cross.PointNum));
	wavehilbert(1:cross.PointNum) = hilbert(FilteredWave2(1:cross.PointNum));
	cross.group2(1:cross.PointNum, numt) = abs(wavehilbert(1:cross.PointNum));
    cross.group2(1:cross.PointNum, numt) = cross.group2(1:cross.PointNum, numt)/max(cross.group2(1:cross.PointNum, numt));
%% ---  Cross-correlation
       
    for k = 1:ShiftNum
        shift = MinShiftNum + k - 1;
        if shift >= 0

            CrossCorrImg(k,numt) = dot(FilteredWave2((1 + shift):cross.PointNum), FilteredWave1(1:(cross.PointNum - shift)));
        else
            CrossCorrImg(k,numt) = dot(FilteredWave2(1:(cross.PointNum + shift)), FilteredWave1((1 - shift):cross.PointNum));
        end
    end
    
    CrossMaxAmp = max(CrossCorrImg(1:ShiftNum,numt));
    if CrossMaxAmp > 0
        CrossCorrImg(1:ShiftNum,numt) = CrossCorrImg(1:ShiftNum,numt)/CrossMaxAmp;
    end
    clear FilteredWave1 FilteredWave2 wavehilbert WinWave1 WinWave2 MovingWin1 MovingWin2 filterData
end

%% --- Plot group velocity spectrum 
if get(handles.Silent_cycle_mode,'Value')~=1
   for numt = 1:cross.NumCtrT   
      cross.GroupVImg1(1:G_VPoint(1).Num, numt) = interp1( TravPtV(1).v_inv ,cross.group1(1:cross.PointNum,numt),G_VPoint(1).v );
      cross.GroupVImg2(1:G_VPoint(2).Num, numt) = interp1( TravPtV(2).v_inv ,cross.group2(1:cross.PointNum,numt),G_VPoint(2).v );
   end 
   group.ArrPoint1 = zeros(1, cross.NumCtrT);
   group.ArrPoint2 = zeros(1, cross.NumCtrT);
 
   [MaxAmp,group.ArrPoint1(1:cross.NumCtrT)] = max(cross.GroupVImg1(1:G_VPoint(1).Num,1:cross.NumCtrT));
   [MaxAmp,group.ArrPoint2(1:cross.NumCtrT)] = max(cross.GroupVImg2(1:G_VPoint(2).Num,1:cross.NumCtrT)); 
   group.ArrV1(1:cross.NumCtrT)= G_VPoint(1).v( group.ArrPoint1(1:cross.NumCtrT) );
   group.ArrV2(1:cross.NumCtrT)= G_VPoint(2).v( group.ArrPoint2(1:cross.NumCtrT) );    
   % The nearer station, 1 is the nearer station
   set(gcf,'CurrentAxes',h3);
   hold(h3,'off');
   minamp = min(min(cross.GroupVImg1(1:G_VPoint(1).Num,1:cross.NumCtrT ) ));
   imagesc(TPoint, G_VPoint(1).v, cross.GroupVImg1(1:G_VPoint(1).Num,1:cross.NumCtrT ),[minamp,1]);
   set(gca,'ydir','normal');
   ylabel('Group Velocity (km/s)','FontName','Times New Roman','Fontweight','bold');
   hold(h3,'on');
   plot(TPoint,group.ArrV1, 'g');
   % the farer station
   set(gcf,'CurrentAxes',h4);
   hold(h4,'off');
   % imagesc(TPoint(1:cross.NumCtrT),Ypoint,cross.group2(1:cross.PointNum,1:cross.NumCtrT),[0,1]);    
   minamp = min(min(cross.GroupVImg2(1:G_VPoint(2).Num,1:cross.NumCtrT ) ));
   imagesc(TPoint, G_VPoint(2).v, cross.GroupVImg2(1:G_VPoint(2).Num,1:cross.NumCtrT ),[minamp,1]);
   set(gca,'ydir','normal')
   ylabel('Group Velocity (km/s)','FontName','Times New Roman','Fontweight','bold');
   xlabel('Period (s)','FontName','Times New Roman','Fontweight','bold','VerticalAlignment','middle');
   hold(h4,'on');
   plot(TPoint,group.ArrV2, 'g');
end % silent
%% -- For dispertion search: cross.PhasVImg
% VPoint,ShiftPtV the revse order 
% ShiftPtV(i) = StaDistance/ShiftPtT(i);[MinTravT:1:MaxTravT]

VPoint = cross.VMin:cross.DeltaV:cross.VMax;
CImgPt = size(VPoint,2);
cross.PhasVImg = zeros(CImgPt, cross.NumCtrT); 
for i = 1:cross.NumCtrT
   cross.PhasVImg(1:CImgPt, i) = interp1(ShiftPtV, CrossCorrImg(1:ShiftNum,i), VPoint, 'spline');
   MaxAmpPhaseVImg = max(cross.PhasVImg(1:CImgPt, i));
   if MaxAmpPhaseVImg > 0
       cross.PhasVImg(1:CImgPt, i) = cross.PhasVImg(1:CImgPt, i)/MaxAmpPhaseVImg;
   end       
end


% --- Step 3.A Seim-automatic Pickup Disper Module    
% --- Semiautomatic pick
function IsDispGood = PhaseVDisper(hObject, handles)
global cross ClickPoint sta src PrevStaLonLat
global PhaseDispDirectory
global g_refphasedisp
h1 = handles.axes1;
h2 = handles.axes2;

SaveOrNot = questdlg('Do you want to get dispersion curve?','Dispersion Curve','Yes','No','Yes');
axes(h2);

if strcmp(SaveOrNot, 'Yes') == 1%  
%     k = 1;
%     ClickPoint = zeros(2,3);
%     while k ~= 0
%         set(handles.MsgEdit,'String','Please left click on lower figure to select disperion curve!');
%         k = waitforbuttonpress;
%         x = get(h2,'XLim');
%         y = get(h2,'YLim');
%         if ClickPoint(1,1) >= x(1,1) && ClickPoint(1,1) <= x(1,2) && ClickPoint(1,2) >= y(1,1) && ClickPoint(1,2) <= y(1,2)
%             hold on
%             plot(h2, ClickPoint(1,1), ClickPoint(1,2),'*g');
%             k = 0;
%             set(handles.MsgEdit, 'String', 'Successful click!');
%         else      
%             k = 1;
%             set(handles.MsgEdit, 'String', 'Wrong click! Please click on the lower figure!');
%         end
%     end
   k = 1;
   while k ~= 0
        set(handles.MsgEdit,'String','Left button click for first point and right for last point to select disperion curve!');
        % Initially, the list of points is empty.
        xy = [];
        n = 0;
        % Loop, picking up the points.
%         disp('Left mouse button picks points.')
%         disp('Right mouse button picks last point.')
        button = 1;
        x = get(h2,'XLim');
        y = get(h2,'YLim');        
        while button == 1
            [xi,yi,button] = ginput(1);
            if xi > x(2) || xi < x(1) || yi > y(2) || yi < y(1)
                set(handles.MsgEdit, 'String', 'Wrong click! Please click on the lower image!');
            else
                hold on; plot(xi,yi,'g+','MarkerSize',8)
                n = n+1;
                xy(:,n) = [xi;yi];  
            end
%             if n == 1
%                 break
%             end
        end
        hold on; plot(xi,yi,'m+','MarkerSize',8)

        if n == 0
            k = 1;
            set(handles.MsgEdit, 'String', 'Wrong click! Please click on the lower image!');
        else
            k = 0;
        end
    end
    
    
    DisperStartT = cross.StartT;
    DisperEndT = cross.EndT;

    StaDistance = abs(sta(1).GCDkm - sta(2).GCDkm);
    TPoint = cross.StartT:cross.DeltaT:cross.EndT;
%     InitialT = floor((ClickPoint(1,1) - cross.StartT)/cross.DeltaT + 1);
%     InitialV = ClickPoint(1,2);
%     InitialY = floor((InitialV - cross.VMin)/cross.DeltaV + 1);
    InitialT = round((xy(1,:) - cross.StartT)/cross.DeltaT + 1);
    InitialY = round((xy(2,:) - cross.VMin)/cross.DeltaV + 1);
    nSelectPt = length(InitialT);
    % sort picked points according to increasing periods
    [InitialT, II] = sort(InitialT);
    InitialY = InitialY(II);
    
    DispPt = zeros(1, cross.NumCtrT); 
    PhaseVDisp = zeros(1, cross.NumCtrT); 
    CImgPt = (cross.VMax - cross.VMin)/cross.DeltaV + 1;
    % DispPt(1:cross.NumCtrT) = AutoSearch(InitialY, InitialT, cross.PhasVImg(1:CImgPt,1:cross.NumCtrT));
    if nSelectPt == 1
        DispPt(1:cross.NumCtrT) = AutoSearch(InitialY, InitialT, cross.PhasVImg(1:CImgPt,1:cross.NumCtrT));
    elseif nSelectPt > 1
        DispPt(1:cross.NumCtrT) = AutoSearchMultiplePoints(InitialY, InitialT, cross.PhasVImg(1:CImgPt,1:cross.NumCtrT));
    end       
    
    PhaseVDisp(1:cross.NumCtrT) = cross.VMin + (DispPt - 1)*cross.DeltaV;
    hold(h2,'on')
    plot(h2, TPoint(1:cross.NumCtrT), PhaseVDisp(1:cross.NumCtrT), 'k-', 'LineWidth', 3);
    % plot the reference
    if get(handles.RefC_Disper,'Value')==1
       refc_StartNum=find( TPoint(1) == g_refphasedisp(:,1) );
       refc_EndNum=find( TPoint(end) == g_refphasedisp(:,1) );
       refc_TNum=refc_StartNum:refc_EndNum;
       hold on; plot(TPoint, g_refphasedisp(refc_TNum,2),'r--');
       hold on; plot(TPoint, g_refphasedisp(refc_TNum,3),'r--');
    end
    
    IsDispGood = 2;

    InputOrNotIndex = 1;
    while InputOrNotIndex == 1  
      set(handles.MsgEdit,'String','Left button click for start period point and end period point to select disperion curve!');    
      for loni=1:2
          [xpoint(loni),ypoint(loni),button] = ginput(1);
          plot([xpoint(loni),xpoint(loni)],[ypoint(loni)-2,ypoint(loni)+2],'c-');
      end
      DisperStartT = xpoint(1);
      DisperEndT =  xpoint(2);
      if isempty(DisperStartT) || isempty(DisperEndT)
         InputOrNotIndex = 1;
         questdlg('Wrong click! Please click on the right place','OK');
      elseif DisperStartT >= cross.StartT && DisperEndT <= cross.EndT && DisperStartT<DisperEndT
         InputOrNotIndex = 0;
         StartTIndex = round((DisperStartT - cross.StartT)/cross.DeltaT) + 1;
         EndTIndex = round((DisperEndT - cross.StartT)/cross.DeltaT) + 1;    
      else
         InputOrNotIndex = 1;
         questdlg('Wrong click! Please click on the right place','OK');
      end
    end
%     %key in the period you wantto save
%     while InputOrNotIndex == 1    
%         prompt = {'Enter Start Period:','Enter End Period:'};
%         title = ['Set start and end period (s) for saving dispersion data'];
%         line = 2; 
%         def = {num2str(cross.StartT), num2str(cross.EndT)};
%         DisperPeriod = inputdlg(prompt,title,line, def);
%         if size(DisperPeriod,1)~=line
%             InputOrNotIndex = 1;
%         else
%             DisperStartT = str2num(DisperPeriod{1});
%             DisperEndT = str2num(DisperPeriod{2});
%             if isempty(DisperStartT) || isempty(DisperEndT)
%                 InputOrNotIndex = 1;
%             else
%                 StartTIndex = round((DisperStartT - cross.StartT)/cross.DeltaT) + 1;
%                 EndTIndex = round((DisperEndT - cross.StartT)/cross.DeltaT) + 1;    
%                 if DisperStartT >= cross.StartT && DisperEndT <= cross.EndT
%                     InputOrNotIndex = 0;
%                 else
%                     InputOrNotIndex = 1;
%                 end
%             end
%         end
%     end
       

    %write T-V to file 
    if PhaseDispDirectory == 0
        PhaseDispDirectory = pwd;
    end
    
    YYChar = num2str(src.YY);
    
    MonChar = num2str(src.Month);
    if length(MonChar) == 1
        MonChar = strcat('0',MonChar);
    end
    
    DDChar = num2str(src.Day);
    if length(DDChar) == 1
        DDChar = strcat('0',DDChar);
    end   
    
    HHChar = num2str(src.HH);
    if length(HHChar) == 1
        HHChar = strcat('0',HHChar);
    end
    
    MMChar = num2str(src.MM);
    if length(MMChar) == 1
        MMChar = strcat('0',MMChar);
    end   
    
    
    SrcTime = strcat(YYChar(3:4),MonChar,DDChar,'_',HHChar,MMChar);
    
    for i = 1:max(size(sta(1).Name, 2), size(sta(2).Name, 2))
        if double(sta(1).Name(i)) > double(sta(2).Name(i))
            StaNamePair = strcat(sta(2).Name,'-', sta(1).Name);
            index1 = 2;
            index2 = 1;
            break
        elseif double(sta(1).Name(i)) < double(sta(2).Name(i))
            StaNamePair = strcat(sta(1).Name,'-', sta(2).Name);
            index1 = 1;
            index2 = 2;
            break
        else
            continue
        end
    end
    
    if sum(sum(PrevStaLonLat)) ~= 0
        set(gcf,'CurrentAxes',h1);
        hold(h1,'on');
        plot(PrevStaLonLat(1:2,1), PrevStaLonLat(1:2,2), 'b-');
        plot(h1, PrevStaLonLat(1:2,1), PrevStaLonLat(1:2,2), 'r^', 'MarkerSize',6, 'MarkerFaceColor','r');
    end

    
    PrevStaLonLat = [sta(index1).Lon, sta(index1).Lat;sta(index2).Lon, sta(index2).Lat ];

    
    dispfilename = strcat('C' ,SrcTime, StaNamePair, '.dat');
    TVfile = fullfile(PhaseDispDirectory, dispfilename);
    
    ftv = fopen(TVfile,'w');
    fprintf(ftv,'%f     ', sta(index1).Lon);
    fprintf(ftv,'%f\n', sta(index1).Lat);
    fprintf(ftv,'%f     ', sta(index2).Lon);
    fprintf(ftv,'%f\n', sta(index2).Lat);
    DataExistIndex = 0;
    for i = StartTIndex:EndTIndex
       wavelength = PhaseVDisp(i)*TPoint(i);
       if wavelength <= 2*StaDistance
           DataExistIndex = 1;
           fprintf(ftv,'%4.1f   ',TPoint(i));
           fprintf(ftv,'%4.3f\n',PhaseVDisp(i));
           set(gcf,'CurrentAxes',h2);
           hold(h2,'on');
           plot(h2, TPoint(i), PhaseVDisp(i), 'co','MarkerSize',4, 'MarkerFaceColor','c');
       else
           break
       end
    end
    fclose(ftv);
    
    
    ReviseOrNot = questdlg('Want to revise the saved dispersion?','Revise dispersion data','Yes','No','No');
    if strcmp(ReviseOrNot, 'Yes') == 1
        IsDispGood = 1;
    else
        if DataExistIndex ~= 1
            delete(TVfile)
            set(handles.MsgEdit, 'String', 'No Disperion Data Written! Disper File Deleted!');
        else
            set(gcf,'CurrentAxes',h1);
            hold(h1,'on');
            plot([sta(index1).Lon, sta(index2).Lon], [sta(index1).Lat, sta(index2).Lat], 'r-');
            plot(h1, [sta(index1).Lon, sta(index2).Lon], [sta(index1).Lat, sta(index2).Lat], 'k^', 'MarkerSize',6, 'MarkerFaceColor','g');
        end
    end

else
    IsDispGood = 2;
end

% --- Automatically search arrival time line on a image 
% --- only pick outone point
function ArrPt = AutoSearch(InitialY, InitialX, ImageData)
% Input: InitialY  InitlaX   ImageData
% OutPut: ArrPt

YSize = size(ImageData, 1);
XSize = size(ImageData, 2);

ArrPt = zeros(1,XSize);

% Center_T search up
step = 3;
point_left = 0;
point_right = 0;

for i = InitialX:XSize
	index1 = 0;
	index2 = 0; 
	point_left = InitialY;
	point_right = InitialY;
	while index1 == 0
        point_left_new = max(1, point_left - step);
	    if ImageData(point_left,i) < ImageData(point_left_new,i)
		  point_left = point_left_new;
%          point_right = point_right - step;
		else
		   index1 = 1;
		   point_left = point_left_new;
        end
    end
    while index2 == 0
        point_right_new = min(point_right + step, YSize);
		if ImageData(point_right,i) < ImageData(point_right_new,i)
		   point_right = point_right_new;
%           point_left = point_left + step;
		else
           index2=1;
		   point_right = point_right_new;
		end
    end

    [MaxAmp, index_max] = max(ImageData(point_left:point_right,i));
    ArrPt(i) = index_max + point_left - 1;
    InitialY = ArrPt(i);
        
end  %end for

% Center_T search down

InitialY = ArrPt(InitialX);
for i = (InitialX - 1):(-1):1
	index1 = 0;
	index2 = 0; 
    point_left = InitialY;
	point_right = InitialY;

	while index1 == 0
        point_left_new = max(1, point_left - step);
	    if ImageData(point_left,i) < ImageData(point_left_new,i)
		  point_left = point_left_new;
%          point_right = point_right - step;
		else
		   index1 = 1;
		   point_left = point_left_new;
        end
    end
    while index2 == 0
        point_right_new = min(point_right + step, YSize);
		if ImageData(point_right,i) < ImageData(point_right_new,i)
		   point_right = point_right_new;
%           point_left = point_left + step;
		else
           index2=1;
		   point_right = point_right_new;
		end
    end
    
    [MaxAmp, index_max] = max(ImageData(point_left:point_right,i));
    ArrPt(i) = index_max + point_left - 1;
    InitialY = ArrPt(i);
    
end  %end for

% --- Automatically search arrival time line on a image with multiple
function ArrPt = AutoSearchMultiplePoints(ptY, ptX, ImageData)
% Input: InitialY  InitlaX ImageData
% OutPut: ArrPt
% DispPt(1:cross.NumCtrT) = AutoSearch(InitialY, InitialT, cross.PhasVImg(1:CImgPt,1:cross.NumCtrT));
% sort ptX and ptY according to the increasing of ptX
[ptX, II] = sort(ptX);
ptY = ptY(II);

YSize = size(ImageData, 1);
XSize = size(ImageData, 2);
nPt = length(ptX);  % number of input searching points

ArrPt = zeros(1,XSize);

% X searching up for the point with maximum X
step = 3;
point_left = 0;
point_right = 0;

InitialX = ptX(nPt);
InitialY = ptY(nPt);
for i = ptX(nPt):XSize
	index1 = 0;
	index2 = 0; 
	point_left = InitialY;
	point_right = InitialY;
	while index1 == 0
        point_left_new = max(1, point_left - step);
	    if ImageData(point_left,i) < ImageData(point_left_new,i)
		  point_left = point_left_new;
%          point_right = point_right - step;
		else
		   index1 = 1;
		   point_left = point_left_new;
        end
    end
    while index2 == 0
        point_right_new = min(point_right + step, YSize);
		if ImageData(point_right,i) < ImageData(point_right_new,i)
		   point_right = point_right_new;
%           point_left = point_left + step;
		else
           index2=1;
		   point_right = point_right_new;
		end
    end

    [MaxAmp, index_max] = max(ImageData(point_left:point_right,i));
    ArrPt(i) = index_max + point_left - 1;
    InitialY = ArrPt(i);
        
end  %end for

% X searching down for the point with maximum X. There will other points
% with smaller X which will act as internal constraints for the searching
% process

InitialX = ptX(nPt);
InitialY = ArrPt(ptX(nPt));
midX = ptX(nPt-1);
midY = ptY(nPt-1);
kk = 0;
for i = ptX(nPt):(-1):1
	index1 = 0;
	index2 = 0;
    
    if i == midX
        InitialY = midY;
        kk = kk + 1;
        if (nPt - kk) > 1
            midX = ptX(nPt - kk - 1);
            midY = ptY(nPt - kk - 1);
        end
    end
            
    point_left = InitialY;
	point_right = InitialY;

	while index1 == 0
        point_left_new = max(1, point_left - step);
	    if ImageData(point_left,i) < ImageData(point_left_new,i)
		  point_left = point_left_new;
%          point_right = point_right - step;
		else
		   index1 = 1;
		   point_left = point_left_new;
        end
    end
    while index2 == 0
        point_right_new = min(point_right + step, YSize);
		if ImageData(point_right,i) < ImageData(point_right_new,i)
		   point_right = point_right_new;
%           point_left = point_left + step;
		else
           index2=1;
		   point_right = point_right_new;
		end
    end
    
    [MaxAmp, index_max] = max(ImageData(point_left:point_right,i));
    ArrPt(i) = index_max + point_left - 1;
    InitialY = ArrPt(i);
    
end  %end for


% --- Step 3.B Automatic Pickup Disper Module   
% --- Automatic pick
function IsDispGood = AutoPhaseDisperProcessing(hObject, handles)
% The use of product of function crossculation: cross.PhasVImg
global g_refphasedisp sta cross

TPoint = cross.StartT:cross.DeltaT:cross.EndT;
VPoint = cross.VMin:cross.DeltaV:cross.VMax;
% PhaseVImg(1:VImgPt,1:NumCtrT)
VImgPt = length(VPoint); 
NumCtrT = length(TPoint); 
dc = VPoint(2) - VPoint(1);
dT = TPoint(2) - TPoint(1);


cRef_low = interp1(g_refphasedisp(:,1), g_refphasedisp(:,2), TPoint, 'pchip');
cRef_high = interp1(g_refphasedisp(:,1), g_refphasedisp(:,3), TPoint, 'pchip');
cRef = (cRef_low + cRef_high)/2;

% half-wave  criterion 
% find the approximate maximum period for dispersion analysis which
% satisfies interstation-distance > minlamdaRatio*c*T
StaDistance = abs(sta(1).GCDkm - sta(2).GCDkm);
minlamdaRatio =2;
% minlamdaRatio =str2num(get(handles.minWavelength, 'String'));
lamda = cRef.*TPoint;
II = find((lamda*minlamdaRatio) >= StaDistance );
if length(II) > 0
    nMaxT = II(1);
else
    nMaxT = NumCtrT;
end

%% -- step 1. Search the dispersion curves from the specifici (period,velocity) points
% set tial T and c for searching dispersion curves
T_search=[1 5]/10;
cmax_ref = max(cRef_high);
if cmax_ref>cross.VMax
   cmax_ref=cross.VMax;
end
cmin_ref = min(cRef_low); 
% choose the search staring period and the number of it is experimental
T_try_index = round(T_search*(TPoint(nMaxT) - TPoint(1))/dT)+1;  
T_try = TPoint(1) + (T_try_index - 1)*dT;  % trial T for searching dispersion index
c_try_index = round(((cmin_ref - VPoint(1)) + (cmax_ref - cmin_ref)*(1:9)/10)/dc) + 1;
c_try = VPoint(1) + (c_try_index - 1)*dc;

% search disperison curves
PhaseVDisp_try = zeros(length(T_try)*length(c_try), NumCtrT);
k = 0;
for i = 1:length(T_try)
% hold on; plot(T_try(i)*ones(1, length(c_try)), c_try, 'w*');
    for j = 1:length(c_try)
        Initialc = c_try_index(j);
        InitialT = T_try_index(i);
        DispPt = AutoSearch(Initialc, InitialT, cross.PhasVImg);%using the relative position
        
        if k == 0
            k =  k + 1;
            PhaseVDisp_try(k,:) = VPoint(1) + (DispPt - 1)*dc;
            % hold on; plot(TPoint, PhaseVDisp_try(k,:), 'g');
        else
            tempDisp = VPoint(1) + (DispPt - 1)*dc;
            SameDisperIndex = 0;
            for nn = 1:k
                if sum(abs(PhaseVDisp_try(nn,:) - tempDisp)) < 1.0e-4
                    SameDisperIndex = 1;
                end
            end
            if SameDisperIndex == 0
                k = k + 1;
                PhaseVDisp_try(k,:) = tempDisp;
%               hold on; plot(TPoint, PhaseVDisp_try(k,:), 'g');
            end
        end                  
    end
end

%% -- step 2. Pick up the only one dispersion curve from the autosearching
% determine the quality of the dispersion curve by looking at how many
% disperion points fall in the reference dispersion range so the reference
% is very important for dispersion picking

NumDispCurve = k;% number of the disperison curves
InRangePt = zeros(1,NumDispCurve); % number of the poinits falling into the ref
for i = 1:NumDispCurve
    GoodIndex = sign(PhaseVDisp_try(i,1:nMaxT) - cRef_low(1:nMaxT)) + sign(cRef_high(1:nMaxT) - PhaseVDisp_try(i,1:nMaxT));
    II = find(GoodIndex == 2);
    InRangePt(i) = length(II); 
end
maxpt = max(InRangePt);
meanpt = mean(InRangePt);
% find the best several dispersion curves within reference range
II = find(InRangePt >= (2*maxpt+meanpt)/3); 
if length(II) == 1   
   PhaseVDisp = PhaseVDisp_try(II,:);

else
   RefObsDispDiff = zeros(1, length(II));
   ObsSumAbsDiff = zeros(1, length(II));
   for i = 1:length(II)
       RefObsDispDiff(i) = sum(abs(PhaseVDisp_try(II(i),1:nMaxT) - (cRef_low(1:nMaxT) + cRef_high(1:nMaxT))/2));% lowest difference dispersion curve
       ObsSumAbsDiff(i) = sum(abs(diff(PhaseVDisp_try(II(i),1:nMaxT))));% smoothest dispersion
   end
   [mindiff, index1] = min(RefObsDispDiff); % find lowest difference dispersion curve with respect to reference
   [minabs, index2] = min(ObsSumAbsDiff); % find smoothest dispersion curve
    
   if index1 == index2
       PhaseVDisp = PhaseVDisp_try(II(index1),:);
   else
      BestTwoDiff = abs(PhaseVDisp_try(II(index1),1:nMaxT) - PhaseVDisp_try(II(index2),1:nMaxT));
      if length(find(BestTwoDiff < 1.0e-3)) > 0.67*nMaxT  % 2/3 of two best dispersion curves are overlapping
         PhaseVDisp = PhaseVDisp_try(II(index2),:); % choose the smoothest one
      else
         PhaseVDisp = PhaseVDisp_try(II(index1),:); % choose the smaller difference one if 2/3 dispersion are different
      end
   end        
end
%%% plot the dispersion curve in step two
if get(handles.Silent_cycle_mode,'Value')==0
   hold on; 
% 1st,plot the final selected dispersion curve
   plot(TPoint, PhaseVDisp,'g', 'LineWidth', 1.5); 
   refc_StartNum=find( TPoint(1) == g_refphasedisp(:,1) );
   refc_EndNum=find( TPoint(end) == g_refphasedisp(:,1) );
   refc_TNum=refc_StartNum:refc_EndNum;
   plot(TPoint, g_refphasedisp(refc_TNum,2),'r--');
   plot(TPoint, g_refphasedisp(refc_TNum,3),'r--');
end

%% -- Step 3. extracts the reasonable part of a whole disper curve between two stations.

RawDisper = [TPoint' PhaseVDisp'];
MinTWidth = dT*25;   % set minimum length of dispersion seg to be saved
% here: dt=cross.DeltaT
if  cross.WaveType == 2
    NewDisper = AutomaticDisperLove(RawDisper,1/dT,MinTWidth,TPoint(nMaxT));
else                    
    NewDisper = AutomaticDisperRayl(RawDisper,1/dT,MinTWidth, TPoint(nMaxT));
end
if get(handles.Silent_cycle_mode,'Value')==0
   plot(NewDisper(:,1), NewDisper(:,2), 'go','MarkerSize',7,'MarkerFaceColor','g');
end

%% -- Step 4. reference cut

% whether NewDisper falls in PhaseV range or not
% NewDisper = [RawDisper, ones(size(RawDisper,1),1)];
% find reasonable dispersion points;in EFG processing this can be relative
% SRN, so here still need to be improved.
% NewDisper:[T; V; quality]
for ii = 1:size(NewDisper,1)
    if NewDisper(ii,2) > cRef_high(ii) || NewDisper(ii,2) < cRef_low(ii)
       NewDisper(ii,2)=0;
       NewDisper(ii,3)=0;
    end
end
%% -- Step 5. save result
MinGap=2;
NewDisper2=Absorb_DispSeg(NewDisper,1/dT,MinGap);
GoodIndex = NewDisper2(:,3);
II = find(NewDisper2(:,3) == 1);

if length(II) >= MinTWidth
   if get(handles.Silent_cycle_mode,'Value')==0
      plot(NewDisper(II,1), NewDisper(II,2), 'co','MarkerSize',5,'MarkerFaceColor','c');
   end
   % pick dispersion
   % the data:PhaseVDisp is from the 2nd Step;GoodIndex is got from the 3rd
   % step.
   SaveAutoPick_PhaseDisper(hObject, handles, PhaseVDisp, GoodIndex);
   IsDispGood = 2;
else 
   IsDispGood = 0;% no outputs
end

% --- for aborb the out range points
function NewDisper2=Absorb_DispSeg(NewDisper,SampleF,MinGap)

% Set minist acceptable gap 
% MinGap=5;
% NewDisper:[T; V; quality]
% To affirm the segmentation
% Nonzero=[Tstart;Tend;length]
Nonzero=zeros(1,3);
NonzeroIndex=0;
for ii=1:size(NewDisper,1)
    if NewDisper(ii,3)==1 
       if ii==1
          NonzeroIndex=NonzeroIndex+1;
          Nonzero(NonzeroIndex,1)=NewDisper(ii,1);
       elseif NewDisper(ii-1,3)<1  % open up the next
          NonzeroIndex=NonzeroIndex+1;
          Nonzero(NonzeroIndex,1)=NewDisper(ii,1);
       end
       Nonzero(NonzeroIndex,3)=Nonzero(NonzeroIndex,3)+1; 
    end
end
Nonzero(:,2)=Nonzero(:,1)+Nonzero(:,3)/SampleF;
Nonzero=sortrows(Nonzero,1);

% MinSegLength=MinTWidth*SampleF;  % MinSegLength=dT*30*1/dT=30;
% AcceptedSeg=find(Nonzero(:,3)>=MinSegLength);
% SegNum=length(AcceptedSeg);
SegNum=NonzeroIndex;


if SegNum<=1
   NewDisper2=NewDisper;
   return
end

if SegNum>1
   NewDisper2=NewDisper;
   %sort_seq=[starT;endT;length;gap(along with the largest Seg)]
   sort_seq=zeros(SegNum,4);
   for loni=1:SegNum
       sort_seq(loni,1:3)=Nonzero(loni,1:3);
   end
   % sort in the reverse order of the length,so the first is the longest
   % one,and we think it is the stablest 
   sort_seq=sortrows(sort_seq,-3);
   % desgin the gap value(the distance from the first Seg) 
   for loni=2:SegNum
       if sort_seq(loni,2)<sort_seq(1,1) 
          sort_seq(loni,4)=sort_seq(1,1)-sort_seq(loni,2);        
       elseif sort_seq(loni,1)>sort_seq(1,2)
          sort_seq(loni,4)=sort_seq(loni,1)-sort_seq(1,2);
       end
   end
   % then sort in order of the gap value
   sort_seq=sortrows(sort_seq,4);
   PrevTstartIX=find(NewDisper(:,1)==sort_seq(1,1));
   PrevTendIX=find(NewDisper(:,1)==sort_seq(1,2)); 
   % If the gap <= MinGap, the result absorb the segment out of the range
   % of reference.
   for ii=2:SegNum
       TstartIX = find( NewDisper(:,1)==sort_seq(ii,1) );
       TendIX = find( NewDisper(:,1)==sort_seq(ii,2) ); 
       if sort_seq(ii,4)<=MinGap
          TstartIX = find( NewDisper(:,1)==sort_seq(ii,1) );
          TendIX = find( NewDisper(:,1)==sort_seq(ii,2) ); 
          if TendIX<PrevTstartIX 
              NewDisper2(TendIX-1:PrevTstartIX+1,3)=1;
          elseif PrevTendIX<TstartIX
              NewDisper2(PrevTendIX-1:TstartIX+1,3)=1;
          end
          % exchange
          if PrevTstartIX>TstartIX
             PrevTstartIX = TstartIX; 
          end
          if PrevTendIX<TendIX
             PrevTendIX = TendIX; 
          end
       else 
         NewDisper2(TstartIX:TendIX,3)=0;
       end
   end
       
       
end

% --- for Love waves
function NewDisper = AutomaticDisperLove(RawDisper,SampleF,MinTWidth, MaxTCal)
% For Love wave automatic pick
% This function extracts the reasonable part of a whole disper curve between two stations.
% e.g. MinTWidth=5;  % choose segments containing at least 5s data
% Rewrite RawDisper
% RawDisper and NewDispers
% 1     2      3
% T   PhaseV quality
% Range of dcdt: if needed, this part can be put in input parameters.
% 1       2          3     
% T   upper-limit  lower-limit
r_dcdt=[  0  0.30  -0.025
       10  0.10  -0.025
       30  0.06  -0.010
       45  0.03  -0.006
      200  0.02  -0.003];  
    
r_dcdt=r_dcdt/SampleF;% SampleF=1/cross.DeltaT


% Create smooth range of dcdt: 
Smooth_r_dcdt = zeros(size(RawDisper,1)-1,3);
Smooth_r_dcdt(:,1) = RawDisper(1:end-1,1);
Smooth_r_dcdt(:,2) = interp1(r_dcdt(:,1),r_dcdt(:,2),Smooth_r_dcdt(:,1),'pchip');
Smooth_r_dcdt(:,3) = interp1(r_dcdt(:,1),r_dcdt(:,3),Smooth_r_dcdt(:,1),'pchip');

% dcdt: time derivative of phase velocity
% 1     2     3
% T    dcdt  g/b(qualities) 
dcdt=zeros(size(RawDisper,1)-1,3);
dcdt(:,1)=RawDisper(1:end-1,1);
dcdt(:,2)=RawDisper(2:end,2)-RawDisper(1:end-1,2);

%%%%%%%%%%%%%%%%%%%%% Quanlity Definition 1. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if lower-limit <= dcdt < upper-limit : good = 1
% if 2 * lower-limit < dcdt < 1.5 * upper-limit  : half-good =1/2
% other cases : bad = 0
for ii=1:size(dcdt,1)
    upper_limit=Smooth_r_dcdt(ii,2);
    lower_limit=Smooth_r_dcdt(ii,3);
    
    if dcdt(ii,2)>0 && dcdt(ii,2)<=upper_limit
        dcdt(ii,3)=1;
    elseif dcdt(ii,2)>= 1.5 * lower_limit && dcdt(ii,2)<=upper_limit*1.5
        dcdt(ii,3)=0.5;
    end
end
%%%%%%%%%%%%%%%%%%%% Quanlity Definition 2.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Recheck data quality
% if d(dc/dt)/dt changes fast => bad 
for ii = 2:size(dcdt,1)
    ddcdtt = abs(dcdt(ii,2)-dcdt(ii-1,2));
    if ddcdtt > Smooth_r_dcdt(ii,2)-Smooth_r_dcdt(ii,3)
       dcdt(ii,3)=0;
    end
end
%%%%%%%%%%%%%%%%%%%% Quanlity Definition 3.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% if there is less than 2 half-good pts in the middle of good pts, just
% accepted them 
% here had bertter SampleF=1
if SampleF == 1
   for ii=2:size(dcdt,1)-2
       if dcdt(ii,3)==0.5 && dcdt(ii-1,3)==1 && dcdt(ii+1,3)==1
          dcdt(ii,3)=1;
       elseif dcdt(ii,3)==0.5 && dcdt(ii-1,3)==1 && dcdt(ii+1,3)>0 && dcdt(ii+2,3)==1
          dcdt(ii,3)=1;
       end
   end
elseif SampleF == 2
    for ii = 2:size(dcdt,1)-SampleF
        if dcdt(ii,3)==0 && dcdt(ii-1,3)>0 && dcdt(ii+1,3)>0
            dcdt(ii,3)=0.5;
        end
    end
    for ii = 2:size(dcdt,1)-SampleF
        if dcdt(ii,3) == 0.5 && dcdt(ii-1,3)==1
            QaulityAfter=dcdt(ii+1:end,3);
            NearestOne=min(QaulityAfter==1);
            if sum(QaulityAfter(1:NearestOne))/NearestOne > 0.5
                dcdt(ii,3)=1;
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%% Link up the good points %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the length of each good-data segment: dcdt(ii,3)==1
% Nonzero:
%   1      2       3 
% Tstart  Tend   length
Nonzero=zeros(size(dcdt,1),3);
NonzeroIndex=0;
for ii=1:size(dcdt,1)
    if dcdt(ii,3)==1 % check every good quality point 
       if ii==1
          NonzeroIndex=NonzeroIndex+1;
          Nonzero(NonzeroIndex,1)=dcdt(ii,1);
       elseif dcdt(ii-1,3)<1  % open up the next
          NonzeroIndex=NonzeroIndex+1;
          Nonzero(NonzeroIndex,1)=dcdt(ii,1);
       end
       Nonzero(NonzeroIndex,3)=Nonzero(NonzeroIndex,3)+1; % add to the length
    end
end
Nonzero(:,2)=Nonzero(:,1)+Nonzero(:,3)/SampleF;
Nonzero=sortrows(Nonzero,1);

MinSegLength=MinTWidth*SampleF;  % MinSegLength=dT*30*1/dT=30;
AcceptedSeg=find(Nonzero(:,3)>=MinSegLength);
SegNum=length(AcceptedSeg);

%%%%%%%%%%%%%%%%%%%% The process of acception of Seg %%%%%%%%%%%%%%%%%%%%%% 
% Analysis wheater the different AcceptedSeg belong to the Same curve 
% NewDisper the same format as RawDisper
if SegNum<1
   NewDisper=zeros(size(RawDisper,1),size(RawDisper,2));
   NewDisper(:,1) = RawDisper(:,1);
   return
end
% in fact, here is if SegNum=1;,the only one and the longest
if SegNum==1
   TstartIX=find( RawDisper(:,1)==Nonzero(AcceptedSeg(1),1) );
   TendIX=find( RawDisper(:,1)==Nonzero(AcceptedSeg(1),2) );
   NewDisper=RawDisper(TstartIX:TendIX,:);
end

if SegNum>=2
   % begin the judge from the longest Seg
   % there will be two or more maxium but matlab max only return the first
   % index
   [longth,longestSeg]=max(Nonzero(:,3));
   PrevTstartIX=find(RawDisper(:,1)==Nonzero(longestSeg,1));
   PrevTendIX=find(RawDisper(:,1)==Nonzero(longestSeg,2));
   NewDisper=RawDisper(PrevTstartIX:PrevTendIX,:);

   % we built sort_seq to replace Nonzero
   %   1    2     3     4
   % starT endT length gap(with the largest Seg)
   sort_seq=zeros(SegNum,4);
   % assign for sort_seq 
   for loni=1:SegNum
       sort_seq(loni,1:3)=Nonzero(AcceptedSeg(loni),:);
   end
   % sort in the reverse order of the length,so the first is the longest
   % one,and we think it is the stablest 
   sort_seq=sortrows(sort_seq,-3);
   % desgin the gap value(the distance from the first Seg) 
   for loni=2:SegNum
       if sort_seq(loni,2)<sort_seq(1,1) 
          sort_seq(loni,4)=sort_seq(1,1)-sort_seq(loni,2);        
       elseif sort_seq(loni,1)>sort_seq(1,2)
          sort_seq(loni,4)=sort_seq(loni,1)-sort_seq(1,2);
       end
   end
   % then sort in order of the gap value
   sort_seq=sortrows(sort_seq,4);
   %  built the output array New
   for ii=2:SegNum
       if sort_seq(ii,1) > MaxTCal - MinTWidth 
          break;
       end

       TstartIX = find( RawDisper(:,1)==sort_seq(ii,1) );
       TendIX = find( RawDisper(:,1)==sort_seq(ii,2) );
       % RawDisper: [TPoint PhaseVDisper]
       % whether this segment and its previous one both belong to the same
       % dispersion curve
       % The dispersion curve should increases in a reasonable range
       if TendIX<PrevTstartIX % in front of the previous Seg
%         disp('Present Period IN FRONT');
          PhaseVGap =RawDisper(PrevTstartIX,2)- RawDisper(TendIX,2);
          TGap =RawDisper(PrevTstartIX,1)- RawDisper(TendIX,1);
          upper_PhaseVGap1 = mean(dcdt(max(TendIX-4,TstartIX):TendIX-1,2))*TGap*2;% for MinTWidth>4;Tstart<TendTX
          upper_PhaseVGap2 = PhaseVGapIntegral(Smooth_r_dcdt,RawDisper(TendIX,1),RawDisper(PrevTstartIX,1)); 
          upper_PhaseVGap = max(upper_PhaseVGap1,upper_PhaseVGap2);
          lower_PhaseVGap = PhaseVGapIntegral(Smooth_r_dcdt,RawDisper(TendIX,1),RawDisper(PrevTstartIX,1))/5+min(r_dcdt(:,3))/2*TGap;
       elseif PrevTendIX<TstartIX
          PhaseVGap = RawDisper(TstartIX,2) - RawDisper(PrevTendIX,2);
          TGap = RawDisper(TstartIX,1) - RawDisper(PrevTendIX,1);
          upper_PhaseVGap1 = mean(dcdt(max(PrevTendIX-4,PrevTstartIX):PrevTendIX-1,2))*TGap*2;
          upper_PhaseVGap2 = PhaseVGapIntegral(Smooth_r_dcdt,RawDisper(PrevTendIX,1),RawDisper(TstartIX,1)); 
          upper_PhaseVGap = max(upper_PhaseVGap1,upper_PhaseVGap2);
          lower_PhaseVGap = PhaseVGapIntegral(Smooth_r_dcdt,RawDisper(PrevTendIX,1),RawDisper(TstartIX,1))/5+min(r_dcdt(:,3))/2*TGap;
       end
    
       if PhaseVGap < upper_PhaseVGap && PhaseVGap > lower_PhaseVGap
          NewDisper=[NewDisper;RawDisper(TstartIX:TendIX,:)];
          
          if PrevTstartIX>TstartIX
             PrevTstartIX = TstartIX; 
          end
          if PrevTendIX<TendIX
             PrevTendIX = TendIX; 
          end
       end
   end
end
    
% Rewrite NewDisper
% NewDisper: T  V  quality
FullNewDisper = zeros(size(RawDisper,1),3);
FullNewDisper(:,1) = RawDisper(:,1);
for ii=1:size(FullNewDisper,1)
    T = FullNewDisper(ii,1);
    TIX = find( NewDisper(:,1) == T );
    if length(TIX) == 1
        FullNewDisper(ii,2) = NewDisper(TIX,2);
        FullNewDisper(ii,3) = 1;
    end
end
clear NewDisper
NewDisper = FullNewDisper;

% --- for Rayleigh waves
function NewDisper = AutomaticDisperRayl(RawDisper,SampleF,MinTWidth,MaxTCal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for rayleigh wave 
% This function extracts the reasonable part of a whole disper curve
% between two stations.
% RawDisper
% 1     2     3
% T   PhaseV  SNR
% Range of dcdt: if needed, this part can be put in input parameters.
% 1       2          3     
% T   upper-limit  lower-limit
r_dcdt=[  0  0.30  -0.20
         10  0.20  -0.15
         20  0.15  -0.10
         30  0.10  -0.01
         40  0.06  -0.005
         50  0.04  -0.003
        200  0.02  -0.003];   
    
r_dcdt=r_dcdt/SampleF;


% Create smooth range of dcdt
Smooth_r_dcdt = zeros(size(RawDisper,1)-1,3);
Smooth_r_dcdt(:,1) = RawDisper(1:end-1,1);
Smooth_r_dcdt(:,2) = interp1(r_dcdt(:,1),r_dcdt(:,2),Smooth_r_dcdt(:,1),'pchip');
Smooth_r_dcdt(:,3) = interp1(r_dcdt(:,1),r_dcdt(:,3),Smooth_r_dcdt(:,1),'pchip');

% dcdt: time derivative of phase velocity
% 1     2     3
% T    dcdt  g/b 
%
% if lower-limit <= dcdt < upper-limit : good = 1
% if 2 * lower-limit < dcdt < 1.5 * upper-limit  : half-good =1/2
% other cases : bad = 0

dcdt=zeros(size(RawDisper,1)-1,3);
dcdt(:,1)=RawDisper(1:end-1,1);
dcdt(:,2)=RawDisper(2:end,2)-RawDisper(1:end-1,2);
%%% quality control 1.
for ii=1:size(dcdt,1)
    upper_limit=Smooth_r_dcdt(ii,2);
    lower_limit=Smooth_r_dcdt(ii,3);
    
    if dcdt(ii,2)>lower_limit && dcdt(ii,2)<=upper_limit
        if dcdt(ii,1) >= 30  % for period less than 30 s, phase velocity of Rayleigh wave usually increases as T increases
            if dcdt(ii,2)<=0
                dcdt(ii,3)=0.5;
            else
                dcdt(ii,3)=1;
            end
        else  % for period < 30s
            dcdt(ii,3)=1;
        end     
    elseif dcdt(ii,2)>= 1.5*lower_limit && dcdt(ii,2)<=upper_limit*1.5
        dcdt(ii,3)=0.5;
    end
end
%%% quality control 2.
% Recheck data quality
% if d(dc/dt)/dt changes fast => bad 
for ii = 2:size(dcdt,1)
    ddcdtt = abs(dcdt(ii,2)-dcdt(ii-1,2));
    if ddcdtt > Smooth_r_dcdt(ii,2)-Smooth_r_dcdt(ii,3)
        dcdt(ii,3)=0;
    end
end
%%% quality control 3.    
% if there is less than 2 half-good pts in the middle of good pts, just
% accepted them
if SampleF == 1
    for ii=2:size(dcdt,1)-2
        if dcdt(ii,3)==0.5 && dcdt(ii-1,3)==1 && dcdt(ii+1,3)==1
            dcdt(ii,3)=1;
        elseif dcdt(ii,3)==0.5 && dcdt(ii-1,3)==1 && dcdt(ii+1,3)>0 && dcdt(ii+2,3)==1
            dcdt(ii,3)=1;
        end
    end
elseif SampleF == 2
    for ii = 2:size(dcdt,1)-SampleF
        if dcdt(ii,3)==0 && dcdt(ii-1,3)>0 && dcdt(ii+1,3)>0
            dcdt(ii,3)=0.5;
        end
    end
    for ii = 2:size(dcdt,1)-SampleF
        if dcdt(ii,3) == 0.5 && dcdt(ii-1,3)==1
            QaulityAfter=dcdt(ii+1:end,3);
            NearestOne=min(QaulityAfter==1);
            if sum(QaulityAfter(1:NearestOne))/NearestOne > 0.5
                dcdt(ii,3)=1;
            end
        end
    end
end
%%%% creat Nonzero 
% Calculate the length of each good-data segment: dcdt(ii,3)==1
% Nonzero:
%   1      2       3 
% Tstart  Tend   length
Nonzero=zeros(size(dcdt,1),3);
NonzeroIndex=0;
for ii=1:size(dcdt,1)
    if dcdt(ii,3)==1
        if ii==1
            NonzeroIndex=NonzeroIndex+1;
            Nonzero(NonzeroIndex,1)=dcdt(ii,1);
        elseif dcdt(ii-1,3)<1
            NonzeroIndex=NonzeroIndex+1;
            Nonzero(NonzeroIndex,1)=dcdt(ii,1);
        end
        Nonzero(NonzeroIndex,3)=Nonzero(NonzeroIndex,3)+1;
    end
end
Nonzero(:,2)=Nonzero(:,1)+Nonzero(:,3)/SampleF;

MinSegLength=MinTWidth*SampleF;  % choose segments containing at least MinTWidth data 

AcceptedSeg=find(Nonzero(:,3)>=MinSegLength);
SegNum=length(AcceptedSeg); % the nomber of the dispersion curve
% fprintf('the SegNum is %f \n',SegNum);
%%%% creat NewDisper to sotre the dispersion curve 
%%% NewDisper: T  C (phase velocity) and no SNR in 3rd column SNR, the same
%%% as the RawDisp array
%%% analysis every Seg belong to the same dispersion curve
if SegNum<1
    NewDisper=zeros(size(RawDisper,1),2);
    NewDisper(:,1) = RawDisper(:,1);
    return
end
% in fact, here is if SegNum=1;,the only one and the longest,
%FOLLOWING is changed
 if SegNum==1
  TstartIX=find( RawDisper(:,1)==Nonzero(AcceptedSeg(1),1) );
  TendIX=find( RawDisper(:,1)==Nonzero(AcceptedSeg(1),2) );
  NewDisper=RawDisper(TstartIX:TendIX,:);
 end
 if SegNum>=2
   % begin the judge from the longest Seg
   [longth,longestSeg] = max(Nonzero(:,3));
   PrevTstartIX=find(RawDisper(:,1) == Nonzero(longestSeg,1));
   PrevTendIX=find(RawDisper(:,1) == Nonzero(longestSeg,2));
   NewDisper=RawDisper(PrevTstartIX:PrevTendIX,:);

   % built an array of search sequnence: sort_seq
   % starT  endT  length  gap
   % in fact we built sort_seq to replace Nonzero
   sort_seq=zeros(SegNum,4);
   % assign for sort_seq 
   for loni=1:SegNum
       sort_seq(loni,1:3)=Nonzero(AcceptedSeg(loni),:);
   end
   % sort in the reverse order of the length,so the first is the longest
   % one,and we think it is the stablest 
   sort_seq=sortrows(sort_seq,-3);
   % desgin the gap value(the distance from the first Seg) for the 4th column
   for loni=2:SegNum
       if sort_seq(loni,2)<sort_seq(1,1) %the prior period
          sort_seq(loni,4)=sort_seq(1,1)-sort_seq(loni,2); % gap             
       elseif sort_seq(loni,1)>sort_seq(1,2)
          sort_seq(loni,4)=sort_seq(loni,1)-sort_seq(1,2);
       end
   end
   % then sort in order of the gap value(4th column) 
   sort_seq=sortrows(sort_seq,4);
   %  built the output array New
   for ii=2:SegNum
     if sort_seq(ii,1) > MaxTCal - MinTWidth %以参考速度在半波长法则下得到的MaxTCal最大周期，MinTWidth：最短周期段长度，一般情况不会出现
       continue;
     end

    TstartIX = find( RawDisper(:,1)==sort_seq(ii,1) );% 
    TendIX = find( RawDisper(:,1)==sort_seq(ii,2) );
    % whether this segment and its previous one both belong to the same
    % dispersion curve
    % The dispersion curve should increases in a reasonable range
    if TendIX<PrevTstartIX % in front of the previous Seg
%     disp('Present Period IN FRONT!!');
    PhaseVGap =RawDisper(PrevTstartIX,2)- RawDisper(TendIX,2);
    TGap =RawDisper(PrevTstartIX,1)- RawDisper(TendIX,1);
    upper_PhaseVGap1 = mean(dcdt(max(TendIX-4,TstartIX):TendIX-1,2))*TGap*2;%这里因为筛选长度是大于10的，所以，Tstart 肯定小于TendIX-4%
    upper_PhaseVGap2 = PhaseVGapIntegral(Smooth_r_dcdt,RawDisper(TendIX,1),RawDisper(PrevTstartIX,1)); 
    upper_PhaseVGap = max(upper_PhaseVGap1,upper_PhaseVGap2);
    % lower_PhaseVGap = PhaseVGapIntegral(Smooth_r_dcdt,RawDisper(PrevTendIX,1),RawDisper(TstartIX,1))/5;
    lower_PhaseVGap = PhaseVGapIntegral(Smooth_r_dcdt,RawDisper(TendIX,1),RawDisper(PrevTstartIX,1))/5;%
    elseif PrevTendIX<TstartIX
%     disp('Present Period BEHIND!!');
    PhaseVGap = RawDisper(TstartIX,2) - RawDisper(PrevTendIX,2);%本次符合条件的起始周期与上一连续周期段的截止周期的相速度差
    TGap = RawDisper(TstartIX,1) - RawDisper(PrevTendIX,1);%本次符合条件的起始周期与上一连续周期段的截止周期的周期跨度
    upper_PhaseVGap1 = mean(dcdt(max(PrevTendIX-4,PrevTstartIX):PrevTendIX-1,2))*TGap*2;
    upper_PhaseVGap2 = PhaseVGapIntegral(Smooth_r_dcdt,RawDisper(PrevTendIX,1),RawDisper(TstartIX,1)); %计算上下两连续周期段，之间周期的正斜率之和
    upper_PhaseVGap = max(upper_PhaseVGap1,upper_PhaseVGap2);
    lower_PhaseVGap = PhaseVGapIntegral(Smooth_r_dcdt,RawDisper(PrevTendIX,1),RawDisper(TstartIX,1))/5;  
    end
    
    if PhaseVGap < upper_PhaseVGap && PhaseVGap > lower_PhaseVGap%%%%%%
      NewDisper=[NewDisper;RawDisper(TstartIX:TendIX,:)];%
      % new assginment
      if PrevTstartIX>TstartIX
         PrevTstartIX = TstartIX; %
      end
      if PrevTendIX<TendIX
         PrevTendIX = TendIX;  % 
      end
    end
    
   end

 end  
% Rewrite NewDisper
FullNewDisper = zeros(size(RawDisper,1),3);
FullNewDisper(:,1) = RawDisper(:,1);
for ii=1:size(FullNewDisper,1)
    T = FullNewDisper(ii,1);
    TIX = find( NewDisper(:,1) == T );
    if length(TIX) == 1
        FullNewDisper(ii,2) = NewDisper(TIX,2);
        FullNewDisper(ii,3) = 1;
    end
end
clear NewDisper
NewDisper = FullNewDisper;%the return array 

% --- for PhaseVGapIntegral 
%  This function calculte the integral of a known function
function r_PhaseVGap=PhaseVGapIntegral(Smooth_r_dcdt,Tstart,Tend)

TstartIX = find(Smooth_r_dcdt(:,1)==Tstart);
TendIX = find(Smooth_r_dcdt(:,1)==Tend);
r_PhaseVGap = sum(Smooth_r_dcdt(TstartIX:TendIX-1,2));

% --- for save the automatic picked phase velocity disper
function SaveAutoPick_PhaseDisper(hObject, handles, VDisp,GoodIndex) 
global PhaseDispDirectory 
global cross src sta PrevStaLonLat

if PhaseDispDirectory == 0
    PhaseDispDirectory = pwd;
end
h1 = handles.axes1;
h2 = handles.axes2;
 
% make the file name
YYChar = num2str(src.YY);
MonChar = num2str(src.Month);
if length(MonChar) == 1
   MonChar = strcat('0',MonChar);
end
DDChar = num2str(src.Day);
if length(DDChar) == 1
   DDChar = strcat('0',DDChar);
end   
HHChar = num2str(src.HH);
if length(HHChar) == 1
   HHChar = strcat('0',HHChar);
end
MMChar = num2str(src.MM);
if length(MMChar) == 1
   MMChar = strcat('0',MMChar);
end   
SrcTime = strcat(YYChar(3:4),MonChar,DDChar,'_',HHChar,MMChar);
% stationpar name
for i = 1:max(size(sta(1).Name, 2), size(sta(2).Name, 2))
    if double(sta(1).Name(i)) > double(sta(2).Name(i))
       StaNamePair = strcat(sta(2).Name,'-', sta(1).Name);
       index1 = 2;
       index2 = 1;
    break
    elseif double(sta(1).Name(i)) < double(sta(2).Name(i))
           StaNamePair = strcat(sta(1).Name,'-', sta(2).Name);
           index1 = 1;
           index2 = 2;
           break
     else
     continue
     end
end
% tie the processing stations
if get(handles.Silent_cycle_mode,'Value')==0
   if sum(sum(PrevStaLonLat)) ~= 0
      set(gcf,'CurrentAxes',h1);
      hold(h1,'on');
      plot(PrevStaLonLat(1:2,1), PrevStaLonLat(1:2,2), 'b-');
      plot(h1, PrevStaLonLat(1:2,1), PrevStaLonLat(1:2,2), 'r^', 'MarkerSize',6, 'MarkerFaceColor','r');

   end
end
PrevStaLonLat = [sta(index1).Lon, sta(index1).Lat;sta(index2).Lon, sta(index2).Lat ];
% make sure station distance must exceed half of the wavelength
DataExistIndex = 0;
TPoint = cross.StartT:cross.DeltaT:cross.EndT;
StaDistance = abs(sta(1).GCDkm - sta(2).GCDkm);

for i = 1:cross.NumCtrT
    minlamdaRatio=2;    
    wavelength = VDisp(i)*TPoint(i);
    if (wavelength<=minlamdaRatio*StaDistance) && GoodIndex(i) == 1
       DataExistIndex = DataExistIndex+ 1;    
    else
       GoodIndex(i) = 0;     
    end
end
II  = find(GoodIndex == 1);

% If the existing date points is more than 5, then save them.
if DataExistIndex >= 10
   dispfilename = strcat('C' ,SrcTime, StaNamePair, '.dat');
   TVfile = fullfile(PhaseDispDirectory, dispfilename);
   
   ftv = fopen(TVfile,'w');
   fprintf(ftv,'%f     ', sta(index1).Lon);
   fprintf(ftv,'%f\n', sta(index1).Lat);
   fprintf(ftv,'%f     ', sta(index2).Lon);
   fprintf(ftv,'%f\n', sta(index2).Lat);
   for i = 1:cross.NumCtrT
       fprintf(ftv,'%4.1f    %4.3f\n',[TPoint(i)  VDisp(i)*GoodIndex(i)]);
   end
   fclose(ftv);
% plot the saveed dispersion curve point finally
   if get(handles.Silent_cycle_mode,'Value')~=1
      set(gcf,'CurrentAxes',h2);
      hold(h2,'on');
      plot(h2,TPoint(II), VDisp(II),'ro','MarkerSize',4,'MarkerFaceColor','r'); 
      
      set(gcf,'CurrentAxes',h1);
      hold(h1,'on');   
      plot([sta(index1).Lon, sta(index2).Lon], [sta(index1).Lat, sta(index2).Lat], 'r-');
      plot(h1, [sta(index1).Lon, sta(index2).Lon], [sta(index1).Lat, sta(index2).Lat], 'k^', 'MarkerSize',6, 'MarkerFaceColor','g');
   end
   words=['Output Phase velocity dispersion file: ',dispfilename];
   if get(handles.Semi_single_mode,'Value')==1 || get(handles.Auto_single_mode,'Value')==1
      disp(words);
   elseif get(handles.Visible_cycle_mode,'Value')==1 || get(handles.Silent_cycle_mode,'Value')==1
       disp(words);
       loglog(words);
   end
end
   

%%%------------------------------------------------------------------------
%%%                 Process Control Module  
%%%------------------------------------------------------------------------ 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%













% Hints: get(hObject,'String') returns contents of deltaPhaseV as text
%        str2double(get(hObject,'String')) returns contents of deltaPhaseV as a double
