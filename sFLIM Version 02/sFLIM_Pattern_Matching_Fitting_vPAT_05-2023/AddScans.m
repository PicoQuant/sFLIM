function varargout = AddScans(varargin)
% AddScans MATLAB code for AddScans.fig
%      AddScans, by itself, creates a new AddScans or raises the existing
%      singleton*.
%
%      H = AddScans returns the handle to a new AddScans or the handle to
%      the existing singleton*.
%
%      AddScans('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AddScans.M with the given input arguments.
%
%      AddScans('Property','Value',...) creates a new AddScans or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AddScans_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AddScans_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AddScans

% Last Modified by GUIDE v2.5 18-Jun-2013 11:02:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AddScans_OpeningFcn, ...
                   'gui_OutputFcn',  @AddScans_OutputFcn, ...
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


% --- Executes just before AddScans is made visible.
function AddScans_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AddScans (see VARARGIN)

% Choose default command line output for AddScans
handles.output = hObject;

handles.par_file = {};
handles.par_IRFfile = '';
handles.filename = {};

set(handles.text2,'String','Ready');

wt = num2cell(25);
wt(1) = {548};
wt(2) = {60};
set(handles.uitable1,'ColumnWidth',wt);

col_name = {'Files','weight'};
set(handles.uitable1,'Data',handles.filename');
set(handles.uitable1,'ColumnName',col_name);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AddScans wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AddScans_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Load_button.
function Load_button_Callback(hObject, eventdata, handles)
% hObject    handle to Load_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

n = numel(handles.par_file);

[FileName,PathName] = uigetfile({'*.ht3;*.pt3','t3r-files';'*.*','All Files' },'Select data files ...','MultiSelect', 'on');

if ~iscell(FileName)
    if FileName ~= 0    
        handles.par_file{n+1} = [PathName FileName];
        handles.filename{n+1} = FileName;
    end
else
    for k = 1:numel(FileName)
        handles.par_file{n+k} = [PathName FileName{k}];
        handles.filename{n+k} = FileName{k};
    end
end
for n = 1:numel(handles.filename)
    data{n,1} = handles.filename{n};
    data{n,2} = 1;

    handles.f_weight(n) = 1;

end

set(handles.uitable1,'Data',data);
guidata(hObject, handles);


% --- Executes on button press in SetParams_button.
function SetParams_button_Callback(hObject, eventdata, handles)
% hObject    handle to SetParams_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.figure1,'Visible','off');
Set_Params;
set(handles.figure1,'Visible','on');


% --- Executes on button press in Start_button.
function Start_button_Callback(hObject, eventdata, handles)
% hObject    handle to Start_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text2,'String','Processing data...');
k = numel(handles.par_file);    

for n = 1:k
    
    W  = handles.par_file{n};
    s = sprintf('Processing file %s  (%d of %d)...', handles.filename{n}, n, k);        
    
    set(handles.text2,'String', s );

    pause(0.2)
    
    N1 = [W(1:end-4) '_FLIM.mat'];
    tmp = dir(N1);

    if (numel(tmp)==0)        
        Process_file(W, handles.par_IRFfile);
    end
    
end

for n = 1:k
    
    W = handles.par_file{n};
    
    load([W(1:end-4) '_DATA.mat'], 'head', 'tag', 'tcspc', 'IRF', 'AP', 'DC');
    
    nx(n)      = size(tag,1);
    ny(n)      = size(tag,2);    
    nch(n)     = size(tcspc,1);
    nbin(n)    = size(tcspc,2);
    num_PIE(n) = size(tcspc,3);
    
end

ind = (nx~=nx(1))|(ny~=ny(1));
if sum(ind)==0
    ind = (nch~=nch(1))|(nbin~=nbin(1))|(num_PIE~=num_PIE(1));
    if sum(ind)==0
        nx      = nx(1);
        ny      = ny(1);
        nch     = nch(1);
        nbin    = nbin(1);
        num_PIE = num_PIE(1);
        
        numpix  = nx*ny;
        numbin  = nch*nbin*num_PIE;
        numimg  = k;
        
        s = sprintf('Computation in progress...');
        
    else
        s = sprintf('TCSPC formats are not matching! Aborted.');
    end
else
    s = sprintf('Image sizes are not matching! Aborted.'); 
end

set(handles.text2,'String', s );
pause(0.2)

p = handles.f_weight./sum(handles.f_weight);

sname = [];

raw_data = zeros(nx,ny,numbin,numimg);

for n = 1:numimg
    
    W = handles.par_file{n};
    
    load([W(1:end-4) '_FLIM.mat'], 'stim')

    stim = reshape(stim, [nx ny numbin]);             % rearrange the data for convenience
    raw_data(:,:,:,n) = stim;
    clear stim;
    
    tmp = sprintf('%s (%5.2f%%) ',handles.filename{n}(1:end-4),100.*p(n));
    
    sname = [sname tmp];
end

raw_data(raw_data<0) = 0;

% p = round(squeeze(sum(raw_data,3)).*repmat(p,[nx ny 1]));

stim = zeros(nx,ny,numbin);

h = waitbar(0,'Progress:');

for n1 = 1:nx
    
    for n2 = 1:ny
        
        tcspc = zeros(1,numbin);
        
        for k = 1:numimg
            
            tmp = round(raw_data(n1,n2,:,k));
            
            n_ph = sum(tmp);
            
            p_n  = zeros(1,n_ph);
            bin  = 1:numbin;
            bin(tmp==0) = [];
            
            a = 1;
            for l = 1:numel(bin)
                e = a+tmp(bin(l));
                p_n(a:e-1) = bin(l);
                a = e;
            end
            
            for j = 1:10
                p_n = p_n(randperm(n_ph));
            end
            
            for j = 1:round(n_ph*p(k))
                bin = p_n(ceil(n_ph.*rand(1)));
                tcspc(bin) = tcspc(bin) +1;
            end
        end
        stim(n1,n2,:) = tcspc;
    end
    waitbar(n1/nx,h);
    drawnow
end
close(h);

stim  = reshape(stim, [nx ny nch nbin num_PIE]);
tcspc = shiftdim(sum(sum(stim,1),2),2);
tcspc = tcspc./repmat(head.tauw,[nch 1 num_PIE]);
tag   = squeeze(sum(stim,4));

pathname = '';
tmp = strfind(handles.par_file{numimg},'\');
if ~isempty(tmp)
    pathname = handles.par_file{numimg}(1:tmp(end));
end

head.CreatorName = sprintf('%-18s','AddScans');
head.CreatorVersion = sprintf('%-12s','1.0.0.0');
tmp = clock;
head.FileTime = sprintf('%02d/%02d/%02d %02d:%02d:%02d ',tmp(3),tmp(2),tmp(1)-2000,tmp(4),tmp(5),round(tmp(6)));

save([pathname sname(1:end-1) '_DATA.mat'], 'head', 'tag', 'tcspc', 'IRF', 'AP', 'DC');
save([pathname sname(1:end-1) '_FLIM.mat'], 'stim');

fid = fopen([pathname sname(1:end-1) '.pt3'],'w');
fwrite(fid, sprintf('%-16s',head.Ident), 'char');                     %   Ident (1-..)                  4
fwrite(fid, head.FormatVersion, 'char');                     %   FormatVersion (1-..)         20
fwrite(fid, head.CreatorName, 'char');           %   CreatorName (1-18)           26
fwrite(fid, head.CreatorVersion, 'char');             %   CreatorVersion (1-7)         44
fwrite(fid, head.FileTime, 'char');              %   FileTime    (1-17)           56

a = [13 10];
fwrite(fid, a, 'char');                     %   'CR/LF'                      74

fwrite(fid, head.Comment, 'char');                     %   Comment (256)                76

fwrite(fid, head.NCurves, 'int32');               %   NCurves (1)                 332
fwrite(fid, head.BitsPerRecord, 'int32');               %   BitsPerRecord (1)           336
fwrite(fid, head.RoutingChannels, 'int32');               %   RoutingChannels (1)         340
fwrite(fid, head.NumberOfBoards, 'int32');               %   NumberOfBoards (1)          344
fwrite(fid, head.ActiveCurve, 'int32');               %   ActiveCurve (1)             348
fwrite(fid, head.MeasMode, 'int32');               %   MeasMode (1)                352
fwrite(fid, head.SubMode, 'int32');               %   SubMode  (1)                356
fwrite(fid, head.RangeNo, 'int32');               %   RangeNo (1)                 360
fwrite(fid, head.Offset, 'int32');               %   Offset (1)                  364
fwrite(fid, 0, 'int32');               %   TAcq (1)                    368
fwrite(fid, head.StopAt, 'int32');               %   StopAt (1)                  372
fwrite(fid, head.StopOnOvfl, 'int32');               %   StopOnOvfl (1)              376
fwrite(fid, head.Restart, 'int32');               %   Restart (1)                 380
fwrite(fid, head.LinLog, 'int32');               %   DisplayLinLog (1)           384
fwrite(fid, head.MinAx, 'int32');               %   DisplayTimeAxMin (1)        388
fwrite(fid, head.MaxAx, 'int32');               %   DisplayTimeAxMax (1)        392
fwrite(fid, head.MinAxCnt, 'int32');               %   DisplayCountAxMin (1)       396
fwrite(fid, head.MaxAxCnt, 'int32');               %   DisplayCountAxMax (1)       400

for x = 1:8
    fwrite(fid, head.DispCurves(2*x-1), 'int32');                %   DisplayCurve[x].MapTo (8x1) 404 412 420 428 436 444 452 460
    fwrite(fid, head.DispCurves(2*x  ), 'int32');                %   DisplayCurve[x].Show  (8x1) 408 416 424 432 440 448 456 464
end;

for x = 1:3
    fwrite(fid, head.Params(3*x-2), 'float');                %   Param[x].Start        (3x1) 468 480 492
    fwrite(fid, head.Params(3*x-1), 'float');                %   Param[x].Step         (3x1) 472 484 496
    fwrite(fid, head.Params(3*x  ), 'float');                %   Param[x].End          (3x1) 476 488 500
end;

fwrite(fid, head.RepeatMode, 'int32');                    %   RepeatMode (1)              504
fwrite(fid, head.RepeatsPerCurve, 'int32');                    %   RepeatsPerCurve (1)         508
fwrite(fid, head.RepeatTime, 'int32');                    %   RepeatTime (1)              512
fwrite(fid, head.RepeatWaitTime, 'int32');                    %   RepeatWaitTime (1)          516

fwrite(fid, head.ScriptName, 'char');       %   ScriptName (1-13)           520

for x = 1:head.NumberOfBoards
    fwrite(fid, head.BoardIdent(:,x), 'char');                     %   HardwareIdent (1-..)        540
    fwrite(fid, head.BoardVersion(:,x), 'char');                     %   HardwareIdent (1-..)        540
    fwrite(fid, head.BoardSerial(x), 'int32');             %   BordSerial (1)              564
    fwrite(fid, head.SyncDiv(x), 'int32');              %   SyncDivider                 568
    fwrite(fid, head.CFDZeroCross0(x), 'int32');            %   CFD_Zero 0 (1)              572
    fwrite(fid, head.CFDLevel0(x), 'int32');            %   CFDLevel 0 (1)              576
    fwrite(fid, head.CFDZeroCross1(x), 'int32');            %   CFD_Zero 1 (1)              580
    fwrite(fid, head.CFDLevel1(x), 'int32');            %   CFDLevel 1 (1)              584
    fwrite(fid, head.Resolution(x), 'float');        %   Resolution (1)              588
    fwrite(fid, head.RouterModel(x), 'int32');        
    fwrite(fid, head.RouterEnabled(x), 'int32');        
    
    fwrite(fid, head.RtCh1_InputType(x), 'int32');        
    fwrite(fid, head.RtCh1_InputLevel(x), 'int32');        
    fwrite(fid, head.RtCh1_InputEdge(x), 'int32');        
    fwrite(fid, head.RtCh1_CFDPresent(x), 'int32');        
    fwrite(fid, head.RtCh1_CFDLevel(x), 'int32');        
    fwrite(fid, head.RtCh1_CFDZeroCross(x), 'int32');        
    fwrite(fid, head.RtCh2_InputType(x), 'int32');        
    fwrite(fid, head.RtCh2_InputLevel(x), 'int32');        
    fwrite(fid, head.RtCh2_InputEdge(x), 'int32');        
    fwrite(fid, head.RtCh2_CFDPresent(x), 'int32');        
    fwrite(fid, head.RtCh2_CFDLevel(x), 'int32');        
    fwrite(fid, head.RtCh2_CFDZeroCross(x), 'int32');        
    fwrite(fid, head.RtCh3_InputType(x), 'int32');        
    fwrite(fid, head.RtCh3_InputLevel(x), 'int32');        
    fwrite(fid, head.RtCh3_InputEdge(x), 'int32');        
    fwrite(fid, head.RtCh3_CFDPresent(x), 'int32');        
    fwrite(fid, head.RtCh3_CFDLevel(x), 'int32');        
    fwrite(fid, head.RtCh3_CFDZeroCross(x), 'int32');        
    fwrite(fid, head.RtCh4_InputType(x), 'int32');        
    fwrite(fid, head.RtCh4_InputLevel(x), 'int32');        
    fwrite(fid, head.RtCh4_InputEdge(x), 'int32');        
    fwrite(fid, head.RtCh4_CFDPresent(x), 'int32');        
    fwrite(fid, head.RtCh4_CFDLevel(x), 'int32');        
    fwrite(fid, head.RtCh4_CFDZeroCross(x), 'int32');        
end

fwrite(fid,  head.ExternalDev, 'int32');                   %   External Devices (1)        592
fwrite(fid,  head.Reserved1, 'int32');                   %   Dummy (1)                   596
fwrite(fid,  head.Reserved2, 'int32');                   %   Dummy (1)                   600

fwrite(fid,  0, 'int32');                   %   Syncrate (1)                604
fwrite(fid,  0, 'int32');                   %   Countrate (1)               608
fwrite(fid,  0, 'int32');                   %   TTTRStopAfter (1)           612
fwrite(fid,  0, 'int32');                   %   TTTRStopReason (1)          616
fwrite(fid,  0, 'uint32');                  %   NCounts (1)                 620

fwrite(fid, head.SpecHeaderLength, 'uint32');                  %   SpecHeaderLength (1)        624    
fwrite(fid,  head.ImgHdr.Dimensions, 'int32');                   %   Dimensions (3=area)    (1)  628
fwrite(fid,  head.ImgHdr.Ident, 'int32');                   %   Ident  (1=E710)        (1)  632
fwrite(fid,  1+log2(head.ImgHdr.Frame), 'int32');     
fwrite(fid,  1+log2(head.ImgHdr.LineStart), 'int32');
fwrite(fid,  1+log2(head.ImgHdr.LineStop), 'int32');
fwrite(fid,  head.ImgHdr.Pattern, 'int32'); 
fwrite(fid,  head.ImgHdr.PixX, 'int32');           
fwrite(fid,  head.ImgHdr.PixY, 'int32');           

fclose(fid);

set(handles.text2,'String','Done');


% --- Executes on button press in Quit_button.
function Quit_button_Callback(hObject, eventdata, handles)
% hObject    handle to Quit_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes when entered data in editable cell(s) in uitable1.
function uitable1_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

data = get(handles.uitable1,'Data');
if eventdata.Indices(2) == 2
    if isnan(eventdata.NewData)
        data{eventdata.Indices(1),2} = 0;
        set(handles.uitable1,'Data',data);
    end
end

for k = 1:size(data,1)
    handles.f_weight(k) = data{k,2};
end

guidata(hObject, handles);


