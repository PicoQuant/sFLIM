function varargout = sFLIM(varargin)
% SFLIM M-file for sFLIM.fig
%      SFLIM, by itself, creates a new SFLIM or raises the existing
%      singleton*.
%
%      H = SFLIM returns the handle to a new SFLIM or the handle to
%      the existing singleton*.
%
%      SFLIM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SFLIM.M with the given input arguments.
%
%      SFLIM('Property','Value',...) creates a new SFLIM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sFLIM_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sFLIM_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sFLIM

% Last Modified by GUIDE v2.5 15-Sep-2011 09:35:37
%#function Set_Params 

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sFLIM_OpeningFcn, ...
                   'gui_OutputFcn',  @sFLIM_OutputFcn, ...
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


% --- Executes just before sFLIM is made visible.
function sFLIM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sFLIM (see VARARGIN)

% Choose default command line output for sFLIM
handles.output = hObject;

% Update handles structure
set(handles.text13,'String','Ready');

handles.par_file = '';
handles.par_IRFfile = '';

guidata(hObject, handles);


% UIWAIT makes sFLIM wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = sFLIM_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function input_filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in FastFLIM.
function FastFLIM_Callback(hObject, eventdata, handles)
% hObject    handle to FastFLIM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%#function FastFLIM 

set(handles.text13,'String','Busy');
pause(0.05)
evalin('base', 'clear all');

W  = handles.par_file;
N1 = [W(1:end-4) '_Moments.mat'];

tmp = dir(N1);

if (numel(tmp)==0)
    s = sprintf('MomentAnalysis(\''%s\'');', handles.par_file);
    evalin('base', s);
end

s = sprintf('FastFLIM(\''%s\'');', handles.par_file);
evalin('base', s);
set(handles.text13,'String','Ready');


% --- Executes on button press in PatternAnalysis
function pushbutton30_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%#function MomentAnalysis Results_4 

set(handles.text13,'String','Busy');
pause(0.05)
W  = handles.par_file;
N1 = [W(1:end-4) '_Moments.mat'];

tmp = dir(N1);

if (numel(tmp)==0)
    s = sprintf('MomentAnalysis(\''%s\'');', handles.par_file);
    evalin('base', s);
end

s = sprintf('Results_4(\''file\'', \''%s\'');', handles.par_file);
evalin('base', s);
set(handles.text13,'String','Ready');


% --- Executes on button press in Process_file.
function Process_file_Callback(hObject, eventdata, handles)
% hObject    handle to Process_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~iscell(handles.par_file)
    n = 1;
    W  = handles.par_file;
    N1 = [W(1:end-4) '_DATA.mat'];

    tmp = dir(N1);

    if (numel(tmp)==0)
        set(handles.text13,'String','Processing data...');
        pause(0.05)
        Process_file(W, handles.par_IRFfile);
    end
else
    for n = 1:numel(handles.par_file)
        W  = handles.par_file{n};
        N1 = [W(1:end-4) '_DATA.mat'];

        tmp = dir(N1);

        if (numel(tmp)==0)
            set(handles.text13,'String','Processing data...');
            pause(0.05)
            Process_file(W, handles.par_IRFfile);
        end
    end
end

if n>1
    tmp = questdlg([{'Please select recording scheme:'}, {'Same excitation, different emission channel.'}, {'Different excitaition, same emission channel.'}],'Combine images','Same excitation','Same emission','Cancel','Yes');
    if ~strcmp(tmp,'Cancel')

        FileName = get(handles.text10, 'string');
        FileName = [FileName(1:end-7) '_combined.ht3'];        
        set(handles.text10, 'string', FileName);

        W  = handles.par_file{1};
        load([W(1:end-4) '_DATA.mat']);

        tmp_DC     = DC;
        tmp_AP     = AP;
        tmp_IRF    = IRF;
        tmp_tcspc  = tcspc;
        tmp_tag    = tag;
        
        for n = 2:numel(handles.par_file)
            W  = handles.par_file{n};
            load([W(1:end-4) '_DATA.mat']);
            if strcmp(tmp,'Same emission')                
                tmp_DC    = cat(3, tmp_DC, DC);
                tmp_AP    = cat(3, tmp_AP, AP);
                tmp_IRF   = cat(3, tmp_IRF, IRF);
                tmp_tcspc = cat(3, tmp_tcspc, tcspc);
                tmp_tag   = cat(4, tmp_tag, tag);
            else
                tmp_DC    = cat(2, tmp_DC, DC);
                tmp_AP    = cat(2, tmp_AP, AP);
                tmp_IRF   = cat(1, tmp_IRF, IRF);
                tmp_tcspc = cat(1, tmp_tcspc, tcspc);
                tmp_tag   = cat(3, tmp_tag, tag);                
            end
        end
        
        DC    = tmp_DC;
        AP    = tmp_AP;
        IRF   = tmp_IRF;
        tcspc = tmp_tcspc;
        tag   = tmp_tag;
        
        save([handles.par_file{1}(1:end-4) '_combined_DATA.mat'], 'head', 'DC', 'AP', 'IRF', 'tcspc', 'tag');
        
        clear('head', 'DC', 'AP', 'IRF', 'tcspc', 'tag', 'tmp_DC', 'tmp_AP', 'tmp_IRF', 'tmp_tcspc', 'tmp_tag');
        
        load([handles.par_file{1}(1:end-4) '_FLIM.mat']);

        tmp_stim = stim;
        for n = 2:numel(handles.par_file)
            W  = handles.par_file{n};
            load([W(1:end-4) '_FLIM.mat']);
            if strcmp(tmp,'Same emission')                
                tmp_stim = cat(5, tmp_stim, stim);
            else
                tmp_stim = cat(3, tmp_stim, stim);
            end
        end
        
        stim = tmp_stim;

        save([handles.par_file{1}(1:end-4) '_combined_FLIM.mat'], 'stim');

        clear('stim','tmp_stim');
        
        handles.par_file = [handles.par_file{1}(1:end-4) '_combined.ht3']; 
    else
       handles.par_file = handles.par_file{1}; 
    end
end

set(handles.text13,'String','Ready');
guidata(hObject, handles);


% --- Executes on button press in FLIM_graph.
function FLIM_graph_Callback(hObject, eventdata, handles)
evalin('base', 'fastflim_graph');
% hObject    handle to FLIM_graph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function input_filename_Callback(hObject, eventdata, handles)
% hObject    handle to input_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
name = get(hObject, 'String');
% Hints: get(hObject,'String') returns contents of input_filename as text
%        str2double(get(hObject,'String')) returns contents of input_filename as a double


% --- Executes on button press in Select_ht3file.
function Select_ht3file_Callback(hObject, eventdata, handles)
% hObject    handle to Select_ht3file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.text10, 'string', '');
handles.par_file = '';

[FileName,PathName] = uigetfile({'*.ht3;*.pt3;*.ptu','t3r-files';'*.*','All Files'}, 'MultiSelect', 'on');

if ~iscell(FileName)
    set(handles.text10, 'string', FileName);
    handles.par_file = [PathName FileName];
elseif numel(FileName) > 1    
    if strcmp(questdlg([{'You have selceted more than one file!'}, {'Processing will combine the data in one image. Continue?'}],'Confirm selection','Yes','No','No'),'Yes')
        set(handles.text10, 'string', [FileName{1} '...']);
        for n = 1:numel(FileName)
            handles.par_file{n} = [PathName FileName{n}];
        end
    end   
else
   warndlg('No files selected!')
end

guidata(hObject, handles);


% --- Executes on button press in select_IRF_file.
function select_IRF_file_Callback(hObject, eventdata, handles)
% hObject    handle to select_IRF_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.text14, 'string', '');
handles.par_IRFfile = '';

[FileName,PathName] = uigetfile({'*.ht3;*.pt3;*.ptu','tr3-files';'*.mat','mat-files' });

if isempty(FileName)
    warndlg('No file selected!')
else
    if ~iscell(FileName)
        set(handles.text14, 'string', FileName);
        handles.par_IRFfile = [PathName FileName];
    elseif numel(FileName) > 1
        set(handles.text14, 'string', FileName{1});
        handles.par_IRFfile = [PathName FileName{1}];
    end
end
guidata(hObject, handles);


% --- Executes on button press in Clear_files.
function Clear_files_Callback(hObject, eventdata, handles)
% hObject    handle to Clear_files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isempty(handles.par_file)    
   if strcmp(questdlg('This will delete all associated data files. Continue?','Confirm deletion','Yes','No','No'),'Yes')
       delete([handles.par_file(1:end-4) '_*.*']);
   end
else
  helpdlg('Please select sFLIM data-file.','No data-file chosen');  
end


% --- Executes on button press in Close.
function Close_Callback(hObject, eventdata, handles)
% hObject    handle to Close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close all


% --- Executes on button press in Set_Params.
function Set_Params_Callback(hObject, eventdata, handles)
% hObject    handle to Set_Params (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%#function Set_Params 

set(handles.figure1,'Visible','off');

%evalin('base', 'Set_Params')
Set_Params;
set(handles.figure1,'Visible','on');


% --- Executes on button press in pushbutton27.
function pushbutton27_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
evalin('base', 'clear all')


% --- Executes on button press in pushbutton28.
function pushbutton28_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text13,'String','Bye');
close;
delete(handles.figure1);


% --- Executes during object creation, after setting all properties.
function text13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'String','Please select PTU-file');

% --- Executes during object deletion, before destroying properties.
function text10_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to text10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'String','');


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);
