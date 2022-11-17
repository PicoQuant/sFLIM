function varargout = TheProcessor(varargin)
% THEPROCESSOR MATLAB code for TheProcessor.fig
%      THEPROCESSOR, by itself, creates a new THEPROCESSOR or raises the existing
%      singleton*.
%
%      H = THEPROCESSOR returns the handle to a new THEPROCESSOR or the handle to
%      the existing singleton*.
%
%      THEPROCESSOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in THEPROCESSOR.M with the given input arguments.
%
%      THEPROCESSOR('Property','Value',...) creates a new THEPROCESSOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TheProcessor_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TheProcessor_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TheProcessor

% Last Modified by GUIDE v2.5 07-Sep-2017 14:22:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TheProcessor_OpeningFcn, ...
                   'gui_OutputFcn',  @TheProcessor_OutputFcn, ...
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


% --- Executes just before TheProcessor is made visible.
function TheProcessor_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TheProcessor (see VARARGIN)

% Choose default command line output for TheProcessor
handles.output = hObject;

handles.par_file = {};
handles.par_IRFfile = '';
handles.filename = {};

set(handles.text2,'String','Ready');

wt = num2cell(25);
wt(1) = {510};
wt(2) = {0};
set(handles.uitable1,'ColumnWidth',wt);

col_name = {'Files'};
set(handles.uitable1,'Data',handles.filename');
set(handles.uitable1,'ColumnName',col_name);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TheProcessor wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = TheProcessor_OutputFcn(hObject, eventdata, handles) 
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

[FileName,PathName] = uigetfile({'*.ht3;*.pt3;*.ptu','t3r-files';'*.*','All Files' },'Select data files ...','MultiSelect', 'on');

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
set(handles.uitable1,'Data',handles.filename');
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

parfor n = 1:k
    
    W  = handles.par_file{n};
    if isempty(handles.par_IRFfile)
        s = sprintf('Processing file %s  (%d of %d)...', handles.filename{n}, n, k);        
    else
        s = sprintf('Processing file %s (%d of %d) using IRF "%s"...', handles.filename{n}, n, k, handles.par_IRFfile);
    end
    
    set(handles.text2,'String', s );

    pause(0.2)
    
    N1 = [W(1:end-4) '_DATA.mat'];
    tmp = dir(N1);

    if (numel(tmp)==0)        
        Process_file(W, handles.par_IRFfile);
    end
    
end

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


% --- Executes on button press in SelectIRF_button.
function SelectIRF_button_Callback(hObject, eventdata, handles)
% hObject    handle to SelectIRF_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.par_IRFfile = '';
handles.IRFfile = '';

[FileName, PathName] = uigetfile({'*.ht3;*.pt3','t3r-files';'*.*','All Files' });

if FileName ~= 0
    handles.par_IRFfile = [PathName FileName];
    handles.IRFfile = FileName;
end
guidata(hObject, handles);
