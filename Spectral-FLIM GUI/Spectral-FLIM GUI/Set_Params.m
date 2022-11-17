function varargout = Set_Params(varargin)
% SET_PARAMS M-file for Set_Params.fig
%      SET_PARAMS, by itself, creates a new SET_PARAMS or raises the existing
%      singleton*.
%
%      H = SET_PARAMS returns the handle to a new SET_PARAMS or the handle to
%      the existing singleton*.
%
%      SET_PARAMS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SET_PARAMS.M with the given input arguments.
%
%      SET_PARAMS('Property','Value',...) creates a new SET_PARAMS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Set_Params_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Set_Params_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Set_Params

% Last Modified by GUIDE v2.5 06-Mar-2014 13:18:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Set_Params_OpeningFcn, ...
                   'gui_OutputFcn',  @Set_Params_OutputFcn, ...
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


% --- Executes just before Set_Params is made visible.
function Set_Params_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Set_Params (see VARARGIN)

% Choose default command line output for Set_Params
handles.output = hObject;

% Update handles structure

lam_start = 490;
lam_step  = 8;
tau_min   = 1;
tau_max   = 4;
cum2_min  = -5;
cum2_max  = 5;
clam_min  = 490;
clam_max  = 650;
wlam_min  = 0;
wlam_max  = 30;
num_PIE   = 1;
Threshold = 30;
pixwidth  = 0;
binning   = 1;
BG_corr   = 1;
BG_AP     = 0.0;
lamPIE    = [475 530 560 640];

load('Get_Params.mat');

handles.par_lam_start = lam_start;
handles.par_lam_step  = lam_step;
handles.par_tau_min   = tau_min;
handles.par_tau_max   = tau_max;
handles.par_cum2_min  = cum2_min;
handles.par_cum2_max  = cum2_max;
handles.par_clam_min  = clam_min;
handles.par_clam_max  = clam_max;
handles.par_wlam_min  = wlam_min;
handles.par_wlam_max  = wlam_max;
handles.par_num_PIE   = num_PIE;
handles.par_Threshold = Threshold;
handles.par_binning   = binning;
handles.par_pixwidth  = pixwidth;
handles.par_BG_corr   = BG_corr;
handles.par_BG_AP     = 100*BG_AP;
handles.par_lamPIE    = lamPIE;

guidata(hObject, handles);

set(handles.edit1 ,'String',sprintf('%d',round(handles.par_lam_start)));
set(handles.edit2 ,'String',sprintf('%d',round(handles.par_lam_step)));
set(handles.edit3 ,'String',sprintf('%d',round(handles.par_num_PIE)));
set(handles.edit5 ,'String',sprintf('%d',round(handles.par_wlam_max)));
set(handles.edit6 ,'String',sprintf('%d',round(handles.par_wlam_min)));
set(handles.edit7 ,'String',sprintf('%4.1f',round(10*handles.par_tau_min)/10));
set(handles.edit8 ,'String',sprintf('%4.1f',round(10*handles.par_tau_max)/10));
set(handles.edit9 ,'String',sprintf('%d',round(handles.par_cum2_min)));
set(handles.edit10,'String',sprintf('%d',round(handles.par_clam_min)));
set(handles.edit11,'String',sprintf('%d',round(handles.par_clam_max)));
set(handles.edit12,'String',sprintf('%d',round(handles.par_cum2_max)));
set(handles.edit13,'String',sprintf('%d',round(handles.par_Threshold)));
set(handles.edit14,'String',sprintf('%4.2f',(handles.par_pixwidth)));
set(handles.edit21,'String',sprintf('%d',round(handles.par_binning)));
set(handles.edit15,'String',sprintf('%4.2f',handles.par_BG_AP));
set(handles.edit16,'String',sprintf('%d',handles.par_lamPIE(1)));
set(handles.edit18,'String',sprintf('%d',handles.par_lamPIE(2)));
set(handles.edit19,'String',sprintf('%d',handles.par_lamPIE(4)));
set(handles.edit20,'String',sprintf('%d',handles.par_lamPIE(3)));
set(handles.checkbox1,'Value',(handles.par_BG_corr>0));
if handles.par_BG_corr == 0
    set(handles.edit15,'Enable','off');
    set(handles.uipanel1,'Visible','off');
else
    set(handles.edit15,'Enable','on');
    set(handles.uipanel1,'Visible','on');    
    if handles.par_BG_corr == 2
        set(handles.radiobutton3 ,'Value',1);
    end
end
if handles.par_num_PIE == 4
    set(handles.edit19,'enable','on');
end
if handles.par_num_PIE > 2
    set(handles.edit20,'enable','on');
end
if handles.par_num_PIE > 1
    set(handles.edit18,'enable','on');
end
% UIWAIT makes Set_Params wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Set_Params_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
close;


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
tmp = round(str2double(get(hObject,'String')));
tmp = max([300 tmp]);
tmp = min([700 tmp]);
handles.par_lam_start = tmp;
guidata(hObject, handles);
set(hObject,'String',sprintf('%d',round(handles.par_lam_start)));

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
tmp = round(str2double(get(hObject,'String')));
tmp = max([  1 tmp]);
tmp = min([100 tmp]);
handles.par_lam_step = tmp;
guidata(hObject, handles);
set(hObject,'String',sprintf('%d',round(handles.par_lam_step)));


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
tmp = round(str2double(get(hObject,'String')));
tmp = max([ 1 tmp]);
tmp = min([ 4 tmp]);
handles.par_num_PIE = tmp;

set(handles.edit18,'enable','off');
set(handles.edit19,'enable','off');
set(handles.edit20,'enable','off');

if tmp == 4
    set(handles.edit19,'enable','on');
end
if tmp > 2
    set(handles.edit20,'enable','on');
end
if tmp > 1
    set(handles.edit18,'enable','on');
end

guidata(hObject, handles);
set(hObject,'String',sprintf('%d',round(handles.par_num_PIE)));


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
tmp = round(str2double(get(hObject,'String')));
tmp = max([handles.par_wlam_min+1 tmp]);
tmp = min([100 tmp]);
handles.par_wlam_max = tmp; 
guidata(hObject, handles);

set(hObject,'String',sprintf('%d',round(handles.par_wlam_max)));


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double
tmp = round(str2double(get(hObject,'String')));
tmp = min([handles.par_wlam_max-1 tmp]);
tmp = max([0 tmp]);
handles.par_wlam_min = tmp;
guidata(hObject, handles);

set(hObject,'String',sprintf('%d',round(handles.par_wlam_min)));


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double
tmp = round(10*str2double(get(hObject,'String')))/10;
tmp = min([handles.par_tau_max-1 tmp]);
tmp = max([0 tmp]);
handles.par_tau_min = tmp; 
guidata(hObject, handles);

set(hObject,'String',sprintf('%4.1f',(handles.par_tau_min)));


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double

tmp = round(10*str2double(get(hObject,'String')))/10;
tmp = max([handles.par_tau_min+1 tmp]);
tmp = min([20 tmp]);
handles.par_tau_max = tmp; 
guidata(hObject, handles);

set(hObject,'String',sprintf('%4.1f',(handles.par_tau_max)));


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double

tmp = round(str2double(get(hObject,'String')));
tmp = min([handles.par_cum2_max-3 tmp]);
tmp = max([0 tmp]);
handles.par_cum2_min = tmp; 
guidata(hObject, handles);

set(hObject,'String',sprintf('%d',round(handles.par_cum2_min)));


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double

tmp = round(str2double(get(hObject,'String')));
tmp = min([handles.par_clam_max-25 tmp]);
tmp = max([0 tmp]);
handles.par_clam_min = tmp; 
guidata(hObject, handles);

set(hObject,'String',sprintf('%d',round(handles.par_clam_min)));


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double

tmp = round(str2double(get(hObject,'String')));
tmp = max([handles.par_clam_min+25 tmp]);
tmp = min([900 tmp]);
handles.par_clam_max = tmp; 
guidata(hObject, handles);

set(hObject,'String',sprintf('%d',round(handles.par_clam_max)));


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double

tmp = round(str2double(get(hObject,'String')));
tmp = max([handles.par_cum2_min+3 tmp]);
tmp = min([500 tmp]);
handles.par_cum2_max = tmp; 
guidata(hObject, handles);

set(hObject,'String',sprintf('%d',round(handles.par_cum2_max)));


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double
tmp = round(str2double(get(hObject,'String')));
tmp = max([0 tmp]);
tmp = min([500 tmp]);
handles.par_Threshold = tmp; 
guidata(hObject, handles);

set(hObject,'String',sprintf('%d',round(handles.par_Threshold)));


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double
tmp = round(str2double(get(hObject,'String'))*100)/100;
tmp = max([0.01 tmp]);
tmp = min([100 tmp]);
handles.par_pixwidth = tmp; 
guidata(hObject, handles);

set(hObject,'String',sprintf('%4.2f',(handles.par_pixwidth)));


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double
tmp = round(str2double(get(hObject,'String')));
tmp = max([1 tmp]);
tmp = min([10 tmp]);
handles.par_binning = tmp; 
guidata(hObject, handles);

set(hObject,'String',sprintf('%d',round(handles.par_binning)));


% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double
tmp = round(100*str2double(get(hObject,'String')))/100;
tmp = max([0 tmp]);
tmp = min([10 tmp]);
handles.par_BG_AP = tmp; 
guidata(hObject, handles);

set(hObject,'String',sprintf('%4.2f',handles.par_BG_AP));


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tmp = get(hObject,'Value');
tmp = max([0 tmp]);
tmp = min([1 tmp]);
guidata(hObject, handles);
if tmp == 0
    set(handles.edit15,'Enable','off');
    set(handles.uipanel1,'Visible','off');
else
    set(handles.edit15,'Enable','on');
    set(handles.uipanel1,'Visible','on');
end
tmp1 = get(handles.uipanel1,'Children');
k = 1;
while get(tmp1(k),'Value') == 0
    k = k+1;
end
if strcmp(get(tmp1(k),'String'),'Individual determination')
   handles.par_BG_corr = tmp*2;
else
   handles.par_BG_corr = tmp;
end
set(hObject,'Value',tmp);
guidata(hObject, handles);

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);



% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%delete('Get_Params.m');

%fid = fopen('Get_Params.m','w');

lam_start = handles.par_lam_start;
lam_step  = handles.par_lam_step;

tau_min  = handles.par_tau_min;
tau_max  = handles.par_tau_max;

cum2_min = handles.par_cum2_min;
cum2_max = handles.par_cum2_max;

clam_min = handles.par_clam_min;
clam_max = handles.par_clam_max;

wlam_min = handles.par_wlam_min;
wlam_max = handles.par_wlam_max;

num_PIE  = handles.par_num_PIE;

pixwidth = handles.par_pixwidth;

Threshold = handles.par_Threshold;
binning   = handles.par_binning;

lamPIE    = handles.par_lamPIE;

maxROI   = 16;
ROIcolor = [1 0 0;       ...
            0 1 0;       ...   
            0 0 1;       ...
            1 1 0;       ...
            1 0 1;       ...
            0 1 1;       ...
            1 1 1;       ...
            0 0.5 0;     ...
            0.5 0.5 0;   ...
            0.3 0.3 0.3; ...
            0.0 0.5 0.5; ...
            0.5 0.0 0.5; ...
            0.7 0.5 0.3; ...
            0.3 0.7 0.5; ...
            0.5 0.3 0.7; ...
            0.3 0.5 0.7];
        
ROIlist = zeros(16,1);

BG_corr = handles.par_BG_corr;
BG_AP   = handles.par_BG_AP/100;

BG_file = '';

if BG_corr == 1
    
    [FileName,PathName] = uigetfile({'*.ht3;*.pt3','t3r-files';'*.*','All Files'}, 'MultiSelect', 'off');
    
    if ~isempty(FileName)
        BG_file = [PathName FileName];    
    else
        warndlg('No file selected!')
        BG_corr = 2;
    end
end

figFLIM  = 10;
figMA    = 20;
figfFLIM = 30;
figone   = 40;
figPA    = 50;

save('Get_Params.mat','lam_start','lam_step', 'tau_max', 'tau_min', 'cum2_min', 'cum2_max', ...
                      'clam_min', 'clam_max', 'wlam_min', 'wlam_max', 'num_PIE', ...
                      'Threshold', 'binning', 'maxROI', 'ROIcolor', 'ROIlist', ...
                      'pixwidth', 'BG_corr', 'BG_file', 'BG_AP', 'lamPIE',  ...  
                      'figFLIM', 'figfFLIM', 'figMA', 'figone', 'figPA');

uiresume(handles.figure1);


function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double
tmp = round(str2double(get(hObject,'String')));
tmp = max([ 300 tmp]);
tmp = min([1000 tmp]);
handles.par_lamPIE(1) = tmp; 
guidata(hObject, handles);

set(hObject,'String',sprintf('%d',round(handles.par_lamPIE(1))));


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double
tmp = round(str2double(get(hObject,'String')));
tmp = max([ 300 tmp]);
tmp = min([1000 tmp]);
handles.par_lamPIE(2) = tmp; 
guidata(hObject, handles);

set(hObject,'String',sprintf('%d',round(handles.par_lamPIE(2))));


% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double
tmp = round(str2double(get(hObject,'String')));
tmp = max([ 300 tmp]);
tmp = min([1000 tmp]);
handles.par_lamPIE(4) = tmp; 
guidata(hObject, handles);

set(hObject,'String',sprintf('%d',round(handles.par_lamPIE(4))));


% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double
tmp = round(str2double(get(hObject,'String')));
tmp = max([ 300 tmp]);
tmp = min([1000 tmp]);
handles.par_lamPIE(3) = tmp; 
guidata(hObject, handles);

set(hObject,'String',sprintf('%d',round(handles.par_lamPIE(3))));


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uipanel1.
function uipanel1_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel1 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

if strcmp(get(eventdata.NewValue,'String'),'Individual determination')
   handles.par_BG_corr = 2;
else
   handles.par_BG_corr = 1;
end
guidata(hObject, handles);

