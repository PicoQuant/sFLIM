function varargout = Results_4(varargin)
% RESULTS_4 M-file for Results_4.fig
%      RESULTS_4, by itself, creates a new RESULTS_4 or raises the existing
%      singleton*.
%
%      H = RESULTS_4 returns the handle to a new RESULTS_4 or the handle to
%      the existing singleton*.
%
%      RESULTS_4('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RESULTS_4.M with the given input arguments.
%
%      RESULTS_4('Property','Value',...) creates a new RESULTS_1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Results_1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Results_1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Results_4

% Last Modified by GUIDE v2.5 15-Sep-2011 13:24:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Results_4_OpeningFcn, ...
    'gui_OutputFcn',  @Results_4_OutputFcn, ...
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


% --- Executes just before Results_4 is made visible.
function Results_4_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Results_4 (see VARARGIN)

% Choose default command line output for Results_4
handles.output = hObject;

set(handles.Update_button,'Enable','off');
set(handles.NewROI_button,'Enable','off');
set(handles.popupmenu1,'Enable','off');
set(handles.popupmenu2,'Enable','off');
set(handles.ClearROI_button,'Enable','off');
set(handles.SaveROI_button,'Enable','off');
set(handles.SavePattern_button,'Enable','off');
set(handles.ClearPattern_button,'Enable','off');
set(handles.ShowPattern_button,'Enable','off');
set(handles.FindAmplitudes_button,'Enable','off');
set(handles.PlotResults_button,'Enable','off');
set(handles.SaveResults_button,'Enable','off');
set(handles.Merge_button,'Enable','off');
set(handles.Intersect_button,'Enable','off');
set(handles.Add_Similar_button,'Enable','off');
set(handles.FitPattern_button,'Enable','off');

% Update handles structure
guidata(hObject, handles);

if nargin < 2
    errordlg('No filename given','File Load Error')
elseif (length(varargin) == 2 && ...
        strcmpi(varargin{1},'file') && ...
        (2 == exist([varargin{2}(1:end-4) '_Moments.mat'],'file')))
    handles.filename      = cell2mat(varargin(2));
else
    errordlg('File Not Found','File Load Error')
end

tmp = strfind(handles.filename,'\');
handles.pathname = handles.filename(1:tmp(end));
handles.filename = handles.filename(tmp(end)+1:end);

load Get_Params;
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
handles.par_lam_PIE   = lamPIE;
handles.par_Threshold = Threshold;
handles.par_binning   = binning;
handles.ROIcolor      = ROIcolor;
handles.maxROI        = maxROI;
handles.figFLIM       = figFLIM;
handles.figMA         = figMA;
handles.figPA         = figPA;
handles.figfFLIM      = figfFLIM;
handles.figone        = figone;

handles.pulse         = 1;
handles.match         = 0;
handles.ROIlist       = zeros(maxROI, num_PIE);

[handles.head, handles.tag, handles.tcspc, handles.IRF, handles.timname, handles.AP, handles.DC] = load_data([handles.pathname handles.filename]);

nx      = handles.head.ImgHdr.PixX;
ny      = handles.head.ImgHdr.PixY;
tx      = floor(nx/binning);
ty      = floor(ny/binning);
nch     = size(handles.tcspc,1);
nbin    = size(handles.tcspc,2);
num_PIE = size(handles.tcspc,3);

if binning > 1
    tmp = handles.tag(1:binning*tx,:,:,:);
    tmp = reshape(tmp, binning, tx, ny, nch, num_PIE);
    tmp = squeeze(mean(tmp,1));
    tmp = tmp(:,1:binning*ty,:,:);
    tmp = reshape(tmp, tx, binning, ty, nch, num_PIE);
    handles.tag = squeeze(mean(tmp,2));
end

if handles.head.ImgHdr.PixelSize ==0
    handles.head.ImgHdr.PixelSize = pixwidth;
end

handles.plt_tau_min   = tau_min;
handles.plt_tau_max   = tau_max;
handles.plt_cum2_min  = cum2_min;
handles.plt_cum2_max  = cum2_max;
handles.plt_clam_min  = clam_min;
handles.plt_clam_max  = clam_max;
handles.plt_wlam_min  = wlam_min;
handles.plt_wlam_max  = wlam_max;
handles.plt_x_min     = handles.head.ImgHdr.X0;
handles.plt_x_max     = handles.head.ImgHdr.X0+handles.head.ImgHdr.PixX*handles.head.ImgHdr.PixelSize;
handles.plt_y_min     = handles.head.ImgHdr.Y0;
handles.plt_y_max     = handles.head.ImgHdr.Y0+handles.head.ImgHdr.PixY*handles.head.ImgHdr.PixelSize;

handles.channel       = zeros(handles.par_num_PIE, size(handles.tcspc,1));

load([handles.pathname handles.filename(1:end-4) '_Moments.mat'],'tav');

handles.tav = tav;
handles.pat = zeros(1, size(handles.tcspc,1), size(handles.tcspc,2), handles.par_num_PIE);

handles.s_ind = uint16(zeros(size(handles.tag,1)*size(handles.tag,2),handles.par_num_PIE));

guidata(hObject, handles);

load('c_map.mat','map');
colormap(map);

s = {};

set(handles.Pattern_List,'String',s);

for n = 1:num_PIE
    s(n) = {sprintf('pulse %d',n)};
end;

set(handles.popupmenu1, 'String',s);

lambda = (handles.par_lam_start+handles.par_lam_step.*(0:size(handles.tcspc,1)-1));
tmp    = repmat(handles.par_lam_PIE(1:handles.par_num_PIE)', [1 size(handles.tcspc,1)])+handles.par_lam_step;

handles.channel = repmat(lambda,[handles.par_num_PIE 1])>tmp;

s = {};

for n = 1:size(handles.tcspc,1);
    if handles.channel(handles.pulse,n)
        s(n) = {sprintf('[X] %03d nm',lambda(n))};
    else
        s(n) = {sprintf('[ ] %03d nm',lambda(n))};
    end
end;

set(handles.popupmenu2, 'String',s);

set(handles.axes4,'visible', 'on');
set(handles.axes4, 'PlotBoxAspectRatio',[1 1 1], ...
    'Box','on', ...
    'XDir','normal', ...
    'YDir','normal', ...
    'FontSize',9,...
    'Color',[50 50 50]./255);
ylabel(handles.axes4, '< \tau > / ns','FontWeight','bold','FontSize',10,'FontAngle','italic');
xlabel(handles.axes4, '\lambda_0 / nm','FontWeight','bold','FontSize',10,'FontAngle','italic');
axis(handles.axes4, [tau_min tau_max clam_min clam_max]);

set(handles.axes1, 'DataAspectRatio', [1,1,1], ...
    'PlotBoxAspectRatio',[1 1 1], ...
    'FontSize',9,...
    'XDir','normal', ...
    'YDir','reverse');
xlabel(handles.axes2, 'x / µm','FontWeight','bold','FontSize',10,'FontAngle','italic');
ylabel(handles.axes2, 'y / µm','FontWeight','bold','FontSize',10,'FontAngle','italic');

set(handles.axes2, 'PlotBoxAspectRatio',[1 1 1], ...
    'Box','on', ...
    'XDir','normal', ...
    'YDir','normal', ...
    'FontSize',9,...
    'Color',[50 50 50]./255);
xlabel(handles.axes2, '\Delta \lambda / nm','FontWeight','bold','FontSize',10,'FontAngle','italic');
ylabel(handles.axes2, '\lambda_0 / nm','FontWeight','bold','FontSize',10,'FontAngle','italic');
axis(handles.axes2, [wlam_min wlam_max clam_min clam_max]);


[handles.M_X, handles.M_Y, handles.Z, handles.ZZ] = compute_data(handles);
compute_fFLIM(handles);
guidata(hObject, handles);

plot_graphs(handles)

pause(0.05)
set(handles.Update_button,'Enable','on');
set(handles.NewROI_button,'Enable','on');
set(handles.LoadROI_button,'Enable','on');
set(handles.LoadPattern_button,'Enable','on');
set(handles.popupmenu1,'Enable','on');
set(handles.popupmenu2,'Enable','on');


% --- Outputs from this function are returned to the command line.
function varargout = Results_4_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
    set(hObject,'String',{});
end


% --- Executes during object creation, after setting all properties.
function Pattern_List_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pattern_List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
    set(hObject,'String',{});
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tmp = get(hObject,'Value');
if tmp ~= handles.pulse
    handles.pulse = tmp;
    guidata(hObject, handles);
end;


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2

lambda = (handles.par_lam_start+handles.par_lam_step.*(0:size(handles.tcspc,1)-1));

% for n = 1:size(handles.tcspc,1);
%    s(n) = {sprintf('[X] %03d nm',lambda(n))};
% end;

tmp = get(hObject,'Value');
s   = get(hObject, 'String');

handles.channel(handles.pulse,tmp) = ~handles.channel(handles.pulse,tmp);

if handles.channel(handles.pulse,tmp)
    s(tmp) = {sprintf('[X] %03d nm',lambda(tmp))};
else
    s(tmp) = {sprintf('[ ] %03d nm',lambda(tmp))};
end;
set(hObject, 'String',s);
guidata(hObject, handles);


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

index_selected = get(hObject, 'Value');
if numel(index_selected)>1
    set(handles.Intersect_button,'Enable','on');
    set(handles.Merge_button,'Enable','on');
else
    set(handles.Intersect_button,'Enable','off');
    set(handles.Merge_button,'Enable','off');
end

if numel(index_selected)>0
    set(handles.ClearROI_button,'Enable','on');
    set(handles.SaveROI_button,'Enable','on');
    set(handles.UseROI_button,'Enable','on');
    set(handles.Add_Similar_button,'Enable','on');
else
    set(handles.ClearROI_button,'Enable','off');
    set(handles.SaveROI_button,'Enable','off');
    set(handles.UseROI_button,'Enable','off');
    set(handles.Add_Similar_button,'Enable','off');
end
guidata(hObject, handles);


% --- Executes on selection change in Pattern_List.
function Pattern_List_Callback(hObject, eventdata, handles)
% hObject    handle to Pattern_List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
index_selected = get(hObject, 'Value');
if numel(index_selected)>0
    set(handles.ClearPattern_button,'Enable','on');
    set(handles.SavePattern_button,'Enable','on');
    set(handles.ShowPattern_button,'Enable','on');
    set(handles.FindAmplitudes_button,'Enable','on');
    set(handles.FitPattern_button,'Enable','on');
else
    set(handles.ClearPattern_button,'Enable','off');
    set(handles.SavePattern_button,'Enable','off');
    set(handles.ShowPattern_button,'Enable','off');
    set(handles.FindAmplitudes_button,'Enable','off');    
end
guidata(hObject, handles);


% --- Executes on button press in Update_button.
function Update_button_Callback(hObject, eventdata, handles)
% hObject    handle to Update_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.Update_button,'String','Busy');
set(handles.Update_button,'Enable','off');
pause(0.05)
guidata(hObject, handles);
[handles.M_X, handles.M_Y, handles.Z, handles.ZZ] = compute_data(handles);
lambda = (handles.par_lam_start+handles.par_lam_step.*(0:size(handles.tcspc,1)-1));
s = {};
for n = 1:size(handles.tcspc,1);
    if handles.channel(handles.pulse,n)
        s(n) = {sprintf('[X] %03d nm',lambda(n))};
    else
        s(n) = {sprintf('[ ] %03d nm',lambda(n))};
    end
end;
set(handles.popupmenu2, 'String',s);


guidata(hObject, handles);
plot_graphs(handles)
s = {};
set(handles.listbox2, 'Value',[]);
if sum(handles.ROIlist(:, handles.pulse))>0
    tst = 1;
    for n = 1:handles.maxROI
        if handles.ROIlist(n, handles.pulse)==1
            s(tst) = {sprintf('ROI %d_%d',n,handles.pulse)};
            tst = tst+1;
        end
    end
    set(handles.listbox2,'String',s);
    set(handles.listbox2,'Enable','on');
    set(handles.ClearROI_button,'Enable','on');

    update_fFLIM(handles)
else
    set(handles.listbox2,'String',s);
    set(handles.listbox2,'Enable','off');
    set(handles.ClearROI_button,'Enable','off');
end

set(handles.Update_button,'String','Update');
set(handles.Update_button,'Enable','on');

% --- Executes on button press in zoom_button.
function zoom_button_Callback(hObject, eventdata, handles)
% hObject    handle to zoom_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes1);

[xi,yi,button] = ginput(1);

if button==1

    pt1     = get(handles.axes1, 'CurrentPoint');
    pt2     = get(handles.axes2, 'CurrentPoint');
    pt3     = get(handles.axes3, 'CurrentPoint');
    pt4     = get(handles.axes4, 'CurrentPoint');

    t1 = xi-[pt1(1,1) pt2(1,1) pt3(1,1) pt4(1,1)];
    t2 = yi-[pt1(1,2) pt2(1,2) pt3(1,2) pt4(1,2)];
    p  = find((t1+t2)==0);

    if isempty(p)
        p = abs(t1+t2);
        p = find(p==min(p));
    end

    if p==1
        axes(handles.axes1);
    elseif p==2
        axes(handles.axes2);
    elseif p==3
        axes(handles.axes3);
    else
        axes(handles.axes4);
    end

    hold on

    plot(xi,yi,'+','Color',[1 1 1])
    xy   = [xi;yi];

    [xi,yi,button] = ginput(1);
    if button==1
        plot(xi,yi,'+','Color',[1 1 1])
        xy(:,2) = [xi;yi];
        hold off

        x_min = min(xy(1,:));
        x_max = max(xy(1,:));
        y_min = min(xy(2,:));
        y_max = max(xy(2,:));

        if p==1         % Intensity image
            handles.plt_x_min = x_min;
            handles.plt_x_max = x_max;
            handles.plt_y_min = y_min;
            handles.plt_y_max = y_max;
        elseif p==2     % wavelength diagram
            handles.plt_wlam_min = x_min;
            handles.plt_wlam_max = x_max;
            handles.plt_clam_min = y_min;
            handles.plt_clam_max = y_max;
        elseif p==3      % moments plot
            handles.plt_tau_min = x_min;
            handles.plt_tau_max = x_max;
            handles.plt_cum2_min = y_min;
            handles.plt_cum2_max = y_max;
        elseif p==4      % tau-wavelength diagram
            handles.plt_tau_min = x_min;
            handles.plt_tau_max = x_max;
            handles.plt_clam_min = y_min;
            handles.plt_clam_max = y_max;
        end;
    else
        if p==1         % Intensity image
            handles.plt_x_min = handles.head.ImgHdr.X0;
            handles.plt_x_max = handles.head.ImgHdr.X0+handles.head.ImgHdr.PixX*handles.head.ImgHdr.PixelSize;
            handles.plt_y_min = handles.head.ImgHdr.Y0;
            handles.plt_y_max = handles.head.ImgHdr.Y0+handles.head.ImgHdr.PixY*handles.head.ImgHdr.PixelSize;
        elseif p==2     % wavelength diagram
            handles.plt_wlam_min = handles.par_wlam_min;
            handles.plt_wlam_max = handles.par_wlam_max;
            handles.plt_clam_min = handles.par_clam_min;
            handles.plt_clam_max = handles.par_clam_max;
        elseif p==3      % moments plot
            handles.plt_tau_min  = handles.par_tau_min;
            handles.plt_tau_max  = handles.par_tau_max;
            handles.plt_cum2_min = handles.par_cum2_min;
            handles.plt_cum2_max = handles.par_cum2_max;
        elseif p==4      % tau-wavelength diagram
            handles.plt_tau_min  = handles.par_tau_min;
            handles.plt_tau_max  = handles.par_tau_max;
            handles.plt_clam_min = handles.par_clam_min;
            handles.plt_clam_max = handles.par_clam_max;
        end;
    end
else
    handles.plt_tau_min   = handles.par_tau_min;
    handles.plt_tau_max   = handles.par_tau_max;
    handles.plt_cum2_min  = handles.par_cum2_min;
    handles.plt_cum2_max  = handles.par_cum2_max;
    handles.plt_clam_min  = handles.par_clam_min;
    handles.plt_clam_max  = handles.par_clam_max;
    handles.plt_wlam_min  = handles.par_wlam_min;
    handles.plt_wlam_max  = handles.par_wlam_max;
    handles.plt_x_min     = handles.head.ImgHdr.X0;
    handles.plt_x_max     = handles.head.ImgHdr.X0+handles.head.ImgHdr.PixX*handles.head.ImgHdr.PixelSize;
    handles.plt_y_min     = handles.head.ImgHdr.Y0;
    handles.plt_y_max     = handles.head.ImgHdr.Y0+handles.head.ImgHdr.PixY*handles.head.ImgHdr.PixelSize;
end
 
guidata(hObject, handles)
plot_graphs(handles);

% --- Executes on button press in NewROI_button.
function NewROI_button_Callback(hObject, eventdata, handles)
% hObject    handle to NewROI_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.NewROI_button,'String','Click');
set(handles.NewROI_button,'Enable','off');
pause(0.05)

tmp = 1;
while handles.ROIlist(tmp, handles.pulse)==1
    tmp = tmp +1;
end;

if tmp <= handles.maxROI
    handles.ROIlist(tmp, handles.pulse) = 1;

    guidata(hObject, handles);

    handles.s_ind = select_data(handles, tmp);

    guidata(hObject, handles);

    plot_graphs(handles)

    s = {};
    tst = 1;
    for n = 1:handles.maxROI
        if handles.ROIlist(n, handles.pulse)==1
            s(tst) = {sprintf('ROI %d_%d',n,handles.pulse)};
            tst = tst+1;
        end
    end
    set(handles.listbox2,'String',s);
    set(handles.listbox2,'Enable','on');
    set(handles.ClearROI_button,'Enable','on');
end;

set(handles.NewROI_button,'String','new ROI');
if tmp<handles.maxROI
    set(handles.NewROI_button,'Enable','on');
    set(handles.LoadROI_button,'Enable','on');
else
    set(handles.LoadROI_button,'Enable','off');
end
guidata(hObject, handles)


% --- Executes on button press in ClearROI_button.
function ClearROI_button_Callback(hObject, eventdata, handles)  % clear ROI
% hObject    handle to ClearROI_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tmp = sum(handles.ROIlist(:, handles.pulse));

if tmp>0
    set(handles.ClearROI_button,'String','Busy');
    set(handles.ClearROI_button,'Enable','off');
    pause(0.05);

    list_entries = get(handles.listbox2, 'String');
    index_selected = get(handles.listbox2, 'Value');
    pat_list_entries = get(handles.Pattern_List, 'String');

    for n = 1:numel(index_selected)

        s   = list_entries{index_selected(n)};
        ROI = str2double(s(end-2));

        handles.s_ind(:,handles.pulse) = bitset(handles.s_ind(:,handles.pulse), ROI, 0);
        handles.ROIlist(ROI, handles.pulse) = 0;
        
        tst =  strcmp(pat_list_entries,s);
        handles.pat(tst,:,:,:) = [];
        pat_list_entries(tst) = [];
    end

    set(handles.Pattern_List,'String',pat_list_entries);

    guidata(hObject, handles);
    
    s = {};

    tst = 1;
    for n = 1:handles.maxROI
        if handles.ROIlist(n, handles.pulse)==1
            s(tst) = {sprintf('ROI %d_%d',n,handles.pulse)};
            tst = tst+1;
        end
    end
    set(handles.listbox2,'Value',1);
    set(handles.listbox2,'String',s);
    listbox2_Callback(handles.listbox2, [], handles);

    plot_graphs(handles)
    %     update_fFLIM(handles)

    set(handles.ClearROI_button,'String','clear ROI');
    if tst>1
        set(handles.ClearROI_button,'Enable','on');
    end
    if tst<=handles.maxROI
        set(handles.NewROI_button,'Enable','on');
    end
end
guidata(hObject, handles);



% --- Executes on button press in LoadROI_button.
function LoadROI_button_Callback(hObject, eventdata, handles)  % Load ROI
% hObject    handle to LoadROI_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(hObject,'Enable','off');

ROInum = 1;
while handles.ROIlist(ROInum, handles.pulse)==1
    ROInum = ROInum +1;
end;

if ROInum <= handles.maxROI
    handles.ROIlist(ROInum, handles.pulse)=1;
    [FileName,PathName,FilterIndex] = uigetfile({'*.roi','ROI Files';...
        [handles.pathname '*.*'],'All Files'},'Load ROI', [handles.pathname 'newroi.roi']);
    if FileName ~= 0
        load([PathName FileName],'ind', '-mat');

        handles.s_ind(ind, handles.pulse) = bitset(handles.s_ind(ind, handles.pulse), ROInum,1);

        guidata(hObject, handles);

        s = {};
        tst = 1;
        for n = 1:handles.maxROI
            if handles.ROIlist(n, handles.pulse)==1
                s(tst) = {sprintf('ROI %d_%d',n,handles.pulse)};
                tst = tst+1;
            end
        end
        set(handles.listbox2,'Value',1);
        set(handles.listbox2,'String',s);
        listbox2_Callback(handles.listbox2, [], handles);

        plot_graphs(handles)
        %          update_fFLIM(handles)
    end
end

set(hObject,'String','Load ROI');
if ROInum>=handles.maxROI
    set(handles.NewROI_button,'Enable','off');
end
guidata(hObject, handles);


% --- Executes on button press in SaveROI_button.
function SaveROI_button_Callback(hObject, eventdata, handles)  % Save ROI
% hObject    handle to SaveROI_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(hObject,'Enable','off');

list_entries = get(handles.listbox2, 'String');
index_selected = get(handles.listbox2, 'Value');

for n = 1:numel(index_selected)

    s   = list_entries{index_selected(n)};
    ROI = str2double(s(end-2));

    ind = bitget(handles.s_ind(:, handles.pulse), ROI) == 1;

    [FileName,PathName] = uiputfile({'*.roi','ROI Files';...
        '*.*','All Files' },'Save ROI',...
        [handles.pathname 'newroi.roi']);
    if FileName ~= 0
        save([PathName FileName],'ind');
    end
end


% --- Executes on button press in Merge_button.
function Merge_button_Callback(hObject, eventdata, handles)
% hObject    handle to Merge_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(hObject,'String','Busy');
set(hObject,'Enable','off');

list_entries   = get(handles.listbox2, 'String');
index_selected = get(handles.listbox2, 'Value');

if numel(index_selected) > 1

    s   = list_entries{index_selected(1)};
    ROI = str2double(s(end-2));
    
    for n = 2:numel(index_selected)

        s = list_entries{index_selected(n)};
        k = str2double(s(end-2));

        handles.ROIlist(k, handles.pulse) = 0;

        ind = bitget(handles.s_ind(:, handles.pulse), k) == 1;
        handles.s_ind(ind, handles.pulse) = bitset(handles.s_ind(ind, handles.pulse), k, 0);
        handles.s_ind(ind, handles.pulse) = bitset(handles.s_ind(ind, handles.pulse), ROI, 1);
    end
    guidata(hObject, handles);
end

s = {};
tst = 1;
for n = 1:handles.maxROI
    if handles.ROIlist(n, handles.pulse)==1
        s(tst) = {sprintf('ROI %d_%d',n,handles.pulse)};
        tst = tst+1;
    end
end
set(handles.listbox2,'Value',1);
set(handles.listbox2,'String',s);
listbox2_Callback(handles.listbox2, [], handles);

plot_graphs(handles)
update_fFLIM(handles)
set(hObject,'String','Merge ROIs');


% --- Executes on button press in Intersect_button.
function Intersect_button_Callback(hObject, eventdata, handles)
% hObject    handle to Intersect_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(hObject,'String','Busy');
set(hObject,'Enable','off');

list_entries = get(handles.listbox2, 'String');
index_selected = get(handles.listbox2, 'Value');

if numel(index_selected) > 1

    s   = list_entries{index_selected(1)};
    ROI = str2double(s(end-2));

    ind = bitget(handles.s_ind(:, handles.pulse), ROI) == 1;

    handles.s_ind(:, handles.pulse) = bitset(handles.s_ind(:, handles.pulse), ROI, 0);

    for n = 2:numel(index_selected)

        s = list_entries{index_selected(n)};
        k = str2double(s(end-2));

        ind = ind & bitget(handles.s_ind(:, handles.pulse), k) == 1;

        handles.s_ind(:, handles.pulse) = bitset(handles.s_ind(:, handles.pulse), k, 0);
        handles.ROIlist(k, handles.pulse) = 0;
        
    end

    handles.s_ind(ind, handles.pulse) = bitset(handles.s_ind(ind, handles.pulse), ROI, 1);
   
    guidata(hObject, handles);
end

s = {};
tst = 1;
for n = 1:handles.maxROI
    if handles.ROIlist(n, handles.pulse)==1
        s(tst) = {sprintf('ROI %d_%d',n,handles.pulse)};
        tst = tst+1;
    end
end
set(handles.listbox2,'Value',1);
set(handles.listbox2,'String',s);
listbox2_Callback(handles.listbox2, [], handles);

plot_graphs(handles)
update_fFLIM(handles)
set(hObject,'String','Intersection');


% --- Executes on button press in Add_Similar_button.
function Add_Similar_button_Callback(hObject, eventdata, handles)
% hObject    handle to Add_Similar_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(hObject,'String','Busy');
set(hObject,'Enable','off');

list_entries = get(handles.listbox2, 'String');
index_selected = get(handles.listbox2, 'Value');

if numel(index_selected) >= 1
      
        ROInum = 1;
        while handles.ROIlist(ROInum, handles.pulse)==1
            ROInum = ROInum +1;
        end;
        
        if ROInum <= handles.maxROI

            handles.ROIlist(ROInum, handles.pulse) = 1;
            
            s   = list_entries{index_selected(1)};
            ROI = str2double(s(end-2));

            ind = bitget(handles.s_ind(:, handles.pulse), ROI) == 1;

            ind = ind & ~isnan(handles.M_X) & ~isnan(handles.M_Y) & ~isnan(handles.Z) & ~isnan(handles.ZZ);
            tmp = squeeze(reshape(handles.tag(:,:,:,handles.pulse),size(handles.tag,1)*size(handles.tag,2),1,size(handles.tag,3)));

            for n = 1:size(tmp,2)
                tmp(:,n) = tmp(:,n).*handles.channel(n);
            end
            intens   = sum(tmp,2);

            ind1 = handles.M_X >= min(handles.M_X(ind)) & handles.M_X <= max(handles.M_X(ind)) & ...
                   handles.M_Y >= min(handles.M_Y(ind)) & handles.M_Y <= max(handles.M_Y(ind)) & ...
                   handles.Z   >= min(handles.Z(ind))   & handles.Z   <= max(handles.Z(ind))   & ...
                   handles.ZZ  >= min(handles.ZZ(ind))  & handles.ZZ  <= max(handles.ZZ(ind))  & ...
                   intens      >= min(intens(ind))      & intens      <= max(intens(ind));

            handles.s_ind(ind1, handles.pulse) = bitset(handles.s_ind(ind1, handles.pulse), ROInum, 1);
        end
        guidata(hObject, handles);
end

s = {};
tst = 1;
for n = 1:handles.maxROI
    if handles.ROIlist(n, handles.pulse)==1
        s(tst) = {sprintf('ROI %d_%d',n,handles.pulse)};
        tst = tst+1;
    end
end

set(handles.listbox2,'Value',1);
set(handles.listbox2,'String',s);
listbox2_Callback(handles.listbox2, [], handles);

plot_graphs(handles)
update_fFLIM(handles)
set(hObject,'String','Add similar');


% --- Executes on button press in UseROI_button.
function UseROI_button_Callback(hObject, eventdata, handles)
% hObject    handle to UseROI_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.UseROI_button,'String','Busy');
set(handles.UseROI_button,'Enable','off');
pause(0.05);

list_entries = get(handles.listbox2, 'String');
index_selected = get(handles.listbox2, 'Value');

binning = handles.par_binning;
nx      = handles.head.ImgHdr.PixX;
ny      = handles.head.ImgHdr.PixY;
tx      = floor(nx/binning);
ty      = floor(ny/binning);
nch     = size(handles.tcspc,1);
nbin    = size(handles.tcspc,2);
num_PIE = size(handles.tcspc,3);

load([handles.pathname handles.timname], 'stim')

if binning > 1
    tmp = stim(1:binning*tx,:,:,:,:);
    tmp = reshape(tmp, binning, tx, ny, nch, nbin, num_PIE);
    tmp = squeeze(mean(tmp,1));
    tmp = tmp(:,1:binning*ty,:,:);
    tmp = reshape(tmp, tx, binning, ty, nch, nbin, num_PIE);
    stim = permute(mean(tmp,2), [1 3 4 5 6 2]);
end

stim = reshape(stim, tx*ty,nch, nbin, handles.par_num_PIE);

pat_list_entries = get(handles.Pattern_List, 'String');

for n = 1:numel(index_selected)

    s   = cell2mat(list_entries(index_selected(n)));
    ROI = str2double(s(end-2));
    if ROI == 0
        ROI = 10;
    end
    p   = numel(pat_list_entries)+1 ;
    for k = 1:numel(pat_list_entries)
        if strcmp(pat_list_entries(k), s)
            p = k;
        end
    end
    pat_list_entries(p) = {s};
    
    ind = bitget(handles.s_ind(:, handles.pulse), ROI) == 1;

    pat = squeeze(sum(stim(ind,:,:,:),1));
    handles.pat(p,:,:,:) = pat./sum(sum(sum(pat)));
    
end

set(handles.Pattern_List,'Value',1);
set(handles.Pattern_List,'String',pat_list_entries);
set(handles.UseROI_button,'String','Use ROI');

guidata(hObject, handles);

pause(0.05);


% --- Executes on button press in LoadPattern_button.
function LoadPattern_button_Callback(hObject, eventdata, handles)
% hObject    handle to LoadPattern_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%set(hObject,'Enable','off');

[FileName,PathName] = uigetfile({'*.pat','Pattern';'*.*','All Files'}, ...
         'Load pattern',[handles.pathname 'new.pat'],'MultiSelect', 'on');
    
if ~iscell(FileName)
    if FileName ~= 0
        pat = importdata([PathName FileName], 'pat');
        
        if iscell(pat) && strcmp(FileName(end-4:end),'.vpat')
          handles.patfam = pat(2);
          handles.match = cell2mat(pat(3));
          handles.vName = sprintf('%s (v)',FileName(1:end-5));
          tmp = pat(2);
          save('newpattern.mat', 'tmp')
        
          pat = cell2mat(pat(1));
        end
        
        if ~isempty(pat)
            if strcmp(FileName(end-4:end),'.vpat')
            s   = sprintf('%s (v)',FileName(1:end-5));
            elseif strcmp(FileName(end-3:end),'.pat')
                s   = FileName(1:end-4);
            end
            pat_list_entries = get(handles.Pattern_List, 'String');
            p   = numel(pat_list_entries)+1 ;
            for n = 1:numel(pat_list_entries)
                if strcmp(pat_list_entries(n), s)
                    p = n;
                end
            end
                        
            tmp = size(handles.tcspc);
            if prod(1.*(size(pat) == tmp))
                pat_list_entries(p) = {s};
                
                handles.pat(p,:,:,:) = pat./sum(sum(sum(pat)));
                
                set(handles.Pattern_List,'Value',1);
                set(handles.Pattern_List,'String',pat_list_entries);
            elseif (tmp(1)==size(pat,1))&&(tmp(3)==size(pat,3))
                pat_list_entries(p) = {s};
                if tmp(2)<size(pat,2)
                    pat = pat(:,1:tmp(2),:);
                    handles.pat(p,:,:,:) = pat./sum(sum(sum(pat)));
                else
                    handles.pat(p,:,:,:) = zeros(tmp(1),tmp(2),tmp(3));
                    handles.pat(p,:,1:size(pat,2),:) = pat./sum(sum(sum(pat)));
                end
                set(handles.Pattern_List,'Value',1);
                set(handles.Pattern_List,'String',pat_list_entries);
            else
                msgbox('Selected pattern is not compatible with data set!','Could not use patten','warn')
            end
        else
            msgbox('Selected file does not contain pattern information!','No pattern found','warn')
        end
    end
else
    for k = 1:numel(FileName)

        fname = FileName{k};

        pat = importdata([PathName fname], 'pat');
                
        if ~isempty(pat)
            if strcmp(fname(end-3:end),'.pat')
                s   = fname(1:end-4);
            end
            pat_list_entries = get(handles.Pattern_List, 'String');
            p   = numel(pat_list_entries)+1 ;
            for n = 1:numel(pat_list_entries)
                if strcmp(pat_list_entries(n), s)
                    p = n;
                end
            end
                        
            tmp = size(handles.tcspc);
            if prod(1.*(size(pat) == tmp))
                pat_list_entries(p) = {s};
                
                handles.pat(p,:,:,:) = pat./sum(sum(sum(pat)));
                
                set(handles.Pattern_List,'Value',1);
                set(handles.Pattern_List,'String',pat_list_entries);
            elseif (tmp(1)==size(pat,1))&&(tmp(3)==size(pat,3))                   
                pat_list_entries(p) = {s};
                if tmp(2)<size(pat,2)
                    pat = pat(:,1:tmp(2),:);
                    handles.pat(p,:,:,:) = pat./sum(sum(sum(pat)));
                else
                    handles.pat(p,:,:,:) = zeros(tmp(1),tmp(2),tmp(3));
                    handles.pat(p,:,1:tmp(2),:) = pat./sum(sum(sum(pat)));
                end
                set(handles.Pattern_List,'Value',1);
                set(handles.Pattern_List,'String',pat_list_entries);
            else
                msgbox('Selected pattern is not compatible with data set!','Could not use patten','warn')
            end
        else
            msgbox('Selected file does not contain pattern information!','No pattern found','warn')
        end
        guidata(hObject, handles);        
    end
end

set(handles.LoadPattern_button,'Enable','on');

guidata(hObject, handles);


% --- Executes on button press in SavePattern_button.
function SavePattern_button_Callback(hObject, eventdata, handles)  % Save Pattern
% hObject    handle to SavePattern_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.SavePattern_button,'String','Busy');
set(handles.SavePattern_button,'Enable','off');
pause(0.05);

list_entries = get(handles.Pattern_List, 'String');
index_selected = get(handles.Pattern_List, 'Value');

for n = 1:numel(index_selected)

    m = index_selected(n);
    [FileName,PathName] = uiputfile({'*.pat','Pattern';...
        '*.*','All Files'},'Save pattern',...
        [handles.pathname 'newpattern.pat']);
    if FileName ~= 0
        pat = shiftdim(handles.pat(m,:,:,:),1);
        save([PathName FileName], 'pat');
    end
end

set(handles.SavePattern_button,'String','Save Pattern');
guidata(hObject, handles);
pause(0.05);


% --- Executes on button press in ClearPattern_button.
function ClearPattern_button_Callback(hObject, eventdata, handles)
% hObject    handle to ClearPattern_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.ClearPattern_button,'Enable','off');
pause(0.05);

list_entries = get(handles.Pattern_List, 'String');
index_selected = get(handles.Pattern_List, 'Value');

handles.pat(index_selected,:,:,:) = [];
list_entries(index_selected) = [];

set(handles.Pattern_List,'Value',1);
set(handles.Pattern_List,'String',list_entries);

guidata(hObject, handles);
pause(0.05);


% --- Executes on button press in ShowPattern_button.
function ShowPattern_button_Callback(hObject, eventdata, handles) % Select similar
% hObject    handle to ShowPattern_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.ShowPattern_button,'Enable','off');

list_entries = get(handles.Pattern_List, 'String');
index_selected = get(handles.Pattern_List, 'Value');

scrsz = get(0,'ScreenSize');

fwidth = 1200;
fheight = round(2*fwidth/3);
fx = 0.5*(scrsz(3)-fwidth);
fy = 0.5*(scrsz(4)-fheight);

py1 = 0.06;
py2 = 0.555;
wy  = 0.435;
px1 = 0.06;
px2 = 0.39;
wx  = 0.27;

nch   = size(handles.tcspc,1);
binw  = handles.head.tauw(:)./handles.head.tauw(1);
tau   = handles.head.tau(:)*ones(1,handles.par_num_PIE);
tmp   = ones(size(tau,1),1)*((0:handles.par_num_PIE-1).*tau(end,1));
tau   = reshape(tau+tmp, [numel(tau) 1]);

figure(handles.figPA);
set(gcf,'Position',[fx fy fwidth fheight]);
set(gcf,'Name','Selected Patterns','NumberTitle','off');

%    colormap(map);

pat = handles.pat(index_selected,:,:,:);

pat0 = zeros(size(pat,1),1,size(pat,3),size(pat,4));

pat = cat(2, pat0, pat, pat0);

decay    = sum(pat,2)./permute(repmat(binw, [1 size(pat,1) 1 size(pat,4)]),[2 3 1 4]);
decay    = reshape(decay, [size(decay,1) size(decay,3)*size(decay,4)]);
decay    = decay./max(max(decay));
%decay    = decay./max(decay,[],2);
decay(decay<=0) = 1e-4;

R1 = subplot('Position', [px1 py2 wx wy]);
hold on;
cla
box on;
for n = 1:numel(index_selected)
    plot(tau,log10(decay(n,:)),'-','LineWidth',2,'Color',handles.ROIcolor(n,:));
end
hold off;

set(R1,'FontSize',9,'color', [204 204 204]./255);

axis([0 floor(tau(end)) -4 0.1])

xlabel('time / ns','FontWeight','bold','FontSize',10,'FontAngle','italic');
ylabel('log (rel. frequency)','FontWeight','bold','FontSize',10,'FontAngle','italic');


R2 = subplot('Position', [px1 py1 wx wy]);

hold on;
cla
box on;

spec = sum(sum(pat,4),3);

lam = (handles.par_lam_start+handles.par_lam_step.*(-1:nch));
ts = lam(1): 1: lam(end);

for n = 1:numel(index_selected)
    %plot(ts,lam,'*','LineWidth',2,'Color',handles.ROIcolor(n,:));
    plot(ts,interp1(lam, spec(n,:), ts,'pchip'),'-','LineWidth',2,'Color',handles.ROIcolor(n,:));
end
hold off

set(R2,'FontSize',9,'color', [204 204 204]./255);

axis([max([lam(1) handles.par_clam_min]) min([lam(end) handles.par_clam_max]) 0 1])

xlabel('wavelength of spectral mean / nm','FontWeight','bold','FontSize',10,'FontAngle','italic');
ylabel('log (rel. frequency)','FontWeight','bold','FontSize',10,'FontAngle','italic');

R3 = subplot('Position', [px2 0.1 0.6 0.8]);

hold on;
cla
box on;

pat = pat./permute(repmat(binw,[1 size(pat,1) size(pat,2) size(pat,4)]),[2 3 1 4]);
pat = reshape(pat, [size(pat,1) size(pat,2) size(pat,3)*size(pat,4)]);
pat = pat./max(max(max(pat)));
pat(pat<=1e-3) = 1e-3;

for n = 1:numel(index_selected)
    
    if size(pat,2)>3
        [X,Y] = meshgrid(tau, lam(2:end-1)+n*2);
        sh = waterfall(X, Y, log10(shiftdim(pat(n,2:end-1,:),1)));
    else
        [X,Y] = meshgrid(tau, lam+n*2);
        sh = waterfall(X, Y, log10(shiftdim(pat(n,:,:),1)));
    end
    set(sh,'LineWidth',3,'EdgeColor',[0.2 0.2 0.2],'FaceColor',handles.ROIcolor(n,:));
    set(sh,'FaceAlpha',0.5,'EdgeAlpha',0.5);
end
set(R3,'CameraPosition',[200 -300 4]);
set(R3,'Clipping','on');
hold off

set(R3,'FontSize',9,'color', [204 204 204]./255);
axis([0 floor(tau(end)) max([lam(1) handles.par_clam_min]) min([lam(end) handles.par_clam_max]) -3 0])

xlabel('time / ns','FontWeight','bold','FontSize',10,'FontAngle','italic');
ylabel('wavelength of spectral mean / nm','FontWeight','bold','FontSize',10,'FontAngle','italic');
zlabel('rel. irradiance','FontWeight','bold','FontSize',10,'FontAngle','italic');


set(handles.ShowPattern_button,'String','Show pattern');
set(handles.Pattern_List,'Value',1);
pause(0.05);
guidata(hObject, handles);


% --- Executes on button press in FindAmplitudes_button.
function FindAmplitudes_button_Callback(hObject, eventdata, handles)
% hObject    handle to FindAmplitudes_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%set(handles.FindAmplitudes_button,'Enable','off');

list_entries   = get(handles.Pattern_List, 'String');
index_selected = get(handles.Pattern_List, 'Value');

if numel(index_selected) >= 1
    
    pat    = handles.pat(index_selected, :,:,:);
    list   = list_entries(index_selected);
    params = zeros(numel(index_selected), 1);
    
    %Variable Patterns Anfang 
    
    
     if isfield(handles,'vName')
        if isfield(handles,'patfam') && any(ismember(list,handles.vName))
            patfam = handles.patfam;
            handles.vParam = cell2mat(patfam{1,1}(end));
            tmp = patfam{1,1}(1:end-1);
            handles.fam = cat(4,tmp{1,1},tmp{1,2});
            for n = 1:size(tmp,2)-2
                handles.fam = cat(4,handles.fam,tmp{1,n+2});
            end
            handles.fam = permute(handles.fam,[4,1,2,3]);
            num_vpat = size(handles.fam,1);
        end
     end
    
 if isfield(handles,'fam') %&& any(isfield(handles.patfam, list))
        s = sprintf('Set_PatternOptions(\''file\'',\''Pattern.mat\'',\''%s\'', patfam, stim);',...
            [handles.pathname handles.filename]);
        match = 1;
    save('Pattern.mat', 'pat', 'list', 'params');    
    load([handles.pathname handles.timname], 'stim')
    
    binw   = handles.head.tauw./handles.head.tauw(1);
    DC     = handles.DC;
    aAP    = 1 + handles.AP;
    nx      = size(stim,1);
    ny      = size(stim,2);
    numpix  = nx*ny;
    nch     = size(stim,3);
    nbin    = size(stim,4);
    num_PIE = size(stim,5);
    num_pat = size(pat,1);
    
    stim = reshape(stim, numpix, nch, nbin, num_PIE);

    BG = sum(stim,3); 
    
     uiwait(eval(s));

    clear BG;

    load('Pattern.mat', 'params');
    
    for n = 1:numel(list)
        name = list{n};
        if match == 1 && strcmp(name(end-2:end),'(v)')
            k = n;
            handles.k = k;
        end
    end
    
    res = zeros(numpix, num_pat+1);

    mode = sum(params(1,1:2).*[1 2],2);

    ch_ind = params(:,4:end-1) == 1;

    ch_ind = permute(reshape(ch_ind,[num_PIE num_pat nch]),[2 3 1]);
    ch_ind = permute(repmat(ch_ind, [1 1 1 nbin]), [1 2 4 3]);
    ch_ind = sum(ch_ind)>0;

    if mode == 3                                         % use spectrum and decay
        tmp = stim.*repmat(ch_ind, [numpix 1]);
        tmp = reshape(tmp,[numpix nch*nbin*num_PIE]);
    elseif mode == 1                                      % use decay only
        tmp = sum(stim.*repmat(ch_ind, [numpix 1]),2);
        tmp = reshape(tmp,[numpix nbin*num_PIE]);
    else                                                  % use spectrum only
        tmp = sum(stim.*repmat(ch_ind, [numpix 1]),3);
        tmp = reshape(tmp,[numpix nch*num_PIE]);
    end
    
    if match == 1 && exist('k','var')

        ch_ind = permute(repmat(ch_ind, [num_pat 1 1 1 num_vpat]),[5 1 2 3 4]);
        vParamDist = zeros(numpix,1);

        res = zeros(numpix, num_vpat, num_pat+1);

        M = permute(repmat(pat,[1 1 1 1 num_vpat]),[5 1 2 3 4]);
        handles.fam = permute(repmat(handles.fam,[1 1 1 1 2]),[5 1 2 3 4]);

        for n = 1:num_vpat
            M(n,k,:,:,:) = handles.fam(1,n,:,:,:);
        end
        M = M.*ch_ind;
        
        if mode == 3                                         % use spectrum and decay
            M   = reshape(M,[num_vpat num_pat nch*nbin*num_PIE]);
        elseif mode == 1                                      % use decay only
            M   = sum(M,2);
            M   = reshape(M,[num_vpat num_pat nbin*num_PIE]);
        else                                                  % use spectrum only
            M   = sum(M,3);
            M   = reshape(M,[num_vpat num_pat nch*num_PIE]);
        end
        
        for n = 1:num_vpat
            for j=1:num_pat
                M(n,j,:) = M(n,j,:)/sum(M(n,j,:));  % normalization
            end
        end

        
        %if params(1,end)==1  %      orthogonal filtering

        %    weight = shiftdim(sum(tmp,1),1);

        %    D = diag(sum(weight)/num_pat./weight);

        %    master = (M*D*M')\M*D;

        %    for pix = 1:numpix
        %        res(pix,1:end-1)  = mean(repmat(tmp(pix,:),[num_pat 1]).*master,2);
        %        res(pix,end)      = sum(tmp(pix,:)) - sum(res(pix,:));
        %    end

        %else

        dist = zeros(numpix, num_vpat);
        final = zeros(numpix,num_pat+1);
        for pix = 1:numpix
            for n = 1:num_vpat

                res(pix,n,1:end-1)  = lsqnonneg(squeeze((M(n,:,:)))', tmp(pix,:)');
                %res(pix,n,1:end-1)  = squeeze((M(n,:,:)))'\tmp(pix,:)';
                res(pix,n,end)      = sum(tmp(pix,:)) - sum(res(pix,n,:));

                dist(pix,n) = res(pix,n,end);

            end
            [m,i] = min(abs(dist(pix,:)));
            vParamDist(pix) = handles.vParam(i);
            final(pix,:) = squeeze(res(pix,i,:));
        end
        vParamDist = reshape(vParamDist, [nx ny]);
        handles.vParamDist = vParamDist;
        res = final;
        
 else
    match = 0;
    
    %Variable Patterns Ende
    
    ind = true(size(params));
    for j=1:size(pat,1)
        s = sum(sum(sum(pat(j,:,:,:),4),3),2);
        if s > 0
            pat(j,:,:,:) = pat(j,:,:,:)/s;  % normalization
        else
            ind(j) = false;
        end
    end
    
    pat    = pat(ind,:,:,:);
    list   = list(ind);
    params = params(ind);
    
    %% select options for analysis
    
    save('Pattern.mat', 'pat', 'list', 'params');
    s = sprintf('Set_PatternOptions(\''file\'',\''Pattern.mat\'',\''%s\'');',[handles.pathname handles.filename]);
    
    uiwait(eval(s));
    
    load('Pattern.mat', 'params');
    
    n1 = size(pat,2);
    n2 = size(pat,3);
    
    handles.mode = sum(params(1,1:2).*[1 2],2);
    guidata(hObject, handles);
    
    handles.ch_ind = params(:,2+(1:n1)) == 1;
    handles.ch_ind = permute(repmat(handles.ch_ind,[1 1 n2]),[2 3 1]);
    handles.de_ind = params(:,3+n1:end) == 1;
    handles.de_ind = permute(repmat(handles.de_ind,[1 1 n1]),[3 2 1]);
    
    handles.ch_ind = (handles.ch_ind>0)&(handles.de_ind>0);
    
    e_ch  = squeeze(sum(sum(handles.ch_ind,3),2)>0);
    e_bin = squeeze(sum(sum(handles.ch_ind,3),1)>0);
    e_pie = squeeze(sum(sum(handles.ch_ind,1),2)>0);
    
    if prod([handles.mode sum(e_ch) sum(e_bin) sum(e_pie)]) > 0
                
%         handles.ch_ind = handles.ch_ind(:,e_ch,e_pie);
%         handles.ch_ind = permute(repmat(handles.ch_ind, [1 1 1 size(pat,3)]), [1 2 4 3]);
        
        load([handles.pathname handles.timname], 'stim')
        
        % Save data pre-processing
        
         % save('Pre-processed_FRET Slide.mat','stim','pat','-v7.3'); fg;
        
        stim = stim(:,:,e_ch,e_bin,e_pie);
        pat  = pat(:,e_ch,e_bin,e_pie);
        handles.ch_ind = handles.ch_ind(e_ch,e_bin,e_pie);
        
        binning = handles.par_binning;
        nx      = size(stim,1);
        ny      = size(stim,2);
        tx      = floor(nx/binning);
        ty      = floor(ny/binning);
        
        numpix  = tx*ty;
        nch     = size(stim,3);
        nbin    = size(stim,4);
        num_PIE = size(stim,5);
        num_pat = size(pat,1);
        
        if binning > 1
            tmp = stim(1:binning*tx,:,:,:,:);
            tmp = reshape(tmp, binning, tx, ny, nch, nbin, num_PIE);
            tmp = squeeze(sum(tmp,1));
            tmp = tmp(:,1:binning*ty,:,:);
            tmp = reshape(tmp, tx, binning, ty, nch, nbin, num_PIE);
            stim = permute(sum(tmp,2), [1 3 4 5 6 2]);
        end
        
        clear tmp;
                
        stim = reshape(stim, [numpix nch nbin num_PIE]); % rearrange the transformed signal for convenience
        
        M      = pat;
        M      = reshape(M,[num_pat nch*nbin*num_PIE]);
        handles.ch_ind = reshape(handles.ch_ind, [1 nch*nbin*num_PIE]);
        M      = M.*repmat(handles.ch_ind, [num_pat 1]);
        
        if handles.mode == 3                                         % use spectrum and decay
            
            tmp  = reshape(stim,[numpix nch*nbin*num_PIE]);
            tmp  = tmp.*repmat(handles.ch_ind, [numpix 1]);
            
            tmp  = reshape(tmp, [numpix nch nbin num_PIE]);
            tmp  = permute(tmp, [1 2 4 3]);
            tmp  = reshape(tmp,[numpix nch*num_PIE*nbin]);

            M  = reshape(M, [num_pat nch nbin num_PIE]);
            M  = permute(M, [1 2 4 3]);
            M  = reshape(M, [num_pat nch*num_PIE*nbin]);
            
            ebin = nbin;
            M(M<0) = 0;
            TT = tmp;
            TT(TT<0) = 0;
            T = TT;
            
        elseif handles.mode == 1                                      % use decay only
            
            handles.ch_ind = reshape(handles.ch_ind, [1 nch nbin num_PIE]);
            
            tmp = sum(stim.*repmat(handles.ch_ind, [numpix 1]),2);
            tmp = reshape(tmp,[numpix nbin*num_PIE]);
            
            M    = reshape(M,[num_pat nch nbin num_PIE]);
            M    = sum(M,2);
            M    = reshape(M,[num_pat nbin*num_PIE]);
            ebin = nbin;
            %%
            M(M<0) = 0;
            TT = tmp;
            TT(TT<0) = 0;
            T = TT;
            
        else                                                 % use spectrum only
            % handles.mode  == 2
            handles.ch_ind = reshape(handles.ch_ind, [1 nch nbin num_PIE]);
            
            tmp = sum(stim.*repmat(handles.ch_ind, [numpix 1]),3);
            tmp = reshape(tmp,[numpix nch*num_PIE]);
            
            M    = reshape(M,[num_pat nch nbin num_PIE]);
            M    = sum(M,3);
            M    = reshape(M,[num_pat nch*num_PIE]);
            ebin = 1;
            %%
            M(M<0) = 0;
            TT = tmp;
            TT(TT<0) = 0;
            T = TT;
        end
        
        clear stim TT;
        
%         reduce effective number of spectral channels only when all
%         information is used
        if handles.mode == 3
            
            kkk = [];
            
            if handles.mode > 1
                
                if handles.mode == 2
                    M   = reshape(M,[num_pat nch ebin num_PIE]);
                    MMM = sum(M,3);
                    MMM = reshape(MMM,[num_pat nch*num_PIE]);
                    M  = permute(M, [1 2 4 3]);
                    M  = reshape(M,[num_pat nch*num_PIE*ebin]);
                else
                    MMM = M;
                end
                
                
                for pulse = 1:num_PIE
                    off = (pulse-1)*nch;
                    ind = (1:nch)+off;
                    
                    [dum,b] = sort(sum(MMM(:,ind),2),'descend');
                    
                    mm  = min([numel(b) 5]);
                    A   = MMM(b(1:mm),ind);
                    
                    [dum, k] = sort(A,1,'descend');
                    
                    kk = [0 nch];
                    for nl = 1:mm
                        kk = [kk find(diff(k(nl,:))~=0)];
                    end
                    kk = sort(unique(kk));
                    
                    while numel(kk) < mm+1
                        [dk, ok] = max(diff(kk));
                        kk = sort([kk kk(ok)+round(dk/2)]);
                    end
                    kk  = off + sort([kk kk(1:end-1)+round(diff(kk)./2)]);
                    kkk = [kkk kk];
                end
                
                kk = unique(kkk);
                nk = numel(kk)-1;
                
                if nk>0
                    MM = zeros(num_pat,nk*ebin);
                    TT = zeros(numpix,nk*ebin);
                    off = num_PIE*nch.*((1:ebin)-1);
                    for o = 1:numel(off)
                        for l = 1:nk
                            MM(:,(o-1)*nk+l) = sum(  M(:,off(o)+((kk(l)+1):kk(l+1))),2);
                            TT(:,(o-1)*nk+l) = sum(tmp(:,off(o)+((kk(l)+1):kk(l+1))),2);
                        end
                    end
                else
                    MM = M;
                    TT = tmp;
                end
                
                MM = reshape(MM,[num_pat nk ebin]);
                MM = permute(MM,[1 3 2]);
                MM = reshape(MM,[num_pat ebin*nk]);
                
                TT = reshape(TT,[numpix nk ebin]);
                TT = permute(TT,[1 3 2]);
                TT = reshape(TT,[numpix ebin*nk]);
            else
                MM = M;
                TT = tmp;
                nk = num_PIE;
            end
            
            clear tmp M MMM;
            
            if handles.mode ~=2
                
                kkk = [];
                
                for chan = 1:nk
                    off = (chan-1)*ebin;
                    ind = (1:ebin)+off;
                    
                    [tmp, b] = sort(sum(MM(:,ind(32:end)),2),'descend');
                    
                    mm  = min([numel(b) 3]);
                    A   = MM(b(1:mm),ind(32:end));
                    
                    [tmp, k] = sort(A,1,'descend');
                    
                    kk = [31 nbin];
                    for nl = 1:mm
                        kk = [kk 32+find(diff(k(nl,:))~=0)];
                    end
                    kk = sort(unique(kk));
                    
                    while numel(kk) < 4
                        [dk, ok] = max(diff(kk));
                        kk = sort([kk kk(ok)+round(dk/2)]);
                    end
                    kk  = off + sort([0 kk kk(1:end-1)+round(diff(kk)./3) kk(1:end-1)+round(2*diff(kk)./3)]);
                    kkk = [kkk kk];
                end
                
                kk = unique(kkk);
                nk = numel(kk)-1;
                
                if nk>0
                    M = zeros(num_pat,nk);
                    T  = zeros(numpix,nk);
                    for l = 1:nk
                        M(:,l) = sum(MM(:,(kk(l)+1):kk(l+1)),2);
                        T(:,l) = sum(TT(:,(kk(l)+1):kk(l+1)),2);
                    end
                else
                    M = MM;
                    T = TT;
                end
                
                T(T<0)= 0;
                
            else
                M = MM;
                T = TT;
            end
        end
        
        nk = size(T,2); % number of bins in pattern data
        
        % do fourier transform

        ek1  = ceil(tx/2);
        ek2  = ceil(ty/2);

        EF = EFilter([tx+2*ek1 ty+2*ek2], 0.80, 0.80, 9); % optional low-pass filter
                
        T(T<0)= 0;
        
        T = reshape(T,[tx ty nk]); % tx = ty = numPixX/numPixY
        
        ftmp  = zeros(size(EF,1), size(EF,2), nk);
        
        for k = 1:nk
            tmp   = zeros(size(EF,1), size(EF,2));
            tmp(ek1+(1:tx),ek2+(1:ty)) = T(:,:,k);
            tmp = EF.*(fftshift(fftshift(fft2(tmp),1),2));   % norm of 2D-fourier-transform of singal components
            ftmp(:,:,k) = abs(ifft2(ifftshift(ifftshift(tmp,1),2)));
        end
        
        T  = ftmp(ek1+(1:tx),ek2+(1:ty),:);
        T(T<0)= 0;
        
        otmp = sum(T,3);
        
        T = reshape(T,[numpix nk]);
        MMM = M';
           
        res = KLF(MMM,T); % Linear unmixing step
        
        res(res<0)= 0;
        %res = reshape(res,[tx ty num_pat]);
        %es haben "rest" Reihen am Ende vom res Array gefehlt, ergänzt werden sie jetzt mit Nullen, KRM 11-2022
        %disp(tx);
        %disp(ty);
        %disp(num_pat);
        sizeres = size(res);
        sizeres1 = sizeres(1);
        rest = (tx * ty) - sizeres1;
        ERG = zeros (rest,num_pat);
        res = cat(1,res,ERG); 
        
        res = reshape(res,tx,ty,num_pat);
        
        res(:,:,num_pat+1) =  otmp - sum(res,3);
        
        set(handles.SaveResults_button,'Enable','on');
        handles.results.pat    = pat;
        handles.results.params = params;
        handles.results.amp    = res;
    end
  end    
end

handles.res = res;
guidata(hObject, handles);



% --- Executes on button press in MatchResults_button.
function MatchResults_button_Callback(hObject, eventdata, handles)
% hObject    handle to MatchResults_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Plot_Results(handles)



% --- Executes on button press in FitPattern_button.
function FitPattern_button_Callback(hObject, eventdata, handles)
% hObject    handle to FitPattern_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%#function Pattern_Fit

set(hObject,'Enable','off');

list_entries   = get(handles.Pattern_List, 'String');
index_selected = get(handles.Pattern_List, 'Value');

if numel(index_selected) >= 1

    pat    = handles.pat(index_selected, :,:,:);
    list   = list_entries(index_selected);
    params = zeros(numel(index_selected), 1);
    binw   = handles.head.tauw./handles.head.tauw(1);
    DC     = handles.DC;
    aAP    = 1 + handles.AP;
        
    save('Pattern.mat', 'pat', 'list', 'params');

    handles.filename;
    
    s = sprintf('Pattern_Fit(\''file\'',\''Pattern.mat\'',\''%s\'');',...
        [handles.pathname handles.filename]);
    
    evalin('base', s);
    
end


% --- Executes on button press in SaveResults_button.
function SaveResults_button_Callback(hObject, eventdata, handles)
% hObject    handle to SaveResults_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FileName,PathName] = uiputfile({'*.mat','MatLab-data';...
    '*.*','All Files' },'Save results',...
    [handles.pathname handles.filename(1:end-4) '_results.mat']);

results = handles.results;
results.x = handles.head.ImgHdr.X0+(1:handles.par_binning:handles.head.ImgHdr.PixX)*handles.head.ImgHdr.PixelSize;
results.y = handles.head.ImgHdr.Y0+(1:handles.par_binning:handles.head.ImgHdr.PixY)*handles.head.ImgHdr.PixelSize;
results.image = handles.tag;
results.tau = handles.head.tau;

pat = handles.results.pat;
tcspc = handles.tcspc;
irf   = handles.IRF;

tmp = round(handles.head.tauw./handles.head.tauw(1));
tau = handles.head.tau;
tau = tau(1) + (tau(2)-tau(1)).*(0:sum(tmp)-1);
ind = [0 cumsum(tmp)];
for k = 1:numel(tmp)
    tst_p(:,:,ind(k)+1:ind(k+1),:) = repmat(pat(:,:,k,:), [1 1 tmp(k) 1 1])./tmp(k);    
    tst_t(:,ind(k)+1:ind(k+1),:) = repmat(tcspc(:,k,:), [1 tmp(k) 1 1]);% ./tmp(k);    
    tst_i(:,ind(k)+1:ind(k+1),:) = repmat(irf(:,k,:), [1 tmp(k) 1 1]); % ./tmp(k);    
end

results.tau   = tau;
results.pat   = tst_p;
results.tcspc = tst_t;
results.IRF   = tst_i;

save([PathName FileName],'results');



% --- Executes on button press in PlotResults_button.
function PlotResults_button_Callback(hObject, eventdata, handles)
% hObject    handle to PlotResults_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

load Get_Params;
scrsz = get(0,'ScreenSize');

fwidth = 1200;
fheight = round(2*fwidth/3);
fx = 0.5*(scrsz(3)-fwidth);
fy = 0.5*(scrsz(4)-fheight);

py1 = 0.06;
py2 = 0.555;
wy  = 0.435;
px1 = 0.06;
px2 = 0.39;
px3 = 0.72;
wx  = 0.27;

xl = [399 440 485 540 580 610 630 699];
yl = [[0 0 0]; [0 0 1]; [0 1 1]; [0 1 0]; [1 1 0]; [1 0.65 0]; [1 0 0]; [1 0 1]];
lambda   = 399:3:699;
spectrum = interp1(xl, yl, lambda);

d_lambda = (handles.par_lam_start+handles.par_lam_step.*(0:7));
dl       = abs(ones(size(lambda'))*d_lambda - lambda'*ones(size(d_lambda)));
k_ind    = mod(find(dl==ones(size(lambda'))*min(dl)),numel(lambda));

ntau   = handles.head.tau;

set(handles.PlotResults_button,'String','Busy');
% set(handles.PlotResults_button,'Enable','off');
pause(0.05);

nx      = handles.head.ImgHdr.PixX;
ny      = handles.head.ImgHdr.PixY;
tx      = floor(nx/binning);
ty      = floor(ny/binning);
nch     = size(handles.tcspc,1);
nbin    = size(handles.tcspc,2);


load([handles.pathname 'Moments.mat'],'rtau');
load([handles.pathname handles.timname],'stim');
load('c_map.mat','map');

stim = shiftdim(stim(:,:,:,:,handles.pulse),1); %#ok<COLND>

if binning > 1
    tmp = stim(1:binning*tx,:,:,:);
    tmp = reshape(tmp, binning, tx, ny, nch, nbin);
    tmp = squeeze(mean(tmp,1));
    tmp = tmp(:,1:binning*ty,:);
    tmp = reshape(tmp, tx, binning, ty, nch, nbin);
    stim = permute(mean(tmp,2), [1 3 4 5 6 2]);
end

stim = reshape(stim, size(stim,1)*size(stim,2),size(stim,3),size(stim,4));
rtau = reshape(rtau, size(rtau,1)*size(rtau,2),size(rtau,3),size(rtau,4),size(rtau,5));

list_entries = get(handles.listbox2, 'String');

for n = 1:numel(list_entries)

    s   = list_entries{n};
    ROI = str2double(s(end-2));

    name = sprintf('Results for ROI %d Excitation pulse %d', ROI, handles.pulse);

    figure(figone+n);
    set(gcf,'Position',[fx fy fwidth fheight]);
    set(gcf,'Name',name,'NumberTitle','off');

    colormap(map);

    ind = bitget(handles.s_ind(:, handles.pulse), ROI) == 1;

    %   num_pix = sum(ind);
    data = reshape(handles.tag, size(handles.tag,1)*size(handles.tag,2),size(handles.tag,3),size(handles.tag,4));
    tavg  = reshape(handles.tav, size(handles.tav,1)*size(handles.tav,2),size(handles.tav,3),size(handles.tav,4),size(handles.tav,5));

    data = data(ind,:,handles.pulse);
    tavg = squeeze(tavg(:,:,handles.pulse,:));
    tavg = tavg(ind,:,:);
    atau = rtau(ind,:,:,handles.pulse);

    data(data<handles.par_Threshold) = 0;
    atau(atau<0) = 0;

    V      = squeeze(sum(stim(ind,:,:),1))./(ones(size(stim,2),1)*handles.head.tauw);
    V(V<0) = 0.1;
    V      = V./(max(max(V)));

    R1 = subplot('Position', [px1 py2 wx wy]);
    hold on;
    cla
    box on;

    for ch = 1:size(V,1)
        plot(ntau,log10(V(ch,:)),'-','Color',spectrum(k_ind(ch),:),'LineWidth',2);
    end
    hold off;

    set(R1,'FontSize',9,'color', [204 204 204]./255);

    axis([0 floor(ntau(end)) -3 0.1])

    xlabel('time / ns','FontWeight','bold','FontSize',10,'FontAngle','italic');
    ylabel('log (rel. frequency)','FontWeight','bold','FontSize',10,'FontAngle','italic');

    R2 = subplot('Position', [px2 py2 wx wy]);

    hold on;
    cla
    box on;

    K = 1:size(data,2);
    K = repmat(reshape(K,[1 numel(K)]),[size(data,1) 1]);

    tZ = data.*K;
    tZ = handles.par_lam_start + handles.par_lam_step.*(sum(tZ,2)./sum(data,2)-1.0);

    lind = handles.par_lam_start+(0:2:handles.par_lam_step*(size(data,2)-1));

    tmp2 = sum(data,2);
    [tmp1, ord] = sort(tZ);
    tmp2 = tmp2(ord);

    lam = handles.par_lam_start+(-1:(1+3*size(data,2))).*handles.par_lam_step./3;
    s = zeros(size(lam));
    for k = 2:numel(lam)-1
        ind  = (tmp1>lam(k-1))&(tmp1<lam(k+1));
        s(k) = sum(tmp2(ind));
    end

    tmp = max(s);
    if tmp>0
        s = s./max(s);
    else
        s(:) = 0;
    end;

    ts = 400: 0.5: 700;
    cs = interp1(lam, s, ts,'cubic');

    plot(ts,cs,'-','Color',[1 0 0]);
    hold off

    set(R2,'FontSize',9,'color', [204 204 204]./255);

    % axis([10*floor(lind(1)/10) 10*ceil(lind(end)/10) 0 1.1])
    axis([handles.par_clam_min handles.par_clam_max 0 1.1])

    xlabel('wavelength of spectral mean / nm','FontWeight','bold','FontSize',10,'FontAngle','italic');
    ylabel('rel. spectral irradiance','FontWeight','bold','FontSize',10,'FontAngle','italic');


    R3 = subplot('Position', [px3 py2 wx wy]);

    hold on;
    cla
    box on;
    for ch = 1:size(tavg,2)
        la = d_lambda(ch);
        tmp = sort(tavg(:,ch,1));
        bw(1) = tmp(round(0.5*numel(tmp)));
        bw(2) = tmp(round(0.25*numel(tmp)));
        bw(3) = tmp(round(0.75*numel(tmp)));
        bw(4) = max(tmp(tmp<(bw(1)+1.5*(bw(3)-bw(2)))));
        bw(5) = min(tmp(tmp>(bw(1)-1.5*(bw(3)-bw(2)))));


        line([la la],[bw(4) bw(5)],'LineWidth',5);
        line([la-3 la+3],[bw(1) bw(1)],'Color',[0 0 1],'LineWidth',2);
        line([la-3 la+3 la+3 la-3 la-3],[bw(2) bw(2) bw(3) bw(3) bw(2)],'Color',[0 0 1],'LineWidth',1);
    end;
    hold off

    set(R3,'FontSize',9,'color', [204 204 204]./255);

    axis([10*floor(lam(1)/10) 10*ceil(lam(end)/10) handles.par_tau_min handles.par_tau_max])

    xlabel('detector wavelength / nm','FontWeight','bold','FontSize',10,'FontAngle','italic');
    ylabel('mean fluorecence lifetime / ns','FontWeight','bold','FontSize',10,'FontAngle','italic');

    R4 = subplot('Position', [px1 py1 wx wy]);

    cla

    M_X = data.*squeeze(tavg(:,:,1));
    M_X = squeeze(sum(M_X,2)./sum(data,2));

    M_Y = data.*squeeze(tavg(:,:,2));
    M_Y = squeeze(sum(M_Y,2)./sum(data,2));

    cscatter(gca, M_X, M_Y);

    set(R4,'PlotBoxAspectRatio',[1 1 1], ...
        'Box','on', ...
        'XDir','normal', ...
        'YDir','normal', ...
        'FontSize',9,...
        'Color',[50 50 50]./255);

    axis([handles.par_tau_min handles.par_tau_max handles.par_cum2_min handles.par_cum2_max]);
    xlabel('mean fluorecence lifetime / ns','FontWeight','bold','FontSize',10,'FontAngle','italic');
    ylabel('sqrt(2nd cumulant) / ns','FontWeight','bold','FontSize',10,'FontAngle','italic');

    R5 = subplot('Position', [px2 py1 wx wy]);

    hold on
    cla
    box on

    M_X = data.*squeeze(atau(:,1,:));
    M_X = squeeze(sum(M_X,2)./sum(data,2));

    M_Y = data.*squeeze(atau(:,2,:));
    M_Y = squeeze(sum(M_Y,2));

    amax = max(M_Y);

    scatter(M_X, M_Y, 10, [1 0 0],'diamond','Clipping','on');

    M_X = data.*squeeze(atau(:,3,:));
    M_X = squeeze(sum(M_X,2)./sum(data,2));

    M_Y = data.*squeeze(atau(:,4,:));
    M_Y = squeeze(sum(M_Y,2));

    amax = max([1 amax max(M_Y)]);

    scatter(M_X, M_Y, 10, [0 1 0],'diamond','Clipping','on');

    set(R5,'PlotBoxAspectRatio',[1 1 1], ...
        'Box','on', ...
        'XDir','normal', ...
        'YDir','normal', ...
        'FontSize',9,...
        'Color',[50 50 50]./255);

    axis([handles.par_tau_min 8*handles.par_tau_max 0  ceil(amax)]);
    xlabel('fluorecence lifetime / ns','FontWeight','bold','FontSize',10,'FontAngle','italic');
    ylabel('Amplitude','FontWeight','bold','FontSize',10,'FontAngle','italic');
    hold off


    R6 = subplot('Position', [px3 py1 wx wy]);
    cla

    M_X = data.*squeeze(atau(:,1,:));
    M_X = squeeze(sum(M_X,2)./sum(data,2));

    M_Y = data.*squeeze(atau(:,3,:));
    M_Y = squeeze(sum(M_Y,2)./sum(data,2));

    cscatter(gca, M_X, M_Y);

    set(R6,'PlotBoxAspectRatio',[1 1 1], ...
        'Box','on', ...
        'XDir','normal', ...
        'YDir','normal', ...
        'FontSize',9,...
        'Color',[50 50 50]./255);

    axis([handles.par_tau_min 8*handles.par_tau_max handles.par_tau_min 8*handles.par_tau_max]);

    xlabel('fluorescence lifetime 1 / ns','FontWeight','bold','FontSize',10,'FontAngle','italic');
    ylabel('fluorescence lifetime 2 / ns','FontWeight','bold','FontSize',10,'FontAngle','italic');
end

set(handles.PlotResults_button,'String','Plot ROI results');
set(handles.PlotResults_button,'Enable','on');
pause(0.05);


function [M_X, M_Y, Z, ZZ] = compute_data(handles)

tag = handles.tag(:,:,:,handles.pulse);
tav = handles.tav(:,:,:,handles.pulse,1:2);

for n = 1:size(tag,3)
    tag(:,:,n) = tag(:,:,n).*handles.channel(handles.pulse,n);
end

W = tag;
W(W<handles.par_Threshold) = 0;

N = sum(W,3);

K = ones(size(tag));
for n = 1:size(K,3)
    K(:,:,n,:) = n;
end;

% weighted avarage of first cumulant

M_X = W.*tav(:,:,:,1,1);
M_X = squeeze(sum(M_X,3))./N;

% weigthed average of second cumulant

M_Y = W.*tav(:,:,:,1,2);
M_Y = squeeze(sum(M_Y,3))./N;

% weigthed average of spectral channel

Z = W.*K;
Z = squeeze(sum(Z,3))./N - 1.0;

% spectral width

ZZ = W.*abs(K-repmat(Z,[1 1 size(tag,3)]));
ZZ = squeeze(sum(ZZ,3))./N; % +1;

Z = (handles.par_lam_start+handles.par_lam_step.*Z);
ZZ  = (ZZ.*handles.par_lam_step);

ind = N==0;
M_X(ind) = 0;
M_Y(ind) = 0;
Z(ind)   = 0;
ZZ(ind)  = 0;

num = size(W);
num = prod(num(1:2));

M_X = reshape(M_X,num,1);
M_Y = reshape(M_Y,num,1);
Z   = reshape(  Z,num,1);
ZZ  = reshape( ZZ,num,1);



function compute_fFLIM(handles)

t_m = (handles.tav(:,:,:,:,1));

% ind = sum(sum(handles.tag,4),3) < handles.par_Threshold;
ind = handles.tag >= handles.par_Threshold;
t_m = ind.*t_m;

save([handles.pathname 'fFLIM.mat'], 't_m');


%% plot_graphs

function plot_graphs(handles)

binning = handles.par_binning;

x = handles.head.ImgHdr.X0+(1:binning:handles.head.ImgHdr.PixX)*handles.head.ImgHdr.PixelSize;
y = handles.head.ImgHdr.Y0+(1:binning:handles.head.ImgHdr.PixY)*handles.head.ImgHdr.PixelSize;

ind_x = x >= handles.plt_x_min & x <= handles.plt_x_max;
ind_y = y >= handles.plt_y_min & y <= handles.plt_y_max;

xl = [399 440 485 540 580 610 630 699];
yl = [[0 0 0]; [0 0 1]; [0 1 1]; [0 1 0]; [1 1 0]; [1 0.65 0]; [1 0 0]; [1 0 1]];
lambda   = 399:3:699;
spectrum = interp1(xl, yl, lambda);

lims = [handles.plt_tau_min handles.plt_tau_max];

val = squeeze(handles.tav(:,:,:,handles.pulse,1));

data = squeeze(handles.tag(:,:,:,handles.pulse));
data = data.*(data >= handles.par_Threshold);

for n = 1:size(data,3)
    data(:,:,n) = data(:,:,n).*handles.channel(n);
    val(:,:,n)  = val(:,:,n).*handles.channel(n);
end

val = sum(val.*data,3)./sum(data,3);

data   = squeeze(sqrt(sum(data,3)));
tmp    = max(max(data));
data   = data./tmp;

data(data<0)=0;
data = repmat(data, [1 1 3]);

val = reshape(val,size(val,1)*size(val,2),1);

val(isnan(val))  = 0;
val(val<lims(1)) = 0;
val(val>lims(2)) = lims(2);

k      = 2 + round((val-lims(1))./(lims(2)-lims(1)).*(size(spectrum,1)-2));
k(k<2) = 1;

im = spectrum(k,:);
im = reshape(im,size(data,1),size(data,2),3);
im = im.*data;

axes(handles.axes1)

hold(handles.axes1,'off');
imagesc(x(ind_x), y(ind_y), im(ind_y,ind_x,:)); % data(ind_y,ind_x,:));
hold(handles.axes1,'on');

% colormap(spectrum)
alpha = zeros(size(handles.tag,1), size(handles.tag,2));
ovl   = zeros(numel(alpha),3);
for n = 1:handles.maxROI
    if handles.ROIlist(n,handles.pulse) == 1
        tst = bitget(handles.s_ind(:,handles.pulse), n)==1;
        alpha(tst==1) = 1;
        for k = 1:3
            tc  = repmat(handles.ROIcolor(n,k),[size(data,1) size(data,2) 1]);
            ovl(tst,k) = tc(tst);
        end;
    end;
end;
ovl = data.*reshape(ovl, size(data,1), size(data,2), 3);

iim2 = image(x(ind_x),y(ind_y), ovl(ind_y,ind_x,:));
set(iim2,'AlphaData',alpha(ind_y,ind_x));

set(handles.axes1, ...
    'Box','on', ...
    'XDir','normal', ...
    'YDir','reverse', ...
    'FontSize',9, ...
    'PlotBoxAspectRatio',[1 1 1]);
xlabel(handles.axes1, 'x / µm','FontWeight','bold','FontSize',10,'FontAngle','italic');
ylabel(handles.axes1, 'y / µm','FontWeight','bold','FontSize',10,'FontAngle','italic');

% axes(handles.axes3)

ind = handles.M_X > handles.plt_tau_min & handles.M_X < handles.plt_tau_max & ...
      handles.M_Y > handles.plt_cum2_min & handles.M_Y< handles.plt_cum2_max & ...
      handles.Z   > handles.plt_clam_min & handles.Z  < handles.plt_clam_max &...
      handles.ZZ  > handles.plt_wlam_min & handles.ZZ < handles.plt_wlam_max;


set(handles.axes3,'PlotBoxAspectRatio',[1 1 1], ...
    'Box','on', ...
    'XDir','normal', ...
    'YDir','normal', ...
    'FontSize',9, ...
    'Color',[50 50 50]./255);

cscatter(handles.axes3, handles.M_X(ind), handles.M_Y(ind));
hold(handles.axes3, 'on');

axis(handles.axes3, [handles.plt_tau_min handles.plt_tau_max handles.plt_cum2_min handles.plt_cum2_max]);

set(handles.axes3,'PlotBoxAspectRatio',[1 1 1], ...
    'Box','on', ...
    'XDir','normal', ...
    'YDir','normal', ...
    'FontSize',9, ...
    'Color',[50 50 50]./255);
xlabel(handles.axes3, '< \tau > / ns','FontWeight','bold','FontSize',10,'FontAngle','italic');
ylabel(handles.axes3, 'sqrt(2nd cumulant) / ns','FontWeight','bold','FontSize',10,'FontAngle','italic');

% for n = 1:handles.maxROI
%     if handles.ROIlist(n,handles.pulse) == 1
%         ind_1 = ind & bitget(handles.s_ind(:,handles.pulse), n)==1;
%     end
% end;

% axes(handles.axes2)

cscatter(handles.axes2, handles.ZZ(ind), handles.Z(ind));
hold(handles.axes2, 'on');

set(handles.axes2, 'PlotBoxAspectRatio',[1 1 1], ...
    'Box','on', ...
    'XDir','normal', ...
    'YDir','normal', ...
    'FontSize',9,...
    'Color',[50 50 50]./255);

xlabel(handles.axes2,'\Delta \lambda (nm)','FontWeight','bold','FontSize',10,'FontAngle','italic');
ylabel(handles.axes2,'\lambda_0 (nm)','FontWeight','bold','FontSize',10,'FontAngle','italic');
axis(handles.axes2, [handles.plt_wlam_min handles.plt_wlam_max handles.plt_clam_min handles.plt_clam_max]);

% for n = 1:handles.maxROI
%     if handles.ROIlist(n,handles.pulse) == 1
%         ind_1 = ind & bitget(handles.s_ind(:,handles.pulse), n)==1;
%     end
%     %variables for onegraphs.m
%     %         N(n,:) = histc((lam_start+lam_step.*Z(bitget(s_ind(:,pulse), n)==1,pulse)),edges);
%     %         N(n,:) = N(n,:) / sum(N(n,:));
% end

% axes(handles.axes4)

cscatter(handles.axes4,handles.M_X(ind),handles.Z(ind));
hold(handles.axes4, 'on');
set(handles.axes4,'PlotBoxAspectRatio',[1 1 1], ...
    'Box','on', ...
    'XDir','normal', ...
    'YDir','normal', ...
    'FontSize',9,...
    'Color',[50 50 50]./255);

xlabel(handles.axes4, '< \tau > / ns','FontWeight','bold','FontSize',10,'FontAngle','italic');
ylabel(handles.axes4, '\lambda_0 / nm','FontWeight','bold','FontSize',10,'FontAngle','italic');
axis(handles.axes4, [handles.plt_tau_min handles.plt_tau_max handles.plt_clam_min handles.plt_clam_max]);

for n = 1:handles.maxROI
    if handles.ROIlist(n,handles.pulse) == 1
        ind_1 = ind & bitget(handles.s_ind(:,handles.pulse), n)==1;
        scatter(handles.axes4,handles.M_X(ind_1), handles.Z(ind_1),10,handles.ROIcolor(n,:),'diamond','Clipping','on');
        scatter(handles.axes2,handles.ZZ(ind_1), handles.Z(ind_1),10,handles.ROIcolor(n,:),'diamond','Clipping','on');
        scatter(handles.axes3,handles.M_X(ind_1),handles.M_Y(ind_1),10,handles.ROIcolor(n,:),'diamond','Clipping','on');
    end
    %variables for onegraphs.m
    %         N(n,:) = histc((lam_start+lam_step.*Z(bitget(s_ind(:,pulse), n)==1,pulse)),edges);
    %         N(n,:) = N(n,:) / sum(N(n,:));
end

hold(handles.axes1, 'off');
hold(handles.axes2, 'off');
hold(handles.axes3, 'off');
hold(handles.axes4, 'off');


%% select_data

function [s_ind] = select_data(handles, ROInum)

s_ind = handles.s_ind;

axes(handles.axes1);

[xi,yi] = ginput(1);
pt1     = get(handles.axes1, 'CurrentPoint');
pt2     = get(handles.axes2, 'CurrentPoint');
pt3     = get(handles.axes3, 'CurrentPoint');
pt4     = get(handles.axes4, 'CurrentPoint');


t1 = xi-[pt1(1,1) pt2(1,1) pt3(1,1) pt4(1,1)];
t2 = yi-[pt1(1,2) pt2(1,2) pt3(1,2) pt4(1,2)];
p  = find((t1+t2)==0);

if isempty(p)
    p = abs(t1+t2);
    p = find(p==min(p));
end;

if p==1
    axes(handles.axes1);
elseif p==2
    axes(handles.axes2);
elseif p==3
    axes(handles.axes3);
else
    axes(handles.axes4);
end

hold on

plot(xi,yi,'+','Color',handles.ROIcolor(ROInum,:))
xy   = [xi;yi];
n    = 1;
but  = 1;
while but == 1
    [xi,yi,but] = ginput(1);
    plot(xi,yi,'+','Color',handles.ROIcolor(ROInum,:))
    n = n+1;
    xy(:,n) = [xi;yi];
end
n = n+1;
xy(:,n) = xy(:,1);

s1 = [xy(1,n-1)-xy(1,2); xy(2,n-1)-xy(2,2)];

xy = [-s1 xy -s1];

cs = spline(1:n, xy);
% Interpolate with a spline curve and finer spacing.
ts = 1: 0.1: n;
xys = ppval(cs,ts);

% Plot the interpolated curve.
plot(xys(1,:),xys(2,:),'Color',handles.ROIcolor(ROInum,:));
hold off

% % w: to save ROI
% wp(ROInum) = p;
% wn(ROInum) = n;
% eval(['wxy_' num2str(ROInum) '= xy;']);

if p==1 % Intensity image

    t_ind = reshape(s_ind(:,handles.pulse),size(handles.tag,1),size(handles.tag,2));

    mx = round((xys(1,:)-handles.head.ImgHdr.X0)./handles.head.ImgHdr.PixelSize/handles.par_binning);
    my = round((xys(2,:)-handles.head.ImgHdr.Y0)./handles.head.ImgHdr.PixelSize/handles.par_binning);
    mx(mx<1) = 1;
    my(my<1) = 1;
    mx(mx>size(t_ind,2)) = size(t_ind,2);
    my(my>size(t_ind,1)) = size(t_ind,1);

    dmx = abs(diff([0 mx]));
    dmy = abs(diff([0 my]));
    tmx = mx(1);
    tmy = my(1);
    for n = 2:numel(dmx)-1
        if (dmx(n)+dmy(n))~=0
            if (dmx(n)>1)||(dmy(n)>1)
                if dmx(n)>dmy(n)
                    nx = mx(n-1):sign(mx(n)-mx(n-1)):mx(n);
                    nx = nx(2:end);
                    ny = round(interp1([mx(n-1) mx(n)],[my(n-1) my(n)], nx));
                    tmx = [tmx nx];
                    tmy = [tmy ny];
                else
                    ny = my(n-1):sign(my(n)-my(n-1)):my(n);
                    ny = ny(2:end);
                    nx = round(interp1([my(n-1) my(n)],[mx(n-1) mx(n)], ny));
                    tmx = [tmx nx];
                    tmy = [tmy ny];
                end
            else
                tmx = [tmx mx(n)];
                tmy = [tmy my(n)];
            end;
        end;
    end;

    for n = 1:numel(tmx)
        t_ind(tmy(n),tmx(n)) = bitset(t_ind(tmy(n),tmx(n)), ROInum);
    end;
    y_min = min(tmy);
    y_max = max(tmy);
    x_min = min(tmx)+1;
    x_max = max(tmx)-1;
    for x = x_min:x_max
        tmp = find(bitget(t_ind(y_min:y_max,x), ROInum)==1);
        tmp(diff([tmp; 2*y_max])==1) = [];
        if mod(numel(tmp),2)==0
            while numel(tmp)>1
                t_ind(y_min-1+(tmp(1)+1:tmp(2)-1),x) = bitset(t_ind(y_min-1+(tmp(1)+1:tmp(2)-1),x), ROInum);
                tmp = tmp(3:end);
            end;
        end;
    end;
    y_min = y_min+1;
    y_max = y_max-1;
    x_min = x_min-1;
    x_max = x_max+1;
    for y = y_min:y_max
        tmp = find(bitget(t_ind(y,x_min:x_max), ROInum)==1);
        tmp(diff([tmp 2*x_max])==1) = [];
        if mod(numel(tmp),2)==0
            while numel(tmp)>1
                t_ind(y,x_min-1+(tmp(1)+1:tmp(2)-1)) = bitset(t_ind(y,x_min-1+(tmp(1)+1:tmp(2)-1)), ROInum);
                tmp = tmp(3:end);
            end;
        end;
    end;
    s_ind(:,handles.pulse) = reshape(t_ind, numel(t_ind),1);

elseif p==2 % wavelength diagram

    mx = xys(1,:);
    my = xys(2,:);
    y_min = min(my);
    y_max = max(my);
    x_min = min(mx);
    x_max = max(mx);

    for n=1:numel(handles.M_X)
        if (handles.ZZ(n)>=x_min)&&(handles.ZZ(n)<=x_max)&&(handles.Z(n)>=y_min)&&(handles.Z(n)<=y_max)
            tx = mx - handles.ZZ(n);
            ty = my - handles.Z(n);
            dx = find(diff(sign(tx)));
            dy = find(diff(sign(ty)));
            y  = ty(dx) - tx(dx).*(ty(dx) - ty(dx+1))./(tx(dx)-tx(dx+1));
            x  = tx(dy) - ty(dy).*(tx(dy) - tx(dy+1))./(ty(dy)-ty(dy+1));
            y = sign(y);
            x = sign(x);
            tmp = mod(numel(find(x<0)),2)+mod(numel(find(x>0)),2)+mod(numel(find(y<0)),2)+mod(numel(find(y>0)),2);
            if tmp>2
                s_ind(n,handles.pulse) = bitset(s_ind(n,handles.pulse), ROInum);
            elseif numel(find(x==0))||numel(find(y==0))
                s_ind(n,handles.pulse) = bitset(s_ind(n,handles.pulse), ROInum);
            end
        end
    end;

elseif p==3  % moments plot

    mx = xys(1,:);
    my = xys(2,:);
    y_min = min(my);
    y_max = max(my);
    x_min = min(mx);
    x_max = max(mx);

    for n=1:numel(handles.M_X)
        if (handles.M_X(n)>=x_min)&&(handles.M_X(n)<=x_max)&&(handles.M_Y(n)>=y_min)&&(handles.M_Y(n)<=y_max)
            tx = mx - handles.M_X(n);
            ty = my - handles.M_Y(n);
            dx = find(diff(sign(tx)));
            dy = find(diff(sign(ty)));
            y  = ty(dx) - tx(dx).*(ty(dx) - ty(dx+1))./(tx(dx)-tx(dx+1));
            x  = tx(dy) - ty(dy).*(tx(dy) - tx(dy+1))./(ty(dy)-ty(dy+1));
            y = sign(y);
            x = sign(x);
            tmp = mod(numel(find(x<0)),2)+mod(numel(find(x>0)),2)+mod(numel(find(y<0)),2)+mod(numel(find(y>0)),2);
            if tmp>2
                s_ind(n,handles.pulse) = bitset(s_ind(n,handles.pulse), ROInum);
            elseif numel(find(x==0))||numel(find(y==0))
                s_ind(n,handles.pulse) = bitset(s_ind(n,handles.pulse), ROInum);
            end
        end
    end;

elseif p==4 % tau-wavelength diagram

    mx = xys(1,:);
    my = xys(2,:);
    y_min = min(my);
    y_max = max(my);
    x_min = min(mx);
    x_max = max(mx);

    for n=1:numel(handles.M_X)
        if (handles.M_X(n)>=x_min)&&(handles.M_X(n)<=x_max)&&(handles.Z(n)>=y_min)&&(handles.Z(n)<=y_max)
            tx = mx - handles.M_X(n);
            ty = my - handles.Z(n);
            dx = find(diff(sign(tx)));
            dy = find(diff(sign(ty)));
            y  = ty(dx) - tx(dx).*(ty(dx) - ty(dx+1))./(tx(dx)-tx(dx+1));
            x  = tx(dy) - ty(dy).*(tx(dy) - tx(dy+1))./(ty(dy)-ty(dy+1));
            y = sign(y);
            x = sign(x);
            tmp = mod(numel(find(x<0)),2)+mod(numel(find(x>0)),2)+mod(numel(find(y<0)),2)+mod(numel(find(y>0)),2);
            if tmp>2
                s_ind(n,handles.pulse) = bitset(s_ind(n,handles.pulse), ROInum);
            elseif numel(find(x==0))||numel(find(y==0))
                s_ind(n,handles.pulse) = bitset(s_ind(n,handles.pulse), ROInum);
            end
        end
    end;

end;


%% update_fFLIM

function update_fFLIM(handles)

pulse   = handles.pulse;
binning = handles.par_binning;
nx      = handles.head.ImgHdr.PixX;
ny      = handles.head.ImgHdr.PixY;
tx      = floor(nx/binning);
ty      = floor(ny/binning);
nch     = size(handles.tcspc,1);
nbin    = size(handles.tcspc,2);
num_PIE = size(handles.tcspc,3);
tau     = handles.head.tau;
binw    = handles.head.tauw./handles.head.tauw(1);

x = handles.head.ImgHdr.X0+(1:binning:handles.head.ImgHdr.PixX)*handles.head.ImgHdr.PixelSize;
y = handles.head.ImgHdr.Y0+(1:binning:handles.head.ImgHdr.PixY)*handles.head.ImgHdr.PixelSize;

xl = [399 440 485 540 580 610 630 699];
yl = [[0 0 0]; [0 0 1]; [0 1 1]; [0 1 0]; [1 1 0]; [1 0.65 0]; [1 0 0]; [1 0 1]];
lambda   = 399:3:699;
spectrum = interp1(xl, yl, lambda);

figure(handles.figMA+pulse);
name = sprintf('update FastFLIM %d', pulse);
set(gcf,'Name',name,'NumberTitle','off');

p1 = 0.06;
p2 = 0.555;
w  = 0.435;

load([handles.pathname handles.timname], 'stim');

if binning > 1
    tmp = stim(1:binning*tx,:,:,:,:);
    tmp = reshape(tmp, binning, tx, ny, nch, nbin, num_PIE);
    tmp = squeeze(mean(tmp,1));
    tmp = tmp(:,1:binning*ty,:,:);
    tmp = reshape(tmp, tx, binning, ty, nch, nbin, num_PIE);
    stim = permute(mean(tmp,2), [1 3 4 5 6 2]);
end

t_m = squeeze(handles.tav(:,:,:,handles.pulse,1)); %#ok<COLND>

W = squeeze(handles.tag(:,:,:,handles.pulse));
V = squeeze(stim(:,:,:,:,handles.pulse)); %#ok<COLND>

W(W<handles.par_Threshold) = 0;

for n = 1:size(W,3)
    W(:,:,n)   = W(:,:,n).*handles.channel(n);
    V(:,:,n,:) = V(:,:,n,:).*handles.channel(n);
end

t_m = sum(t_m.*W,3)./sum(W,3);

if size(W,3)>1
    V = squeeze(sum(V,3));
end

lam = (handles.par_lam_start+handles.par_lam_step.*(-1:size(W,3)));
ts  = handles.par_clam_min: 1: handles.par_clam_max;
tmp = reshape(W, size(W,1).*size(W,2), size(W,3));

S3 = subplot('Position', [p2 p1 w w]);

hold on;
cla
box on;

for l = 1:handles.maxROI
    if handles.ROIlist(l,handles.pulse) == 1
        ind  = bitget(handles.s_ind(:,pulse), l)==1;
        if sum(ind)>0
            s   = squeeze(sum(tmp(ind,:),1));
            s   = [0; s(:); 0]./max(s);
            plot(ts,interp1(lam, s, ts,'cubic'),'-','Color',handles.ROIcolor(l,:));
        end;
    end;
end;
hold off

set(S3,'FontSize',9,'color', [204 204 204]./255);

% axis([10*floor(lind(1)/10) 10*ceil(lind(end)/10) 0 1.1])
axis([handles.par_clam_min handles.par_clam_max 0 1.1])

xlabel('wavelength of spectral mean / nm','FontWeight','bold','FontSize',10,'FontAngle','italic');
ylabel('rel. spectral irradiance','FontWeight','bold','FontSize',10,'FontAngle','italic');

S4 = subplot('Position', [p2 p2 w w]);

hold on;
cla
box on;

for n = 1:handles.maxROI
    if handles.ROIlist(n,handles.pulse)==1
        tmp = reshape(V, numel(handles.s_ind(:,pulse)), nbin);
        ind = bitget(handles.s_ind(:,pulse), n)==1;
        tmp = squeeze(sum(tmp(ind,:),1));
        tmp = tmp./(ones(size(tmp,1))*handles.head.tauw);
        tmp(tmp<0) = 1;
        plot(tau,log10(tmp./max(tmp)),'-','Color',handles.ROIcolor(n,:));
        clear tmp;
    end;
end
hold off;

set(S4,'FontSize',9,'color', [204 204 204]./255);

axis([0 floor(tau(end)) -3 0.1])

xlabel('time / ns','FontWeight','bold','FontSize',10,'FontAngle','italic');
ylabel('log (rel. frequency)','FontWeight','bold','FontSize',10,'FontAngle','italic');


S1 = subplot('Position', [p1 p2 w w]);

data   = squeeze(handles.tag(:,:,:,handles.pulse));

for n = 1:size(data,3)
    data(:,:,n) = data(:,:,n).*handles.channel(n);
end

data   = squeeze(sqrt(sum(data,3)));
tmp    = max(max(data));
data   = data./tmp;

data(data<0)=0;

tmp = repmat(data, [1 1 3]);

lims = [handles.par_tau_min handles.par_tau_max];

val = reshape(t_m,numel(t_m),1); %#ok<COLND>

val(isnan(val))  = 0;
val(val<lims(1)) = 0;
val(val>lims(2)) = lims(2);

k      = 2 + round((val-lims(1))./(lims(2)-lims(1)).*(size(spectrum,1)-2));
k(k<2) = 1;

im = spectrum(k,:);
im = reshape(im,size(tmp,1),size(tmp,2),3);
im = im.*tmp;

alpha = zeros(size(handles.tag,1),size(handles.tag,2));
ovl   = zeros(numel(alpha),3);
for n = 1:handles.maxROI
    if handles.ROIlist(n,handles.pulse) == 1
        tst = bitget(handles.s_ind(:,pulse), n)==1;
        alpha(tst==1) = 1;
        for k = 1:3
            tc  = repmat(handles.ROIcolor(n,k),[size(im,1) size(im,2) 1]);
            ovl(tst,k) = tc(tst);
        end;
    end;
end;
ovl = tmp.*reshape(ovl, size(im,1), size(im,2), 3);

image(x,y,im)
hold on;

colormap(spectrum)

set(S1,'DataAspectRatio', [1,1,1], ...
    'PlotBoxAspectRatio',[1 1 1], ...
    'XDir','normal', ...
    'YDir','reverse', ...
    'FontSize',9,...
    'CLim', lims, ...
    'color', [0 0 0]);

xlabel('x / µm','FontWeight','bold','FontSize',10,'FontAngle','italic');
ylabel('y / µm','FontWeight','bold','FontSize',10,'FontAngle','italic');

tst = 2 + round(((lims(1):.5:lims(2))-lims(1))./(lims(2)-lims(1)).*(size(spectrum,1)-3));
n= 1;
for ttst = lims(1):0.5:lims(2)
    stst(n) = cellstr(sprintf('%4.1f',ttst));
    n = n+1;
end

cca = colorbar;
cca.Label.String = 'lifetime / ns';
%colorbar('CLim', lims, 'YTickMode','manual','XColor',[0 0 0],'YColor',[0 0 0],'Box','off','YTick', tst, ...
%    'FontSize',9, 'YTickLabel',stst,'TickDir','out');

iim2 = image(x,y, ovl);
set(iim2,'AlphaData',alpha);

hold off

S2 = subplot('Position', [p1 p1 w w]);

lims = [399 699];

K = 1:size(W,3);
K = repmat(reshape(K,[1 1 numel(K)]),[size(W,1) size(W,2) 1]);

tZ = W.*K;
tZ = handles.par_lam_start + handles.par_lam_step.*(sum(tZ,3)./sum(W,3)-1.0);

val = tZ;
val = reshape(val,size(val,1)*size(val,2),1);

val(isnan(val))  = 0;
val(val<lims(1)) = 0;
val(val>lims(2)) = lims(2);

k      = 2 + round((val-lims(1))./(lims(2)-lims(1)).*(size(spectrum,1)-3));
k(k<2) = 1;

im = spectrum(k,:);
im = reshape(im,size(tmp,1),size(tmp,2),3);
im = im.*tmp;

image(x,y,im)
hold on;

set(S2,'DataAspectRatio', [1,1,1], ...
    'PlotBoxAspectRatio',[1 1 1], ...
    'XDir','normal', ...
    'YDir','reverse', ...
    'FontSize',9,...
    'CLim', lims, ...
    'color', [0 0 0]);
xlabel('x / µm','FontWeight','bold','FontSize',10,'FontAngle','italic');
ylabel('y / µm','FontWeight','bold','FontSize',10,'FontAngle','italic');

tst = 2 + round(((440:40:680)-lims(1))./(lims(2)-lims(1)).*(size(spectrum,1)-3));

ccb = colorbar;
ccb.Label.String = 'wavelength / nm';
%colorbar('CLim', lims, 'XColor',[0 0 0],'YColor',[0 0 0],'Box','off','YTick', tst, ...
%    'FontSize',9, 'YTickLabel',{'440', '480', '520', '560', '600', '640', '680'},'TickDir','out');

iim3 = image(x,y, ovl);
set(iim3,'AlphaData',alpha);

hold off

%%

function [] = Plot_Results(handles)

list_entries   = get(handles.Pattern_List, 'String');
index_selected = get(handles.Pattern_List, 'Value');
list   = list_entries(index_selected);

if isfield(handles, 'k')
    k = handles.k;
end
res = handles.res;

nx      = size(res,1);
ny      = size(res,2);
num_pat = size(res,3)-1;

scrsz = get(0,'ScreenSize');

fwidth = 800;
fheight = round(2*fwidth/3);
fx = 0.5*(scrsz(3)-fwidth);
fy = 0.5*(scrsz(4)-fheight);

py1 = 0.10;
px1 = 0.06;
wy  = 0.88;
wx  = 0.88;

x = handles.head.ImgHdr.X0+(1:handles.par_binning:handles.head.ImgHdr.PixX)*handles.head.ImgHdr.PixelSize;
y = handles.head.ImgHdr.Y0+(1:handles.par_binning:handles.head.ImgHdr.PixY)*handles.head.ImgHdr.PixelSize;

res = reshape(res,[nx*ny num_pat+1]);
lims = [0 max(max(res))];
if lims(2) <= lims(1)
   lims(2) = lims(1)+1;
end

%map  = [0 0 0; 0 0 1; 0 1 1; 0 1 0; 1 1 0; 1 0 0; 0.5 0 0];
%map1 = interp1(linspace(0, 1, 7), map, linspace(0, 1, 256));

map1 = gray;

res = reshape(res,[nx ny num_pat+1]);

for n = 1:num_pat+1
    if n>num_pat
        name = sprintf('Residuals');
    else
        s   = cell2mat(list(n,:));
        name = sprintf('Pattern: %s', s);
    end

    figure(handles.figPA+n);

    set(gcf,'Position',[fx fy fwidth fheight]);
    set(gcf,'Name',name,'NumberTitle','off');
    R1 = subplot('Position', [px1 py1 wx wy]);
    colormap(map1);
    hold on;
    cla
    imagesc(x,y,res(:,:,n));
%    imagesc(x,y,res(:,:,n).*lims(2)./max(max(res(:,:,n))));

    set(R1,'DataAspectRatio', [1,1,1], ...
        'PlotBoxAspectRatio',[1 1 1], ...
        'XDir','normal', ...
        'YDir','reverse', ...
        'FontSize',9,...
        'CLim', [0 max(max(res(:,:,n)))], ... % Change
        'color', [0 0 0]);

    xlabel('x (µm)','FontWeight','bold','FontSize',10,'FontAngle','italic');
    ylabel('y (µm)','FontWeight','bold','FontSize',10,'FontAngle','italic');
    m = n;

    colorbar('YColor',[0 0 0])
end

if exist('k','var')
    if handles.match == 1
        name = sprintf('%s: lifetime fit', cell2mat(list(k,:)));            
    elseif handles.match == 2 || handles.match ==3
        name = sprintf('FRET Efficiency');
    end
    figure(handles.figPA+m+1);

    set(gcf,'Position',[fx fy fwidth fheight]);
    set(gcf,'Name',name,'NumberTitle','off');
    R1 = subplot('Position', [px1 py1 wx wy]);

    xl = [399 460 550 640 699];% 699];
    yl = [[0 0 1]; [0 1 1]; [0 1 0]; [1 1 0]; [1 0 0]];%; [1 0 1]];
    lambda   = 399:3:699;
    spectrum = interp1(xl, yl, lambda);

    %intens = squeeze(sqrt(sum(handles.tag,3)));
    intens = squeeze(squeeze(res(:,:,k)));
    intens = intens./max(max(intens));

    tmp = repmat(intens(:,:), [1 1 3]);

    tmp_lims = [min(handles.vParam) max(handles.vParam)];

    val = handles.vParamDist;

    val(isnan(val))  = 0;
    val(val<tmp_lims(1)) = 0;
    val(val>tmp_lims(2)) = tmp_lims(2);

    g      = 2 + round((val-tmp_lims(1))./(tmp_lims(2)-tmp_lims(1)).*(size(spectrum,1)-2));
    g(g<2) = 1;

    im = spectrum(g,:);
    im = reshape(im,size(tmp,1),size(tmp,2),3);
    im = im.*tmp;

    image(x,y,im)

    colormap(spectrum)

    set(R1,'DataAspectRatio', [1,1,1], ...
        'PlotBoxAspectRatio',[1 1 1], ...
        'XDir','normal', ...
        'YDir','reverse', ...
        'FontSize',9,...
        'CLim', tmp_lims, ...
        'color', [0 0 0]);

    xlabel('x (µm)','FontWeight','bold','FontSize',10,'FontAngle','italic');
    ylabel('y (µm)','FontWeight','bold','FontSize',10,'FontAngle','italic');

    if handles.match == 1
        tst = 2 + round(((tmp_lims(1):.5:tmp_lims(2))-tmp_lims(1))./(tmp_lims(2)-tmp_lims(1)).*(size(spectrum,1)-3));
        n= 1;
        for ttst = tmp_lims(1):0.5:tmp_lims(2)
            stst(n) = cellstr(sprintf('%4.1f',ttst));
            n = n+1;
        end
        colorbar('CLim', tmp_lims, 'YTickMode','manual','XColor',[0 0 0],'YColor',[0 0 0],'Box','off','YTick', tst, ...
            'FontSize', 9,'TickDir','out', 'YTickLabel',stst);
    elseif handles.match == 2 || handles.match == 3
        tst = 2 + round(((tmp_lims(1):5:tmp_lims(2))-tmp_lims(1))./(tmp_lims(2)-tmp_lims(1)).*(size(spectrum,1)-3));
        n= 1;
        for ttst = tmp_lims(1):5:tmp_lims(2)
            stst(n) = cellstr(sprintf('%4.1f',ttst));
            n = n+1;
        end
        colorbar('CLim', tmp_lims, 'YTickMode','manual','XColor',[0 0 0],'YColor',[0 0 0],'Box','off','YTick', tst, ...
            'FontSize', 9,'TickDir','out', 'YTickLabel',stst);
    end
end


if handles.match == 2 || handles.match == 3

    Binding = 100.*res(:,:,1)./(res(:,:,1)+res(:,:,2));

    name = sprintf('Binding %d', cell2mat(list(k,:)));            

    figure(handles.figPA+m+2);

    set(gcf,'Position',[fx fy fwidth fheight]);
    set(gcf,'Name',name,'NumberTitle','off');
    R1 = subplot('Position', [px1 py1 wx wy]);

    xl = [399 460 550 640 699];% 699];
    yl = [[0 0 1]; [0 1 1]; [0 1 0]; [1 1 0]; [1 0 0]];%; [1 0 1]];
    lambda   = 399:3:699;
    spectrum = interp1(xl, yl, lambda);

    %intens = squeeze(sqrt(sum(handles.tag,3)));
    intens = squeeze(res(:,:,1)+res(:,:,2));
    intens = intens./max(max(intens));

    tmp = repmat(intens(:,:), [1 1 3]);

    tmp_lims = [0 100];
    
    Binding = squeeze(Binding);
    Binding(isnan(Binding))  = tmp_lims(2);
    Binding(Binding<tmp_lims(1)) = 0;
    Binding(Binding>tmp_lims(2)) = tmp_lims(2);

    g      = 2 + round((Binding-tmp_lims(1))./(tmp_lims(2)-tmp_lims(1)).*(size(spectrum,1)-2));
    g(g<2) = 1;

    im = spectrum(g,:);
    im = reshape(im,size(tmp,1),size(tmp,2),3);
    im = im.*tmp;

    image(x,y,im)

    colormap(spectrum)

    set(R1,'DataAspectRatio', [1,1,1], ...
        'PlotBoxAspectRatio',[1 1 1], ...
        'XDir','normal', ...
        'YDir','reverse', ...
        'FontSize',9,...
        'CLim', tmp_lims, ...
        'color', [0 0 0]);

    xlabel('x (µm)','FontWeight','bold','FontSize',10,'FontAngle','italic');
    ylabel('y (µm)','FontWeight','bold','FontSize',10,'FontAngle','italic');


    tst = 2 + round(((tmp_lims(1):10:tmp_lims(2))-tmp_lims(1))./(tmp_lims(2)-tmp_lims(1)).*(size(spectrum,1)-3));
    n= 1;
    for ttst = tmp_lims(1):10:tmp_lims(2)
        stst(n) = cellstr(sprintf('%4.1f',ttst));
        n = n+1;
    end
    colorbar('CLim', tmp_lims, 'YTickMode','manual','XColor',[0 0 0],'YColor',[0 0 0],'Box','off','YTick', tst, ...
        'FontSize', 9,'TickDir','out', 'YTickLabel',stst);
end



name = sprintf('Composite');

im = MakeComposite('data',res(:,:,1:num_pat),'list',list);

figure(handles.figPA+num_pat+3);

set(gcf,'Position',[fx fy fwidth fheight]);
set(gcf,'Name',name,'NumberTitle','off');
R1 = subplot('Position', [px1 py1 wx wy]);

hold on;

cla

image(x,y,reshape(im,[nx ny 3]))

set(R1,'DataAspectRatio', [1,1,1], ...
    'PlotBoxAspectRatio',[1 1 1], ...
    'XDir','normal', ...
    'YDir','reverse', ...
    'FontSize',9,...
    'CLim', lims, ...
    'color', [0 0 0]);

xlabel('x (µm)','FontWeight','bold','FontSize',10,'FontAngle','italic');
ylabel('y (µm)','FontWeight','bold','FontSize',10,'FontAngle','italic');

hold off

s = {};
tst = 1;
for n = 1:handles.maxROI
    if handles.ROIlist(n, handles.pulse)==1
        s(tst) = {sprintf('ROI %d_%d',n,handles.pulse)};
        tst = tst+1;
    end
end
set(handles.listbox2,'Value',1);
set(handles.listbox2,'String',s);
listbox2_Callback(handles.listbox2, [], handles);

plot_graphs(handles)
