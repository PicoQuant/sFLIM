function varargout = Set_PatternOptions(varargin)
% SET_PATTERNOPTIONS M-file for Set_PatternOptions.fig
%      SET_PATTERNOPTIONS, by itself, creates a new SET_PATTERNOPTIONS or raises the existing
%      singleton*.
%
%      H = SET_PATTERNOPTIONS returns the handle to a new SET_PATTERNOPTIONS or the handle to
%      the existing singleton*.
%
%      SET_PATTERNOPTIONS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SET_PATTERNOPTIONS.M with the given input arguments.
%
%      SET_PATTERNOPTIONS('Property','Value',...) creates a new SET_PATTERNOPTIONS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Set_PatternOptions_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Set_PatternOptions_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Set_PatternOptions

% Last Modified by GUIDE v2.5 19-May-2014 21:40:09

% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Set_PatternOptions_OpeningFcn, ...
                   'gui_OutputFcn',  @Set_PatternOptions_OutputFcn, ...
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


% --- Executes just before Set_PatternOptions is made visible.
function Set_PatternOptions_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Set_PatternOptions (see VARARGIN)

% Choose default command line output for Set_PatternOptions
handles.output = hObject;

load Get_Params;
handles.par_lam_start = lam_start;
handles.par_lam_step  = lam_step;
handles.par_clam_min  = clam_min;
handles.par_clam_max  = clam_max;
handles.par_num_PIE   = num_PIE;
handles.ROIcolor      = ROIcolor;

handles.channel = 1;


if length(varargin) < 3
    uiwait(errordlg('No pattern given','Data missing'))
elseif (length(varargin) == 3 && strcmpi(varargin{1},'file'))
    handles.file = cell2mat(varargin(2));
    load(cell2mat(varargin(2)), 'pat', 'list', 'params');
    handles.pat      = pat;
    handles.List     = list;
    handles.params   = params;
    handles.filename = cell2mat(varargin(3));
    handles.match    = 0;
elseif (length(varargin) == 4 && strcmpi(varargin{1},'file'))
    handles.file = cell2mat(varargin(2));
    load(cell2mat(varargin(2)), 'pat', 'list', 'params');
    handles.pat      = pat;
    handles.List     = list;
    handles.params   = params;
    handles.filename = cell2mat(varargin(3));
    
    patfam   = varargin(4);
    handles.vTau = cell2mat(patfam{1,1}{1,1}(end));
    tmp = patfam{1,1}{1,1}(1:end-1);
    handles.fam = cat(4,tmp{1,1},tmp{1,2});
    for n = 1:size(tmp,2)-2
    handles.fam = cat(4,handles.fam,tmp{1,n+2});
    end
    handles.fam = permute(handles.fam,[4,1,2,3]);
    handles.match = 1;
end

tmp = strfind(handles.filename,'\');
handles.pathname = handles.filename(1:tmp(end));
handles.filename = handles.filename(tmp(end)+1:end);
[handles.head, handles.tag, handles.tcspc, handles.IRF, handles.timname,...
    handles.AP, handles.DC] = load_data([handles.pathname handles.filename]);

handles.nch   = size(handles.tcspc,1);
handles.binw  = handles.head.tauw(:)./handles.head.tauw(1);
handles.binw  = squeeze(repmat(handles.binw,[num_PIE 1]));
handles.taush = handles.head.tau;
handles.tau   = handles.head.tau(:)*ones(1,handles.par_num_PIE);
handles.tmp   = ones(size(handles.tau,1),1)*((0:handles.par_num_PIE-1).*handles.tau(end,1));
handles.tau   = reshape(handles.tau+handles.tmp, [numel(handles.tau) 1]);
handles.m0    = handles.head.m0;


% Update handles structure
guidata(hObject, handles);

load Get_Params;


col_name(1,1) = {'spectrum'};
col_name(1,2) = {'all'};

lambda = (lam_start+lam_step.*(0:size(handles.pat,2)-1));

for n = 1:size(handles.pat,2)
    col_name(1,2+n) = {sprintf('%03d',lambda(n))};
end

for n = 1:num_PIE
    row_name(1,n) = {sprintf(' Pulse %01d',n)};
end

data = ' [X]';
data = repmat({data}, num_PIE, size(handles.pat,2)+2);

lambda = repmat([inf inf lambda], num_PIE, 1);
tmp = repmat(lamPIE(1:num_PIE)', 1, size(handles.pat,2)+2)+lam_step;

ind = lambda<=tmp;

data(ind) = {' [  ]'};

set(handles.uitable1,'Data',data);
set(handles.uitable1,'RowName',row_name);
set(handles.uitable1,'ColumnName',col_name);
wt = num2cell(25.*ones(1,size(data,2)));
wt(1) = {60};
wt(2) = {30};
set(handles.uitable1,'ColumnWidth',wt);
set(handles.uitable1,'ColumnEditable',false(1, size(data,2)));

for n = 1:size(handles.pat,3)-handles.head.m0+1
    col_name(1,2+n) = {sprintf('%02d',n)};
end

data = ' [X]';
data = repmat({data}, num_PIE, size(handles.pat,3)-handles.head.m0+3);

col_name(1,1) = {'decay'};
set(handles.uitable2,'Data',data);
set(handles.uitable2,'RowName',row_name);
set(handles.uitable2,'ColumnName',col_name);
wt = num2cell(25.*ones(1,size(data,2)));
wt(1) = {50};
wt(2) = {30};
set(handles.uitable2,'ColumnWidth',wt);
set(handles.uitable2,'ColumnEditable',false(1, size(data,2)));


handles.ROIcolor = ROIcolor;

axes(handles.axes1)
hold on;
cla
set(handles.axes1, ...
    'Box','on', ...
    'XDir','normal', ...
    'YDir','normal', ...
    'FontSize',9,...
    'Color',[50 50 50]./255);
xlabel('wavelength of spectral mean / nm','FontWeight','bold','FontSize',10,'FontAngle','italic');
ylabel('rel. spectral irradiance','FontWeight','bold','FontSize',10,'FontAngle','italic');
axis([clam_min clam_max 0 1.1]);

s   = squeeze(sum(sum(pat,4),3));
s(isnan(s)) = 0;
lam = (lam_start+lam_step.*(-1:size(s,2)));
s = [zeros(size(s,1),1) s zeros(size(s,1),1)];
ts = lam_start-lam_step: 1: lam_start+lam_step*(1+size(handles.pat,2));
cs = interp1(lam(:), s', ts','cubic');
for a = 1:size(s,1)
    plot(ts,cs(:,a)./max(max(cs)),'Color',ROIcolor(1+mod(a-1,size(ROIcolor,1)),:),'LineWidth',3); %#ok<COLND>
end
hold off

handles.lam  = lam;
handles.ts   = ts;


axes(handles.axes2)
hold on;
cla
set(handles.axes2, ...
    'Box','on', ...
    'XDir','normal', ...
    'YDir','normal', ...
    'FontSize',9,...
    'Color',[50 50 50]./255);
xlabel('time / ns','FontWeight','bold','FontSize',10,'FontAngle','italic');
ylabel('rel. irradiance','FontWeight','bold','FontSize',10,'FontAngle','italic');
axis([0 ceil(handles.taush(end)) -4 0.1]);

s = squeeze(sum(sum(pat,4),2));
s = s./repmat(handles.head.tauw, [size(s,1) 1]);
s(isnan(s)) = 0;

ts = handles.taush;
for a = 1:size(s,1)
   plot(ts,log10(s(a,:)./max(max(s))),'Color',ROIcolor(1+mod(a-1,size(ROIcolor,1)),:),'LineWidth',3); %#ok<COLND>
end
hold off


guidata(hObject, handles);


% col_name = list;
% row_name = {'rel. Intensity'};
% set(handles.uitable2,'Data',squeeze(handles.res(:,1,1,1))');
% set(handles.uitable2,'RowName',row_name);
% set(handles.uitable2,'ColumnName',col_name);

guidata(hObject, handles);


% UIWAIT makes Set_PatternOptions wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Set_PatternOptions_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Cancel_button.
function Cancel_button_Callback(hObject, eventdata, handles)
% hObject    handle to Cancel_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

params = zeros(size(handles.pat,4), 2+size(handles.pat,2)+size(handles.pat,3));

save(handles.file, 'params', '-APPEND');

guidata(hObject, handles);
close

% --- Executes on button press in OK_button.
function OK_button_Callback(hObject, eventdata, handles)
% hObject    handle to OK_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% load('Get_Params.mat');

data1 = get(handles.uitable1, 'Data');
data2 = get(handles.uitable2, 'Data');

n1 = size(data1,2)-2;
n2 = size(data2,2)-2;

params = zeros(size(handles.pat,4), 2+n1+size(handles.pat,3));

params(:,1) = strcmpi(data1(:,1),{' [X]'});
params(:,2) = strcmpi(data2(:,1),{' [X]'});
params(:,2+(1:n1)) = strcmpi(data1(:,3:end),{' [X]'});
params(:,2+n1+(1:handles.head.m0-1)) = repmat(strcmpi(data2(:,3),{' [X]'}),[1 handles.head.m0-1]);
params(:,1+n1+handles.head.m0+(1:n2)) = strcmpi(data2(:,3:end),{' [X]'});

save(handles.file, 'params', '-APPEND');
guidata(hObject, handles);
close


% --- Executes when selected cell(s) is changed in uitable1.
function uitable1_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)

load Get_Params;

data  = get(hObject,'Data');
data1 = get(handles.uitable2,'Data');

save('data','data')
if numel(eventdata.Indices)>0
    if strcmpi(data(eventdata.Indices(1),eventdata.Indices(2)),{' [X]'})
        if eventdata.Indices(2) == 1
           data(:,eventdata.Indices(2)) = {' [  ]'};
           data1(:,1) = {' [X]'};
        elseif eventdata.Indices(2) == 2
           data(eventdata.Indices(1),2:end) = {' [  ]'};
        else
           data(eventdata.Indices(1),eventdata.Indices(2)) = {' [  ]'};
        end
    else
        if eventdata.Indices(2) < 2
           data(:,eventdata.Indices(2)) = {' [X]'};
        elseif eventdata.Indices(2) == 2
           data(eventdata.Indices(1),2:end) = {' [X]'};
        else
           data(eventdata.Indices(1),eventdata.Indices(2)) = {' [X]'};
        end
    end
    set(hObject,'Selected','off');
    set(hObject,'Data',data);
    set(handles.uitable2,'Data',data1);
end

guidata(hObject, handles);
update_plot(handles)


% --- Executes when selected cell(s) is changed in uitable2.
function uitable2_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitable2 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
load Get_Params;

data = get(hObject,'Data');
data1 = get(handles.uitable1,'Data');

save('data','data')
if numel(eventdata.Indices)>0
    if strcmpi(data(eventdata.Indices(1),eventdata.Indices(2)),{' [X]'})
        if eventdata.Indices(2) == 1
           data(:,1)  = {' [  ]'};
           data1(:,1) = {' [X]'};
        elseif eventdata.Indices(2) == 2
           data(eventdata.Indices(1),2:end) = {' [  ]'};
        else
           data(eventdata.Indices(1),eventdata.Indices(2)) = {' [  ]'};
        end
    else
        if eventdata.Indices(2) < 2
           data(:,eventdata.Indices(2)) = {' [X]'};
        elseif eventdata.Indices(2) == 2
           data(eventdata.Indices(1),2:end) = {' [X]'};
        else
           data(eventdata.Indices(1),eventdata.Indices(2)) = {' [X]'};
        end
    end
    set(hObject,'Selected','off');
    set(hObject,'Data',data);
    set(handles.uitable1,'Data',data1);
end

guidata(hObject, handles);
update_plot(handles)


function update_plot(handles)

style = [{'-'};{'--'};{'-.'};{':'}];

num_Pat = size(handles.pat,1);
num_PIE = size(handles.pat,4);

data1 = get(handles.uitable1,'Data');
data2 = get(handles.uitable2,'Data');

use_spec = strcmpi(data1(:,1),{' [X]'});
use_dec  = strcmpi(data2(:,1),{' [X]'});

s   = permute(sum(handles.pat,3),[1 2 4 3]);

s   = reshape(permute(s,[3 1 2]),[size(s,1)*size(s,3) size(s,2)]);
m   = max(max(s));
ind = repmat(strcmpi(data1(:,3:end),{' [X]'}),[size(handles.pat,1) 1]);
s   = s.*ind;

s = [zeros(size(s,1),1) s zeros(size(s,1),1)];

hold(handles.axes1,'on')
cla(handles.axes1)

for a = 1:num_Pat
    for b = 1:num_PIE
        ind = (a-1)*num_PIE+b;
        color = handles.ROIcolor(1+mod(a-1,size(handles.ROIcolor,1)),:);
        if use_spec(b)
            cs = interp1(handles.lam(:), s(ind,:)', handles.ts,'cubic');
            plot(handles.axes1,handles.ts,cs./m,'Color',color,'LineWidth',3,'LineStyle',cell2mat(style(b))); %#ok<COLND>
        else
            bar(handles.axes1,handles.lam(:)+(a-1).*2, (s(ind,:)>0)*(5-b)./4 ,'BarWidth',0.2*b,'FaceColor','none','LineWidth',2,'EdgeColor',color);
        end
    end
end
hold(handles.axes1,'off')
axis([handles.par_clam_min handles.par_clam_max 0 1.1]);
%  axis([handles.lam(2) handles.lam(end-1) 0 1])


s   = permute(sum(handles.pat,2),[1 3 4 2]); 
if size(handles.pat,1) == 1
    s = shiftdim(s,-1);
end

s   = s./repmat(handles.head.tauw,[size(s,1) 1 size(s,3)]);
s   = reshape(permute(s,[3 1 2]),[size(s,1)*size(s,3) size(s,2)]);
m   = max(max(s));
ind = strcmpi(data2(:,3:end),{' [X]'});
ind = [repmat(ind(:,1),[1 handles.head.m0]) ind(:,2:end)];
ind = repmat(ind,[size(handles.pat,1) 1]);
s   = s.*ind;

hold(handles.axes2,'on')
cla(handles.axes2)

for a = 1:num_Pat
    for b = 1:num_PIE
        ind = (a-1)*num_PIE+b;
        color = handles.ROIcolor(1+mod(a-1,size(handles.ROIcolor,1)),:);
        if use_dec(b)
            plot(handles.axes2,handles.taush,log10(s(ind,:)./m),'Color',color,'LineWidth',3,'LineStyle',cell2mat(style(b))); %#ok<COLND>
        else
            bar(handles.axes2,handles.taush(:)+(a-1)/8, -(s(ind,:)>0).*4 ,'BarWidth',0.8*b,'FaceColor','none','LineWidth',2,'EdgeColor',color);
        end
    end
end
hold(handles.axes2,'off')
axis([0 ceil(handles.taush(end)) -4 0.1]);
