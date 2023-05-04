function varargout = Set_PatternOptionsVar(varargin)
% SET_PATTERNOPTIONSVAR M-file for Set_PatternOptionsVar.fig
%      SET_PATTERNOPTIONSVAR, by itself, creates a new SET_PATTERNOPTIONSVAR or raises the existing
%      singleton*.
%
%      H = SET_PATTERNOPTIONSVAR returns the handle to a new SET_PATTERNOPTIONSVAR or the handle to
%      the existing singleton*.
%
%      SET_PATTERNOPTIONSVAR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SET_PATTERNOPTIONSVAR.M with the given input arguments.
%
%      SET_PATTERNOPTIONSVAR('Property','Value',...) creates a new SET_PATTERNOPTIONSVAR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Set_PatternOptionsVar_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Set_PatternOptionsVar_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Set_PatternOptionsVar

% Last Modified by GUIDE v2.5 13-Jan-2023 16:12:30

% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Set_PatternOptionsVar_OpeningFcn, ...
                   'gui_OutputFcn',  @Set_PatternOptionsVar_OutputFcn, ...
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


% --- Executes just before Set_PatternOptionsVar is made visible.
function Set_PatternOptionsVar_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Set_PatternOptionsVar (see VARARGIN)

% Choose default command line output for Set_PatternOptionsVar
handles.output = hObject;

load Get_Params;
handles.par_lam_start = lam_start;
handles.par_lam_step  = lam_step;
handles.par_clam_min  = clam_min;
handles.par_clam_max  = clam_max;
handles.par_num_PIE   = num_PIE;
handles.ROIcolor      = ROIcolor;

handles.channel = 1;


if length(varargin) < 4
    uiwait(errordlg('No pattern given','Data missing'))
elseif (length(varargin) == 4 && strcmpi(varargin{1},'file'))
    handles.file = cell2mat(varargin(2));
    load(cell2mat(varargin(2)), 'pat', 'list', 'params');
    handles.pat      = pat;
    handles.List     = list;
    handles.params   = params;
    handles.filename = cell2mat(varargin(3));
    handles.stim   = cell2mat(varargin(4));
    handles.match = 0;
elseif (length(varargin) == 5 && strcmpi(varargin{1},'file'))
    handles.file = cell2mat(varargin(2));
    load(cell2mat(varargin(2)), 'pat', 'list', 'params');
    handles.pat      = pat;
    handles.List     = list;
    handles.params   = params;
    handles.filename = cell2mat(varargin(3));
    handles.stim     = cell2mat(varargin(5)); 
    
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

col_name(1,1) = {'decay'};
col_name(1,2) = {'spectrum'};
col_name(1,3) = {'all on/off'};
lambda = (lam_start+lam_step.*(0:size(handles.pat,2)-1));

for n = 1:size(handles.pat,2)
    col_name(1,3+n) = {sprintf('%03d',lambda(n))};
end

for n = 1:size(handles.pat,1)
    row_name(1,1+(n-1)*num_PIE) = handles.List(n);
end

data = ' [X]';
data = repmat({data},size(handles.pat,1)*num_PIE, size(handles.pat,2)+3);

lambda = repmat([inf inf inf lambda], size(handles.pat,1)*num_PIE, 1);
tmp = repmat(lamPIE(1:num_PIE)', size(handles.pat,1), size(handles.pat,2)+3)+lam_step;

%ind = lambda<=tmp; %Only wavelengths greater than tmp are checked with X

ind = lambda<=0; %All wavelengths are checked with X

data(ind) = {' [  ]'};

set(handles.uitable1,'Data',data);
set(handles.uitable1,'RowName',row_name);
set(handles.uitable1,'ColumnName',col_name);
wt = num2cell(25.*ones(1,size(pat,2)+2));
wt(1) = {40};
wt(2) = {50};
wt(3) = {60};
set(handles.uitable1,'ColumnWidth',wt);
set(handles.uitable1,'ColumnEditable',false(1, size(data,2)));

handles.ROIcolor = ROIcolor;

guidata(hObject, handles);

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
lam = (lam_start+lam_step.*(-1:size(s,2)));
s = [zeros(size(s,1),1) s zeros(size(s,1),1)];
ts = lam_start-lam_step: 1: lam_start+lam_step*(1+size(handles.pat,2));
cs = interp1(lam(:), s', ts','cubic');
for a = 1:size(s,1)
    plot(ts,cs(:,a)./max(max(cs)),'Color',ROIcolor(mod(a,size(ROIcolor,1)),:),'LineWidth',3); %#ok<COLND>
end
hold off

handles.lam      = lam;
handles.ts       = ts;


% The used Pattern in handles.list are matched to the data contained in
% handles.stim in the lsqnonneg() sense. For the first a approach stim is
% summed up for all pixels and all channels. If one (one at most) 
% contains a family of patterns the best family member is chosen. 


tmp = {};
lam = (handles.par_lam_start+handles.par_lam_step.*(0:handles.nch-1));
for n = 1 : size(handles.tcspc,1)
    tmp(n) = {sprintf(' %s nm',num2str(lam(n)))};
end
set(handles.wavelength, 'String',tmp);

zoom on

for n = 1:numel(handles.List)
    name = handles.List{n};
    name(name == ' ') = '';
    if handles.match == 1 && strcmp(name(end-2:end),'(v)')
        k = n;
        handles.VaryEl = k;
    end
end

[handles.M,handles.res] = Match_Fun(handles);

col_name = list;
row_name = {'rel. Intensity'};
set(handles.uitable2,'Data',squeeze(handles.res(:,1,1,1))');
set(handles.uitable2,'RowName',row_name);
set(handles.uitable2,'ColumnName',col_name);

guidata(hObject, handles);


% UIWAIT makes Set_PatternOptionsVar wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Set_PatternOptionsVar_OutputFcn(hObject, eventdata, handles) 
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
guidata(hObject, handles);
close

% --- Executes on button press in OK_button.
function OK_button_Callback(hObject, eventdata, handles)
% hObject    handle to OK_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

load Get_Params;

params = zeros(size(handles.pat,1)*num_PIE, size(handles.pat,2)+4);

if get(handles.Is_Ortho,'Value')
    params(:,end) = 1;
end

data = get(handles.uitable1, 'Data');
params(:,1:end-1) = strcmpi(data,{' [X]'});

save(handles.file, 'params', '-APPEND');
guidata(hObject, handles);
close

% --- Executes on checkbox press in Is_Ortho_button.
function Is_Ortho_Callback(hObject, eventdata, handles)
% hObject    handle to Cancel_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tmp = get(hObject,'Value');
set(hObject,'Value',tmp);
guidata(hObject, handles);


% --- Executes when selected cell(s) is changed in uitable1.
function uitable1_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)

load Get_Params;

style = [{'-'};{'--'};{'-.'};{':'}];

data = get(hObject,'Data');
save('data','data')
if numel(eventdata.Indices)>0
    if strcmpi(data(eventdata.Indices(1),eventdata.Indices(2)),{' [X]'})
        if eventdata.Indices(2) < 3
           data(:,eventdata.Indices(2)) = {' [  ]'};
           data(:,1+mod(eventdata.Indices(2),2)) = {' [X]'};
        elseif eventdata.Indices(2) == 3
           data(eventdata.Indices(1),3:end) = {' [  ]'};
        else
           data(eventdata.Indices(1),eventdata.Indices(2)) = {' [  ]'};
        end
    else
        if eventdata.Indices(2) < 3
           data(:,eventdata.Indices(2)) = {' [X]'};
        elseif eventdata.Indices(2) == 3
           data(eventdata.Indices(1),3:end) = {' [X]'};
        else
           data(eventdata.Indices(1),eventdata.Indices(2)) = {' [X]'};
        end
    end
    set(hObject,'Selected','off');
    set(hObject,'Data',data);
end

axes(handles.axes1)
hold on;
cla

s   = squeeze(sum(handles.pat,3));
s   = reshape(permute(s,[3 1 2]),[size(s,1)*size(s,3) size(s,2)]);
m   = max(max(s));
s   = s.*strcmpi(data(:,4:end),{' [X]'});
use_spec = strcmpi(data(:,2),{' [X]'});

s = [zeros(size(s,1),1) s zeros(size(s,1),1)];
for a = 1:size(handles.pat,1)
    for b = 1:num_PIE
        ind = (a-1)*num_PIE+b;
        color = handles.ROIcolor(mod(a,size(handles.ROIcolor,1)),:);
        if use_spec(ind)
            cs = interp1(handles.lam(:), s(ind,:)', handles.ts,'cubic');
            plot(handles.ts,cs./m,'Color',color,'LineWidth',3,'LineStyle',cell2mat(style(b))); %#ok<COLND>
        else
            bar(handles.lam(:)+(a-1)*(lam_step)/8, (s(ind,:)>0)*(5-b)/4 ,'BarWidth',0.2*b,'FaceColor','none','LineWidth',2,'EdgeColor',color);
        end
    end
end
hold off
guidata(hObject, handles);


% --- Executes on button press in Match_button.
function Match_button_Callback(hObject, eventdata, handles)
% hObject    handle to Match_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.M,handles.res] = Match_Fun(handles);

guidata(hObject, handles);




% --- Executes on selection change in wavelength.
function wavelength_Callback(hObject, eventdata, handles)
% hObject    handle to wavelength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

update_plot(handles);

guidata(hObject, handles);

% Hints: contents = cellstr(get(hObject,'String')) returns wavelength contents as cell array
%        contents{get(hObject,'Value')} returns selected item from wavelength


% --- Executes during object creation, after setting all properties.
function wavelength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wavelength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% match

function [M,res] = Match_Fun(handles)

pat = handles.pat;

num_pat = size(pat,1);
nch     = size(pat,2);
nbin    = size(pat,3);
num_PIE = size(pat,4);
numpix  = size(handles.stim,1);

load Get_Params;

% Here the choice given by the uitable1 is considered, giving
% the options of matching with spectrum, decay or selected channels
% information.

params = zeros(size(handles.pat,1)*num_PIE, size(handles.pat,2)+4);
if get(handles.Is_Ortho,'Value')
    params(:,end) = 1;
end
data = get(handles.uitable1, 'Data');
params(:,1:end-1) = strcmpi(data,{' [X]'});
    
mode = sum(params(1,1:2).*[1 2],2);

ch_ind = params(:,4:end-1) == 1;
ch_ind = permute(reshape(ch_ind,[num_PIE num_pat nch]),[2 3 1]);
ch_ind = permute(repmat(ch_ind, [1 1 1 nbin]), [1 2 4 3]);
ch_ind = sum(ch_ind)>0;

if mode == 3                                         % use spectrum and decay
    tmp = handles.stim.*repmat(ch_ind, [numpix 1]);
    tmp = squeeze(sum(tmp,1));
    tmp = reshape(tmp,[1 nch*nbin*num_PIE]);
elseif mode == 1                                      % use decay only
    tmp = sum(handles.stim.*repmat(ch_ind, [numpix 1]),2);
    tmp = squeeze(sum(tmp,1));
    tmp = reshape(tmp,[1 nbin*num_PIE]);
else                                                  % use spectrum only
    tmp = sum(handles.stim.*repmat(ch_ind, [numpix 1]),3);
    tmp = squeeze(sum(tmp,1));
    tmp = reshape(tmp,[1 nch*num_PIE]);
end

ch_ind = repmat(ch_ind, [num_pat 1 1 1]);


if isfield(handles, 'VaryEl')
    k = handles.VaryEl;
    vTau = handles.vTau;
    fam = handles.fam;
    for n = 1:size(fam,1)
        M = pat;
        M(k,:,:,:) = fam(n,:,:,:);
        res = ampfit(M, tmp, mode, ch_ind);
    end

else
    M = pat;
    res = ampfit(M, tmp, mode, ch_ind);
end


res = repmat(res(1,end-1),[num_pat nch nbin num_PIE]);
M = M.*res;
M = reshape(M, [num_pat nch nbin*num_PIE]);
handles.M = M;
handles.res = res;

% show and plot results in UI

set(handles.uitable2,'Data',squeeze(res(:,1,1,1))');

update_plot(handles);

%% pat(tau)

function [res] = ampfit(M, tmp, mode, ch_ind)

num_pat = size(M,1);
nch     = size(M,2);
nbin    = size(M,3);
num_PIE = size(M,4);

res = zeros(1, num_pat+1);

M      = M.*ch_ind;
           
if mode == 3                                         % use spectrum and decay
    M   = reshape(M,[num_pat nch*nbin*num_PIE]);
elseif mode == 1                                      % use decay only
    M   = sum(M,2);
    M   = reshape(M,[num_pat nbin*num_PIE]);
else                                                  % use spectrum only
    M   = sum(M,3);
    M   = reshape(M,[num_pat nch*num_PIE]);
end

for j=1:num_pat
    M(j,:) = M(j,:)/sum(M(j,:));  % normalization
end

res(1,1:end-1)  = lsqnonneg((M)', tmp');
res(1,end)    = sum(tmp) - sum(res);

%% plot

function [] = update_plot(handles)

chan = get(handles.wavelength,'Value');
if chan ~= handles.channel
    handles.channel = chan;
end;
nch = size(handles.M,2);

stim = squeeze(sum(handles.stim,1));
stim = reshape(stim,[size(stim,1) size(stim,2)*size(stim,3)]);

axes(handles.axes2)
cla;
hold on

% plot(handles.tau, log10(squeeze(stim(chan,:))./handles.binw'),'.')
% for n = 1:size(handles.M,1)
%     plot(handles.tau,log10(squeeze(handles.M(n,chan,:,:))./handles.binw),'-','Color',handles.ROIcolor(n,:));
% end
% plot(handles.tau,log10(squeeze(sum(handles.M(:,chan,:,:),1))./handles.binw),'r--','LineWidth',3)

hold off

% axis([0 floor(handles.tau(end)) 0 max(log10(squeeze(stim(round(nch/2),:,:))./handles.binw'))+1])
% xlabel('time / ns','FontWeight','bold','FontSize',10,'FontAngle','italic');
% ylabel('log (rel. frequency)','FontWeight','bold','FontSize',10,'FontAngle','italic');
% legend(['data',handles.List','match'],'FontSize',8)

resid = squeeze(stim(chan,:))' - squeeze(sum(handles.M(:,chan,:,:),1));

axes(handles.axes3)
cla;
hold on

% plot(handles.tau, log10(resid),'r-')

hold off

% axis([0 floor(handles.tau(end)) -5 5])
% xlabel('time / ns','FontWeight','bold','FontSize',10,'FontAngle','italic');
