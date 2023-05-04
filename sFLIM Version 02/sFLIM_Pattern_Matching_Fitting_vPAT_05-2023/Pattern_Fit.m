function varargout = Pattern_Fit(varargin)
% PATTERN_FIT M-file for Pattern_Fit.fig
%      PATTERN_FIT, by itself, creates a new PATTERN_FIT or raises the
%      existing
%      singleton*.
%
%      H = PATTERN_FIT returns the handle to a new PATTERN_FIT or the handle to
%      the existing singleton*.
%
%      PATTERN_FIT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PATTERN_FIT.M with the given input arguments.
%
%      PATTERN_FIT('Property','Value',...) creates a new PATTERN_FIT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Pattern_Fit_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Pattern_Fit_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Pattern_Fit

% Last Modified by GUIDE v2.5 20-Jan-2023 17:04:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Pattern_Fit_OpeningFcn, ...
                   'gui_OutputFcn',  @Pattern_Fit_OutputFcn, ...
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


% --- Executes just before Pattern_Fit is made visible.
function Pattern_Fit_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Pattern_Fit (see VARARGIN)

% Choose default command line output for Pattern_Fit
handles.output = hObject;


load Get_Params;
handles.par_lam_start = lam_start;
handles.par_lam_step  = lam_step;
handles.par_clam_min  = clam_min;
handles.par_clam_max  = clam_max;
handles.par_num_PIE   = num_PIE;
handles.ROIcolor      = ROIcolor;
handles.lamPIE        = lamPIE;

handles.ROI = 1;
handles.channel = 1;
handles.fit = 1;
handles.pulse = 1;
handles.backgr = 0.0;

%if length(varargin) < 2
%    uiwait(errordlg('No pattern given','Data missing'))
%elseif (length(varargin) == 2 && strcmpi(varargin{1},'file'))
%    handles.filename = cell2mat(varargin(3));
%    load(cell2mat(varargin(2)), 'pat', 'list', 'params');
%    handles.pat      = pat;
%    handles.List     = list;
%    handles.params   = params;
%end

load(cell2mat(varargin(2)), 'pat', 'list', 'params');

handles.pat      = pat;
handles.list     = list;
%handles.params   = params;
handles.filename = cell2mat(varargin(3));

tmp = strfind(handles.filename,'\');
handles.pathname = handles.filename(1:tmp(end));
handles.filename = handles.filename(tmp(end)+1:end);
%disp(handles.filename);
[handles.head, handles.tag, handles.tcspc, handles.IRF, handles.timname...
    , handles.AP, handles.DC] = load_data([handles.pathname handles.filename]);

handles.nch   = size(handles.tcspc,1);
handles.binw  = handles.head.tauw(:)./handles.head.tauw(1);
handles.taush = handles.head.tau;
handles.tau   = handles.head.tau(:)*ones(1,handles.par_num_PIE);
handles.tmp   = ones(size(handles.tau,1),1)*((0:handles.par_num_PIE-1).*handles.tau(end,1));
handles.tau   = reshape(handles.tau+handles.tmp, [numel(handles.tau) 1]);
handles.m0    = handles.head.m0;

% Update handles structure

guidata(hObject, handles);

set(handles.search,'Value',1);
set(handles.select_chan,'Value',1);
set(handles.in_spec, 'Value',1);


set(handles.popupmenu1, 'String',list);
set(handles.popupmenu6, 'String',{'Lifetime', 'FLIM FRET const Int. ratio', 'FLIM FRET const Ampl. ratio', 'FLIM FRET mono-exp.'})


tmp = {};
lam = (handles.par_lam_start+handles.par_lam_step.*(0:(handles.nch)));
for n = 1 : size(handles.tcspc,1)
    tmp(n) = {sprintf(' %s nm',num2str(lam(n)))};
end
set(handles.popupmenu2, 'String',tmp);

tmp = {};
PIE = 1:num_PIE;
for n = 1 : size(PIE,2)
    tmp(n) = {sprintf('Pulse %s',num2str(PIE(n)))};
end
set(handles.popupmenu3, 'String',tmp);


set(handles.edit1,'String', '0.5');
set(handles.edit2,'String', '3.1');
set(handles.edit3,'String', '0.28');
set(handles.edit4,'String', '1.8');
set(handles.edit29,'String', '0.0');

set(handles.wave1,'Enable', 'off')
set(handles.wave2,'Enable', 'off')
set(handles.wave3,'Enable', 'off')
set(handles.wave4,'Enable', 'off')
set(handles.pulse1,'Enable', 'off');
set(handles.pulse2,'Enable', 'off');
set(handles.pulse3,'Enable', 'off');
set(handles.pulse4,'Enable', 'off');

for n = 1:num_PIE
    eval(sprintf('set(handles.wave%d, \''Enable\'', \''on\'');', n))
    eval(sprintf('set(handles.pulse%d, \''Enable\'', \''on\'');', n))
    eval(sprintf('set(handles.wave%d, \''String\'', \''%d\'');', n, lamPIE(n))) 
end

set(handles.varycheck, 'Value', 0)
set(handles.vary1, 'Enable', 'off');
set(handles.vary2, 'Enable', 'off');
set(handles.no_of_pat, 'Enable', 'off');
set(handles.from, 'Enable', 'off');
set(handles.to, 'Enable', 'off');
set(handles.popupmenu6, 'Enable', 'off')

zoom on

axes(handles.axes1)

axis([0 floor(handles.taush(end)) -3.5 0.1]);
xlabel('time / ns','FontWeight','bold','FontSize',10,'FontAngle','italic');
ylabel('log (rel. frequency)','FontWeight','bold','FontSize',10,'FontAngle','italic');

axes(handles.axes2)

axis([0 floor(handles.taush(end)) 0 0.01]);
xlabel('time / ns','FontWeight','bold','FontSize',10,'FontAngle','italic');
ylabel('log (rel. frequency)','FontWeight','bold','FontSize',10,'FontAngle','italic');

axes(handles.axes3)

axis([handles.par_clam_min handles.par_clam_max 0 1.1]);
xlabel('wavelength of spectral mean / nm','FontWeight','bold','FontSize',10,'FontAngle','italic');
ylabel('rel. spectral irradiance','FontWeight','bold','FontSize',10,'FontAngle','italic');

axes(handles.axes4)

axis([0 floor(handles.taush(end)) -3.5 0.1]);
xlabel('time / ns','FontWeight','bold','FontSize',10,'FontAngle','italic');
ylabel('log (rel. frequency)','FontWeight','bold','FontSize',10,'FontAngle','italic');

axes(handles.axes5)

%axis([0 floor(handles.tau(end)) -3 0.1]);
xlabel('time / ns','FontWeight','bold','FontSize',10,'FontAngle','italic');
ylabel('log (rel. frequency)','FontWeight','bold','FontSize',10,'FontAngle','italic');

axes(handles.axes6)

xlabel('time / ns','FontWeight','bold','FontSize',10,'FontAngle','italic');



% UIWAIT makes Pattern_Fit wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Pattern_Fit_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in iPattern.
function iPattern_Callback(hObject, eventdata, handles)
% hObject    handle to iPattern (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% argin
selected_ROI = handles.ROI;
binw  = handles.binw;
tau   = handles.tau;
taush = handles.taush;
nch   = handles.nch;
par_num_PIE = handles.par_num_PIE;
pulse = handles.pulse;

% Reads the Pattern of the selected ROI and adds zeros left and right to
% the channel dimension.
pat = handles.pat(selected_ROI,:,:,:);
pat0 = zeros(size(pat,1),1,size(pat,3),size(pat,4));
pat = cat(2, pat0, pat, pat0);

% ???
decay = sum(pat,2)./permute(repmat(binw, [1 size(pat,1) 1 size(pat,4)]),[2 3 1 4]);
decay = squeeze(decay);

% 'decay'
% disp(decay');
% decay2 = decay';
% filen = 'newspeingelesen.mat';
% save('newspeingelesen.mat', 'decay2');


M = max(decay,[],1);
decay = decay./M(pulse);
%Untergrund abziehen
%decay = decay - min(decay)+0.00001;
%decay(decay<=1e-3) = 1e-5;

handles.decay = decay(:,pulse);


% ploting: overall ROI decay of the selected ROI
axes(handles.axes1)
hold on;
cla

plot(taush,log10(decay(:,pulse)),'.','LineWidth',2,'Color',handles.ROIcolor(selected_ROI,:));

hold off;

axis([0 floor(taush(end)) -3.5 0.1])

xlabel('time / ns','FontWeight','bold','FontSize',10,'FontAngle','italic');
ylabel('log (rel. frequency)','FontWeight','bold','FontSize',10,'FontAngle','italic');
% 'taush'
% disp(taush);



% ploting: IRF of each channel
axes(handles.axes2);
cla;
hold on;
for n = 1:nch
    plot(handles.taush, handles.IRF(n,:,pulse).*handles.head.tauw );
end
hold off;

axis auto

xlabel('time / ns','FontWeight','bold','FontSize',10,'FontAngle','italic');
ylabel('log (rel. frequency)','FontWeight','bold','FontSize',10,'FontAngle','italic');


axes(handles.axes3)

hold on;
cla
box on;

spec = sum(pat(1,:,:,pulse),3);
lam = (handles.par_lam_start+handles.par_lam_step.*(-1:nch))';
ts = handles.par_clam_min: 1: handles.par_clam_max;

[M,I] = max(spec);

handles.I = I;
spec = interp1(lam, spec./M, ts,'cubic');
handles.spec = [ts;spec]';
    
plot(ts,spec,'-','LineWidth',2,...
    'Color',handles.ROIcolor(selected_ROI,:));

hold off

axis([handles.par_clam_min handles.par_clam_max 0 1.1])

xlabel('wavelength of spectral mean / nm','FontWeight','bold','FontSize',10,'FontAngle','italic');
ylabel('rel. spectral irradiance','FontWeight','bold','FontSize',10,'FontAngle','italic');

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in iIRF.
function iIRF_Callback(hObject, eventdata, handles)
% hObject    handle to iIRF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% argin
nch   = size(handles.tcspc,1);
selected_ROI = handles.ROI;
pulse = handles.pulse;

% ploting: IRF of each channel
axes(handles.axes2);
cla;
hold on;
for n = 1:nch
    plot(handles.head.tau, handles.IRF(n,:,pulse).*handles.head.tauw );
end
hold off;

axis auto

xlabel('time / ns','FontWeight','bold','FontSize',10,'FontAngle','italic');
ylabel('log (rel. frequency)','FontWeight','bold','FontSize',10,'FontAngle','italic');




% --- Executes on button press in iSpectrum.
function iSpectrum_Callback(hObject, eventdata, handles)
% hObject    handle to iSpectrum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.in_spec,'Value') == 1
    
    set(handles.specname,'BackgroundColor', [0.83 0.82 0.78]);
    set(handles.specname,'String','');
    
   % argin
    selected_ROI = handles.ROI;
    binw  = handles.binw;
    tau   = handles.tau;
    taush = handles.taush;
    nch   = handles.nch;
    num_PIE = handles.par_num_PIE
    pulse = handles.pulse;


    % Reads the Pattern of the selected ROI and adds zeros left and right to
    % the channel dimension.
    pat = handles.pat(selected_ROI,:,:,:);
    pat0 = zeros(size(pat,1),1,size(pat,3),size(pat,4));
    pat = cat(2, pat0, pat, pat0);
    
    absorb = squeeze(sum(sum(pat,2),3));
    absorb = absorb./max(absorb);
    for i = 1:num_PIE
        eval(sprintf('set(handles.pulse%d,\''String\'',absorb(%d))',i,i));
    end
    
    axes(handles.axes3)

    hold on;
    cla
    box on;
    spec = sum(pat(1,:,:,pulse),3);
    lam = (handles.par_lam_start+handles.par_lam_step.*(-1:nch))';
    ts = handles.par_clam_min: 1: handles.par_clam_max;

    [M,I] = max(spec);

    handles.I = I;
    spec = interp1(lam, spec./M, ts,'cubic');
    handles.spec = [ts;spec]';

    plot(ts,spec,'-','LineWidth',2,...
        'Color',handles.ROIcolor(selected_ROI,:));

    hold off

    axis([handles.par_clam_min handles.par_clam_max 0 1.1])

    xlabel('wavelength of spectral mean / nm','FontWeight','bold','FontSize',10,'FontAngle','italic');
    ylabel('rel. spectral irradiance','FontWeight','bold','FontSize',10,'FontAngle','italic');

elseif get(handles.ex_spec,'Value') == 1

    [FileName,PathName] = uigetfile({'*.txt','txt-files';'*.*','All Files' });
    fid = fopen(FileName,'r');
    InputText = textscan(fid,'%f','delimiter','\n');
    
    set(handles.specname,'BackgroundColor', 'white');
    set(handles.specname,'String',FileName);
    
    spec = reshape(InputText{1},[2,length(InputText{1})/2])';
    handles.spec = [spec(:,1) spec(:,2)./max(spec(:,2))];
    fclose('all');

    axes(handles.axes3);
    hold on;
    cla;
    plot(spec(:,1),spec(:,2)./max(spec(:,2)))
    hold off 
end
    
guidata(hObject,handles);


% --- Executes on button press in fit.
function fit_Callback(hObject, eventdata, handles)
% hObject    handle to fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% argin
selected_ROI = handles.ROI;
selected_chan = handles.channel;
binw  = handles.binw;
tau   = handles.taush;
par_num_PIE = handles.par_num_PIE;
taulong = handles.tau;
pulse = handles.pulse;
backgr = 0.0;

% Prepare Data

pat = handles.pat(selected_ROI,:,:,:);

if get(handles.select_chan,'Value') == 1

    pat = pat./permute(repmat(binw,[1 size(pat,1) size(pat,2) size(pat,4)]),[2 3 1 4]);
    pat = pat./max(max(max(pat)));
    %pat(pat<=1e-3) = 1e-3;

    mesu = shiftdim(pat(:,selected_chan,:,pulse),1);
    irf = handles.IRF(selected_chan,:,pulse).*handles.binw';
    
elseif get(handles.select_chan,'Value') == 0

    decay = sum(pat,2)./permute(repmat(binw, [1 size(pat,1) 1 size(pat,4)]),[2 3 1 4]);
    decay = squeeze(decay);
    M = max(decay,[],1);
    decay = decay./M(pulse);
    %decay(decay<=1e-3) = 1e-3;
    
    mesu = decay(:,pulse)';
    
    binw = repmat(handles.binw',[size(handles.IRF,1) 1],1);
    irf = sum(handles.IRF(:,:,pulse).*...
        repmat(handles.binw',[size(handles.IRF,1) 1],1))./selected_chan;
    
end


handles.fit;

% Get Initial Values

backgr = mean(mesu(1:10));

if handles.fit == 1
    
    handles.fit = 0;
    
% generate initial values

backgr = mean(mesu(1:10));
params(5) = backgr;
min1 = min(mesu);
mesu_s = mesu - min1+0.00001;
mesu_s = mesu_s / max(mesu_s);
tmp = log10(mesu_s);
irf = irf/max(irf);
%tmp = (tmp - min(tmp));

[~,I1] = max(tmp);
I1 = I1+5;
% disp(M);
% disp(3*M/4);
% disp(I1);
tmp(1:I1) = [];
% disp(tmp);


for i = 1:length(tmp)
    if tmp(i) <= -0.2
        I2= I1 + i;
        tmp(1:i)=[];
        break
    end
end

for i = 1:length(tmp)
    if tmp(i) <= -2
        I3= I2 + i;
        tmp(1:i)=[];
        break
    end
end



for i = 1:length(tmp)
    if tmp(i) == min(tmp)
        I4 = I3 + i - 10;
        break
    end
end
%I1, I2, I3, I4
tau1 = (tau(I2)-tau(I1))/(log(mesu(I1))-log(mesu(I2)));
tau2 = (tau(I4)-tau(I3))/(log(mesu(I3))-log(mesu(I4)));
%     disp(tau1);
%     disp(tau2);
    A = max(mesu_s);
    A1 = mesu_s(I2);
    A2 = 1.5-mesu_s(I2);

    setparams = [A1,tau1,A2,tau2,backgr];
    
    A1 = str2double(get(handles.edit1,'String'));
    tau1 = str2double(get(handles.edit2,'String'));
    A2 = str2double(get(handles.edit3,'String'));
    tau2 = str2double(get(handles.edit4,'String'));
    backgr = str2double(get(handles.edit29,'String'));
    getparams = [A1,tau1,A2,tau2,backgr];
    
elseif handles.fit == 0
    
    A1 = str2double(get(handles.edit1,'String'));
    tau1 = str2double(get(handles.edit2,'String'));
    A2 = str2double(get(handles.edit3,'String'));
    tau2 = str2double(get(handles.edit4,'String'));
    backgr = str2double(get(handles.edit29,'String'));
    setparams = [A1,tau1,A2,tau2,backgr];
    getparams = setparams;
 
end

handles.backgr = backgr;

% Choice of fixed parameters

a1 = get(handles.checkbox2,'Value');
t1 = get(handles.checkbox3,'Value');
a2 = get(handles.checkbox4,'Value');
t2 = get(handles.checkbox5,'Value');
bgfix = get(handles.checkbox14,'Value');

fix = [a1,t1,a2,t2,bgfix];

antifix = ones(1,length(fix))-fix;  
  
setparams = antifix.*setparams;
getparams = fix.*getparams;

params = antifix.*setparams + fix.*getparams;

%Areas with very low values are not fitted
ind2=1;
for ind1=1:length(tau)
  if mesu(ind1)>1e-4
    tau2(ind2)=tau(ind1);
    mesu2(ind2)=mesu(ind1);
    irf2(ind2)=irf(ind1);
    ind2=ind2+1;
  end
end
tau=tau2; 
mesu=mesu2;
irf=irf2;

% Interpolation of tau to an aquidistant time array

t = tau;

res = diff([0 t]);

t = 0.0:min(res)/2:max(t);

ind = zeros(size(tau));
for k = 1:length(t)
    for i = 1:length(tau)
        if round(t(k)*1000)==round(tau(i)*1000)
         
            ind(1,i) = k;
           
        end
    end
end

%shift of irf curve one point to right
shift = 0;
irf=cat(2,zeros(1,shift),irf(1:end-shift));
  

% %%Untergrund abziehen
% mesu_sort = sort(mesu);
% backgr = mean(mesu_sort(1:10));%-std(mesu_sort(1:10));
% 
% mesu = mesu - backgr;
% mesu_size = size(mesu);
% 
% for index = 1:mesu_size(2)
%   
%   if ((mesu(index) < 1e-5))
%    mesu(index) = 1e-5;   
%   end 
%   %disp(mesu(index));
% end

%%Fit

f = @(setparams)norm(log10(Faltung(antifix.*setparams' + getparams,...
    irf, tau,t,ind))-log10(mesu));
% disp(f);  
  
%tic;
[setparams, Ssq, CNT] = LMFsolve(f,setparams,'XTol',1e-8,'MaxIter',1000);
%toc;
%CNT
params = antifix.*setparams' + fix.*getparams;

% Convolution

C = Faltung(params, irf, tau,t,ind);

% Plot and show results

axes(handles.axes4);
hold on
cla;

% plot(tau,(C),'r')
plot(tau,log10(C),'r')

% plot(tau,log10(C),'.','LineWidth',2,'Color',...
%     handles.ROIcolor(selected_ROI,:))

%plot(tau,log10(handles.IRF(selected_chan,:)),'b')
plot(tau,log10(mesu),'.','LineWidth',2,'Color',...
    handles.ROIcolor(selected_ROI,:))
  
% ploting IRF
mult = max(mesu)/max(irf);

plot(tau,log10(irf*mult));
%plot(tau,(irf*mult)); 
% IRF end

hold off
%axis([0 floor(tau(end)) 0 1])
axis([0 floor(tau(end)) -3.5 0.1])
xlabel('time / ns','FontWeight','bold','FontSize',10,'FontAngle','italic');
ylabel('log (rel. frequency)','FontWeight','bold','FontSize',10,'FontAngle','italic');

set(handles.edit1,'String',num2str(params(1)));
set(handles.edit2,'String',num2str(params(2)));
set(handles.edit3,'String',num2str(params(3)));
set(handles.edit4,'String',num2str(params(4)));
set(handles.edit29,'String',num2str(params(5)));
set(handles.text3,'String',num2str(params(1)/params(3)));
set(handles.text4,'String',num2str((params(1)*params(2))/(params(3)*params(4))));


% plot resiudals
axes(handles.axes6);

hold on
cla;
anf = 2;
abw = log10(C(anf:end))-log10(mesu(anf:end));
plot(tau(anf:end),abw);
dim = size(abw);
plot(tau(anf:end),zeros(dim));%zero line

hold off

%chi = Ssq;
chi = (norm(abw))^2;
set(handles.text23,'String',num2str(chi))

xlim([0 floor(tau(end))])
ylim auto
xlabel('time / ns','FontWeight','bold','FontSize',10,'FontAngle','italic');

guidata(hObject,handles);


% --- Executes on button press in plot_IV.
function plot_IV_Callback(hObject, eventdata, handles)
% hObject    handle to plot_IV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% argin
selected_ROI = handles.ROI;
selected_chan = handles.channel;
binw  = handles.binw;
tau   = handles.taush;
pulse = handles.pulse;


A1 = str2double(get(handles.edit1,'String'));
tau1 = str2double(get(handles.edit2,'String'));
A2 = str2double(get(handles.edit3,'String'));
tau2 = str2double(get(handles.edit4,'String'));
backgr = str2double(get(handles.edit29,'String'));

pat = handles.pat(selected_ROI,:,:,:);

if get(handles.select_chan,'Value') == 1

    pat = pat./permute(repmat(binw,[1 size(pat,1) size(pat,2) size(pat,4)]),[2 3 1 4]);
    pat = pat./max(max(max(pat)));
    %pat(pat<=1e-3) = 1e-3;

    mesu = shiftdim(pat(:,selected_chan,:,pulse),1);
    irf = handles.IRF(selected_chan,:,pulse).*handles.binw';

elseif get(handles.select_chan,'Value') == 0

    mesu = handles.decay;
    irf = sum(handles.IRF(:,:,pulse).*...
        repmat(handles.binw',[size(handles.IRF,1) 1],1))./selected_chan;

end

params = [A1,tau1,A2,tau2,backgr];

% Interpolation of tau to an aquidistant time array
ind = [];

t = tau;

res = diff([0 t]);

t = 0.0:min(res):max(t);

ind = zeros(size(tau));
for k = 1:length(t)
    for i = 1:length(tau)
        if round(t(k)*1000)==round(tau(i)*1000)
            ind(1,i) = k;
        end
    end
end

%ind = round(t.*1000) == round(tau.*1000);

C = Faltung(params, irf, tau,t,ind);

axes(handles.axes4);
hold on
cla;
%plot(tau,log10(handles.IRF(selected_chan,:)),'b')
plot(tau,log10(C),'r')
plot(tau,log10(mesu),'.','LineWidth',2,'Color',...
    handles.ROIcolor(selected_ROI,:))

hold off

axis([0 floor(tau(end)) -3 0.1])
xlabel('time / ns','FontWeight','bold','FontSize',10,'FontAngle','italic');
ylabel('log (rel. frequency)','FontWeight','bold','FontSize',10,'FontAngle','italic');

guidata(hObject,handles);


% --- Executes on button press in getIV - Get inital parameters
function getIV_Callback(hObject, eventdata, handles)
% hObject    handle to getIV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% argin
selected_ROI = handles.ROI;
selected_chan = handles.channel;
binw  = handles.binw;
tau   = handles.taush;
pulse = handles.pulse;

backgr = str2double(get(handles.edit29,'String'));

pat = handles.pat(selected_ROI,:,:,:);

if get(handles.select_chan,'Value') == 1

    pat = pat./permute(repmat(binw,[1 size(pat,1) size(pat,2) size(pat,4)]),[2 3 1 4]);
    pat = pat./max(max(max(pat)));
    %pat(pat<=1e-3) = 1e-3;

    mesu = shiftdim(pat(:,selected_chan,:,pulse),1);
    
elseif get(handles.select_chan,'Value') == 0

    mesu = handles.decay;
    
end

irf = handles.IRF(selected_chan,:,pulse).*handles.binw';

% generate initial values
% backgr = mean(mesu(1:10));
mesu_sort = sort(mesu);
backgr = mean(mesu_sort(1:10));
params(5) = backgr;
min1 = min(mesu);
% mesu = mesu - min1+0.00001;
mesu = mesu / max(mesu);
tmp = log10(mesu);
irf = irf/max(irf);
%tmp = (tmp - min(tmp));

[M,I1] = max(tmp);
I1 = I1+5;
% disp(M);
% disp(3*M/4);
% disp(I1);
tmp(1:I1) = [];
% disp(tmp);


for i = 1:length(tmp)
    if tmp(i) <= -0.2
        I2= I1 + i;
        tmp(1:i)=[];
        break
    end
end

for i = 1:length(tmp)
    if tmp(i) <= -2
        I3= I2 + i;
        tmp(1:i)=[];
        break
    end
end



for i = 1:length(tmp)
    if tmp(i) == min(tmp)
        I4 = I3 + i-10;
        break
    end
end
%I1, I2, I3, I4
tau1 = 0.8 *(tau(I2)-tau(I1))/(log(mesu(I1))-log(mesu(I2)));
tau2 = 0.8 *(tau(I4)-tau(I3))/(log(mesu(I3))-log(mesu(I4)));

A = max(mesu);

params = [(1.2*mesu(I2)),tau1,(1.3-mesu(I2)),tau2,params(5)];

% Interpolation of tau to an aquidistant time array

t = tau;

res = diff([0 t]);

t = 0.0:min(res):max(t);

ind = zeros(size(tau));
for k = 1:length(t)
    for i = 1:length(tau)
        if round(t(k)*1000)==round(tau(i)*1000)
            ind(1,i) = k;
        end
    end
end

% Convolution

C = Faltung(params, irf, tau,t,ind);

% Plot results

axes(handles.axes4);
hold on
cla;

plot(tau,log10(mesu(I1)*exp(-(tau-tau(I1))./tau1)))
plot(tau,log10(mesu(I3)*exp(-(tau-tau(I3))./tau2)))

plot(tau,log10(C),'r')
plot(tau,log10(mesu),'.','LineWidth',2,'Color',...
     handles.ROIcolor(selected_ROI,:))
% plot(tau(I1:I2),log10(mesu(I1:I2)),'-','LineWidth',2,'Color',...
%     'y')
% plot(tau(I3:I4),log10(mesu(I3:I4)),'-','LineWidth',2,'Color',...
%     'g')

hold off

axis([0 floor(tau(end)) -3.5 0.1])
xlabel('time / ns','FontWeight','bold','FontSize',10,'FontAngle','italic');
ylabel('log (rel. frequency)','FontWeight','bold','FontSize',10,'FontAngle','italic');

%params = [(mesu(I2)),tau1,(1-mesu(I2)),tau2,backgr];

set(handles.edit1,'String', num2str(params(1)))
set(handles.edit2,'String', num2str(params(2)))
set(handles.edit3,'String', num2str(params(3)))
set(handles.edit4,'String', num2str(params(4)))
set(handles.edit29,'String', num2str(params(5)))

guidata(hObject,handles);


% --- Executes on button press in copy_par.
function copy_par_Callback(hObject, eventdata, handles)
% hObject    handle to copy_par (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


set(handles.edit5,'String',get(handles.edit1,'String'))
set(handles.edit6,'String',get(handles.edit2,'String'))
set(handles.edit7,'String',get(handles.edit3,'String'))
set(handles.edit8,'String',get(handles.edit4,'String'))


% --- Executes on button press in save_spec.
function save_spec_Callback(hObject, eventdata, handles)
% hObject    handle to save_spec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

spec = handles.spec;

[FileName,PathName] = uiputfile({'*.txt','txt-files';...
        '*.*','All Files' },'Save spectrum',...
        'newspectrum.txt');
if FileName ~= 0
    %spec = interp1(lam, spec./M, ts,'cubic');
    save([PathName FileName], 'spec','-ascii');
end

guidata(hObject, handles);


% --- Executes on button press in create_pat.
function create_pat_Callback(hObject, eventdata, handles)
% hObject    handle to create_pat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%First part: reading the relative intensity from a imported emission
%spectrum for each channel

nch   = handles.nch;
spec  = handles.spec;

lam = (handles.par_lam_start+handles.par_lam_step.*(0:nch-1));

M = zeros(length(lam),1);
%chan = [];

for j = 1:length(lam)
    for i = 1:length(spec(:,1))
        if spec(i,1) == lam(j)
            if spec(i,2) > 1e-3
                M(j) = spec(i,2);
            end 
        end
    end
end

A1 = str2double(get(handles.edit5,'String'));
tau1 = str2double(get(handles.edit6,'String'));
A2 = str2double(get(handles.edit7,'String'));
tau2 = str2double(get(handles.edit8,'String'));
backgr = 0;
%backgr = str2double(get(handles.edit29,'String'));

% Second Part: creating the Pattern

% argin
%selected_chan = handles.channel;
binw  = handles.binw;
tau   = handles.taush;
taul  = handles.tau;
nch = handles.nch;
%selected_ROI = handles.ROI;
par_num_PIE = handles.par_num_PIE;

% vector for relative absorption coefficient

if par_num_PIE == 1
    
    pulse1 = repmat(str2double(get(handles.pulse1,'String')),[nch,length(tau)]);
    
    relInt = pulse1./max(pulse1,2);
    
    
elseif par_num_PIE == 2
    
    pulse1 = repmat(str2double(get(handles.pulse1,'String')),[nch,length(tau)]);
    pulse2 = repmat(str2double(get(handles.pulse2,'String')),[nch,length(tau)]);
    
    relInt = cat(2,pulse1,pulse2);
    
elseif par_num_PIE == 3
    
    pulse1 = repmat(str2double(get(handles.pulse1,'String')),[nch,length(tau)]);
    pulse2 = repmat(str2double(get(handles.pulse2,'String')),[nch,length(tau)]);
    pulse3 = repmat(str2double(get(handles.pulse3,'String')),[nch,length(tau)]);
    
    relInt = cat(2,pulse1,pulse2,pulse3);
    
elseif par_num_PIE == 4
    
    pulse1 = repmat(str2double(get(handles.pulse1,'String')),[nch,length(tau)]);
    pulse2 = repmat(str2double(get(handles.pulse2,'String')),[nch,length(tau)]);
    pulse3 = repmat(str2double(get(handles.pulse3,'String')),[nch,length(tau)]);
    pulse4 = repmat(str2double(get(handles.pulse4,'String')),[nch,length(tau)]);
    
    relInt = cat(2,pulse1,pulse2,pulse3,pulse4); 
end

relInt = relInt./max(max(relInt));

if any(isnan(relInt))
    uiwait(errordlg('Please insert absorbtion coefficients'))
end

% Interpolation of tau to an aquidistant time array
t = tau;
res = diff([0 t]);
t = 0.0:min(res):max(t);

ind = zeros(size(tau));
for k = 1:length(t)
    for i = 1:length(tau)
        if round(t(k)*1000)==round(tau(i)*1000)
            ind(1,i) = k;
        end
    end
end

% Pattern creation
pat = [];
for k = 1:par_num_PIE
    C2 = [];
    for i = 1:nch
        irf = handles.IRF(i,:,par_num_PIE).*binw';
        %shift irf by one channel to the right - better fitting
        %shift of irf curve one point to right
        shift = 0;
        irf=cat(2,zeros(1,shift),irf(1:end-shift));
        irf=irf/max(irf);
        params = [A1,tau1,A2,tau2,backgr];
%         A1,tau1,A2,tau2,backgr
        C1 = Pattern_Faltung(params,irf,tau,t,ind);
        C1 = M(i).*C1/sum(C1,2);
        C2 = cat(1,C2,C1);  
    end
    pat = cat(2,pat,C2);
end
%pat = pat/max(pat);
pat = pat.*relInt;

handles.C = pat;

% generation of a set of patterns varying one tau value from 'mint' to
% 'maxt' in 'no' steps

if get(handles.varycheck, 'Value') == 1
    if handles.match == 2 || handles.match == 3
    set(handles.vary2, 'Enable', 'off')
    set(handles.vary1, 'Value', 0)
    set(handles.vary1, 'Enable', 'off')
    end
    
    if isempty(get(handles.no_of_pat,'String'))
      uiwait(errordlg('Please enter no. of tau values and the ranges of tau values'))
    end  
    no = str2double(get(handles.no_of_pat,'String'));
    mint = str2double(get(handles.from,'String'));
    maxt = str2double(get(handles.to,'String'));
    step = (maxt-mint)/no;
    FRETAmp = repmat(A2, [1 no+1]);
    
    if handles.match == 1
        vPara = mint + (0:no).*step;
        if get(handles.vary1,'Value') == 1  
            tau1 = vPara;
            tau2 = repmat(tau2, [1, length(tau1)]);
        elseif get(handles.vary2,'Value') == 1
            tau2 = vPara;
            tau1 = repmat(tau1, [1, length(tau2)]);
        end
    elseif handles.match == 2 || handles.match == 3
        const = tau1*A1/(tau2*A2);
        
        vPara = mint + (0:no).*step;

        tau2 = tau2.*(1-(vPara./100));
        tau1 = repmat(tau1, [1, length(tau2)]);

        if handles.match == 2
            FRETAmp = tau1.*repmat(A1,[1 no+1])./(const*tau2)
        end
     elseif handles.match == 4
        vPara = mint + (0:no).*step;
        tau2 = tau2.*(1-(vPara./100));
        tau1 = repmat(1., [1., length(tau2)]);
        A1 = 0;
        FRETAmp = repmat(A2, [1., length(tau2)]);         
     end
        
    %tau1.*A1./(tau2.*FRETAmp);

    handles.patfam = cell(1,no+2);
    handles.patfam{end} = vPara;
    patplot = cell(1,no+1);

    for n = 1:no+1
        pat = [];
        for k = 1:par_num_PIE
            C2 = [];
            for i = 1:nch
                irf = handles.IRF(i,:,par_num_PIE).*binw';
                params = [A1,tau1(n),FRETAmp(n), tau2(n),backgr];
                C1 = Pattern_Faltung(params,irf,tau,t,ind);
                C1 = M(i).*C1/sum(C1,2);
                C2 = cat(1,C2,C1);
            end
            pat = cat(2,pat,C2);
        end

        pat = pat.*relInt;
        
        %the first pattern of the pattern family 
        %is stored as single pattern which can opened in pattern matching
        %GUI
        if n == 1 
          handles.C = pat;
        end  
%         'pat'
%         disp(pat/max(pat));
%         'tau'
%         disp(tau);
        patplot{n} = pat; %for display
    
        pat = pat.*permute(repmat(handles.binw,[handles.par_num_PIE size(pat,1)]),[2 1]);
        pat = reshape(pat,[handles.nch,length(handles.taush),handles.par_num_PIE]);
    
        handles.patfam{n} = pat;  %for save file     
    end

    axes(handles.axes5);
    hold on
    cla;
    %normalization, determin maximum
    for i = 1:no
         pat = patplot{i};
         [~,I] = max(max(pat,[],2));
         maxpat = max(pat(I,:));        
    end
         maxpatfam = max(maxpat);
   %prapare pattern for plotting     
    for i = 1:no
         pat = patplot{i}/maxpatfam; %normalization
         [~,I] = max(max(pat,[],2));
%          disp(pat);
         plot(taul,log10(pat(I,:)),'r');  
    end

    hold off

    axis([0 floor(taul(end)) -3.5 0.1])
    xlabel('time / ns','FontWeight','bold','FontSize',10,'FontAngle','italic');
    ylabel('log (rel. frequency)','FontWeight','bold','FontSize',10,'FontAngle','italic');
%     disp(taul);
else
    axes(handles.axes5);
    hold on
    cla;

    [~,I] = max(max(pat,[],2));
    pat = pat/max(pat); %normalization
    plot(taul,log10(pat(I,:)),'r');  

    hold off

    axis([0 floor(taul(end)) -3.5 0.1])
    xlabel('time / ns','FontWeight','bold','FontSize',10,'FontAngle','italic');
    ylabel('log (rel. frequency)','FontWeight','bold','FontSize',10,'FontAngle','italic');
    
end


guidata(hObject,handles)



% --- Executes on button press in save_pat.
function save_pat_Callback(hObject, eventdata, handles)
% hObject    handle to save_pat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pat = handles.C;

if get(handles.varycheck, 'Value') == 0
[FileName,PathName] = uiputfile({'*.pat','Pattern';...
        '*.*','All Files'},'Save pattern',...
        'newpattern');
else 
[FileName,PathName] = uiputfile({'*.vpat','variable Pattern';...
        '*.*','All Files'},'Save pattern',...
        'newpattern');     
end

if FileName ~= 0
    pat = pat.*permute(repmat(handles.binw,[handles.par_num_PIE size(pat,1)]),[2 1]);
    pat = reshape(pat,[handles.nch,length(handles.taush),handles.par_num_PIE]);
    if get(handles.varycheck, 'Value') == 1
        pat = {pat, handles.patfam, handles.match};
    end
    save([PathName FileName], 'pat');
end

% patfam1 = (handles.patfam{2});
% 'patfam'
% disp((patfam1)/max(patfam1));
% disp((patfam1)/1.0);
% handles.match
guidata(hObject, handles);


function pulse1_Callback(hObject, eventdata, handles)
% hObject    handle to pulse1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pulse1 as text
%        str2double(get(hObject,'String')) returns contents of pulse1 as a double


% --- Executes during object creation, after setting all properties.
function pulse1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pulse1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pulse2_Callback(hObject, eventdata, handles)
% hObject    handle to pulse2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pulse2 as text
%        str2double(get(hObject,'String')) returns contents of pulse2 as a double


% --- Executes during object creation, after setting all properties.
function pulse2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pulse2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pulse3_Callback(hObject, eventdata, handles)
% hObject    handle to pulse3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pulse3 as text
%        str2double(get(hObject,'String')) returns contents of pulse3 as a double


% --- Executes during object creation, after setting all properties.
function pulse3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pulse3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pulse4_Callback(hObject, eventdata, handles)
% hObject    handle to pulse4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pulse4 as text
%        str2double(get(hObject,'String')) returns contents of pulse4 as a double


% --- Executes during object creation, after setting all properties.
function pulse4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pulse4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.search, 'Value') == 1
    handles.fit = 1;
end

tmp = get(hObject,'Value');
if tmp ~= handles.ROI
    handles.ROI = tmp;
end;

guidata(hObject, handles);

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.search, 'Value') == 1
    handles.fit = 1;
end

tmp = get(hObject,'Value');
if tmp ~= handles.channel
    handles.channel = tmp;
end;

guidata(hObject, handles);

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


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


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tmp = get(hObject,'Value');
if tmp ~= handles.pulse
    handles.pulse = tmp;
end;

guidata(hObject, handles);

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu6.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.match = get(hObject, 'Value');

tmp = get(hObject,'Value');
if tmp == 1
    set(handles.text47, 'String','no. of tau values')
    set(handles.text49, 'String', 'ns     to   ')
    set(handles.text50, 'String', 'ns')
    set(handles.vary2, 'Enable', 'on')
    set(handles.vary1, 'Enable', 'on')
    set(handles.edit5, 'Enable', 'on')
    set(handles.edit6, 'Enable', 'on')
    set(handles.edit7, 'Enable', 'on')
else
    set(handles.text47, 'String','no. of FRET Eff. values')
    set(handles.text49, 'String', '%      to   ')
    set(handles.text50, 'String', '%')
    set(handles.from, 'String', '30')
    set(handles.to, 'String', '80')
    set(handles.vary2, 'Value', 1)
    set(handles.vary2, 'Enable', 'off')
    set(handles.vary1, 'Value', 0)
    set(handles.vary1, 'Enable', 'off')
    if tmp == 4
      set(handles.edit5, 'Enable', 'off')
      set(handles.edit6, 'Enable', 'off')
      set(handles.edit7, 'Enable', 'off')
    else
      set(handles.edit5, 'Enable', 'on')
      set(handles.edit6, 'Enable', 'on')
      set(handles.edit7, 'Enable', 'on')
    end

end

guidata(hObject, handles);

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu6 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu6


% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in search.
function search_Callback(hObject, eventdata, handles)
% hObject    handle to search (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value') == 1
    handles.fit = 1;
elseif get(hObject,'Value') == 0
    handles.fit = 0;
end

guidata(hObject, handles);

% Hint: get(hObject,'Value') returns toggle state of search


% --- Executes on button press in select_chan.
function select_chan_Callback(hObject, eventdata, handles)
% hObject    handle to select_chan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.fit = 1;

if get(hObject,'Value') == 1
    set(handles.popupmenu2,'Enable', 'on')
elseif get(hObject,'Value') == 0
    set(handles.popupmenu2,'Enable','off')
end

guidata(hObject,handles);
% Hint: get(hObject,'Value') returns toggle state of select_chan


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.search, 'Value') == 1
    handles.fit = 1;
end

guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


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

if get(handles.search, 'Value') == 1
    handles.fit = 1;
end

guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


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

if get(handles.search, 'Value') == 1
    handles.fit = 1;
end

guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


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


function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.search, 'Value') == 1
    handles.fit = 1;
end

guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
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

% --- Executes on button press in varycheck.
function varycheck_Callback(hObject, eventdata, handles)
% hObject    handle to varycheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject, 'Value') == 1
    set(handles.vary1, 'Enable', 'on');
    set(handles.vary2, 'Enable', 'on');
    set(handles.no_of_pat, 'Enable', 'on');
    set(handles.from, 'Enable', 'on');
    set(handles.to, 'Enable', 'on');
    set(handles.popupmenu6, 'Enable', 'on')
    set(handles.edit5, 'Enable', 'on')
    set(handles.edit6, 'Enable', 'on')
    set(handles.edit7, 'Enable', 'on')
    handles.match = get(handles.popupmenu6, 'Value');
elseif get (hObject, 'Value') == 0
    set(handles.vary1, 'Enable', 'off');
    set(handles.vary2, 'Enable', 'off');
    set(handles.no_of_pat, 'Enable', 'off');
    set(handles.from, 'Enable', 'off');
    set(handles.to, 'Enable', 'off');
    set(handles.edit5, 'Enable', 'on')
    set(handles.edit6, 'Enable', 'on')
    set(handles.edit7, 'Enable', 'on')
    set(handles.popupmenu6, 'Enable', 'off')
end

if handles.match == 2 || handles.match == 3 || handles.match == 4
    set(handles.vary2, 'Enable', 'off')
    set(handles.vary1, 'Value', 0)
    set(handles.vary1, 'Enable', 'off')
end
if handles.match == 4
    set(handles.edit5, 'Enable', 'off')
    set(handles.edit6, 'Enable', 'off')
    set(handles.edit7, 'Enable', 'off')
    
end

guidata(hObject, handles);


% Hint: get(hObject,'Value') returns toggle state of varycheck


% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox7


% --- Executes on button press in checkbox8.
function checkbox8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox8


% --- Executes on button press in checkbox9.
function checkbox9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox9


% --- Executes on button press in checkbox10.
function checkbox10_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox10


% --- Executes on button press in in_spec.
function in_spec_Callback(hObject, eventdata, handles)
% hObject    handle to in_spec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of in_spec


% --- Executes on button press in ex_spec.
function ex_spec_Callback(hObject, eventdata, handles)
% hObject    handle to ex_spec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ex_spec

% --- Executes on button press in in_spec.
function vary1_Callback(hObject, eventdata, handles)
% hObject    handle to in_spec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of in_spec


% --- Executes on button press in ex_spec.
function vary2_Callback(hObject, eventdata, handles)
% hObject    handle to ex_spec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ex_spec


function wave1_Callback(hObject, eventdata, handles)
% hObject    handle to wave1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wave1 as text
%        str2double(get(hObject,'String')) returns contents of wave1 as a double


% --- Executes during object creation, after setting all properties.
function wave1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wave1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wave2_Callback(hObject, eventdata, handles)
% hObject    handle to wave2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wave2 as text
%        str2double(get(hObject,'String')) returns contents of wave2 as a double


% --- Executes during object creation, after setting all properties.
function wave2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wave2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wave3_Callback(hObject, eventdata, handles)
% hObject    handle to wave3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wave3 as text
%        str2double(get(hObject,'String')) returns contents of wave3 as a double


% --- Executes during object creation, after setting all properties.
function wave3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wave3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wave4_Callback(hObject, eventdata, handles)
% hObject    handle to wave4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wave4 as text
%        str2double(get(hObject,'String')) returns contents of wave4 as a double


% --- Executes during object creation, after setting all properties.
function wave4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wave4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function no_of_pat_Callback(hObject, eventdata, handles)
% hObject    handle to no_of_pat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(hObject, 'String', num2str(round(str2num(get(hObject, 'String')))));

% Hints: get(hObject,'String') returns contents of no_of_pat as text
%        str2double(get(hObject,'String')) returns contents of no_of_pat as a double


% --- Executes during object creation, after setting all properties.
function no_of_pat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to no_of_pat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function from_Callback(hObject, eventdata, handles)
% hObject    handle to from (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of from as text
%        str2double(get(hObject,'String')) returns contents of from as a double


% --- Executes during object creation, after setting all properties.
function from_CreateFcn(hObject, eventdata, handles)
% hObject    handle to from (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function to_Callback(hObject, eventdata, handles)
% hObject    handle to to (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of to as text
%        str2double(get(hObject,'String')) returns contents of to as a double


% --- Executes during object creation, after setting all properties.
function to_CreateFcn(hObject, eventdata, handles)
% hObject    handle to to (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




%% Convolution of IRF and a double exponential decay

function [C] = Faltung(params, irf, tau, t, ind)

A1 = params(1);
tau1 = params(2);
if tau1 == 0
   tau1 = 1; 
end
A2 = params(3);
tau2 = params(4);
if tau2 == 0
    tau2 = 1;
end
backgr = params(5);
%disp(backgr);
% The IRF is interpolated and normalized

tau(end+1) = 0;
tau = circshift(tau,1);
irf(end+1) = min(irf);
irf = circshift(irf,1);

irf = interp1(tau,irf,t);
irf = irf./sum(irf);

% Computes a exponential decay function.

decay = A1*exp(-t./tau1) + A2*exp(-t./tau2) + backgr;

%max(decay)

% Calculates the convolution of irf and decay adding an offset

C = conv(irf, decay);

%max(C)
tmp = [];
for i = ind
    tmp(end+1) = C(i);
end
C = tmp;

%C = C*ind;
%C(0) = [];

%C(C<=1e-3) = 1e-3; 
C(C<=1e-4) = params(5); % einziger Unterschied zu Pattern_Faltung
C(C<=1e-5) = 1e-5; % if params(5) = 0

%% Convolution of IRF and a double exponential decay

function [C] = Pattern_Faltung(params, irf, tau, t, ind)

A1 = params(1);
tau1 = params(2);
if tau1 == 0
   tau1 = 1; 
end
A2 = params(3);
tau2 = params(4);
if tau2 == 0
    tau2 = 1;
end

backgr = params(5);
% The IRF is interpolated and normalized

tau(end+1) = 0;
tau = circshift(tau,1);
irf(end+1) = min(irf);
irf = circshift(irf,1);

irf = interp1(tau,irf,t);
irf = irf./sum(irf);

% Computes a exponential decay function.

decay = A1*exp(-t./tau1) + A2*exp(-t./tau2) + backgr;
%max(decay)

% Calculates the convolution of irf and decay adding an offset

C = conv(irf, decay);
%max(C)

tmp = [];
for i = ind
    tmp(end+1) = C(i);
end
C = tmp;
C(C<=params(5)) = params(5); %set background
C(C<=1e-5) = 1e-5; % if background is 0.0 then limit 


% --- Executes on button press in checkbox14.
function checkbox14_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox14



function edit29_Callback(hObject, eventdata, handles)
% hObject    handle to edit29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit29 as text
%        str2double(get(hObject,'String')) returns contents of edit29 as a double


% --- Executes during object creation, after setting all properties.
function edit29_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
