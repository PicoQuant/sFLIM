function varargout = KLF(varargin)
% KLF MATLAB code for KLF.fig
%      KLF, by itself, creates a new KLF or raises the existing
%      singleton*.
%
%      H = KLF returns the handle to a new KLF or the handle to
%      the existing singleton*.
%
%      KLF('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in KLF.M with the given input arguments.
%
%      KLF('Property','Value',...) creates a new KLF or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before KLF_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to KLF_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help KLF

% Last Modified by GUIDE v2.5 18-Mar-2013 18:28:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @KLF_OpeningFcn, ...
                   'gui_OutputFcn',  @KLF_OutputFcn, ...
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


% --- Executes just before KLF is made visible.
function KLF_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to KLF (see VARARGIN)

% Choose default command line output for KLF
handles.output = hObject;

handles.n_iter = 0;
handles.conv   = [];

handles.x     = [];
handles.llh   = [];

guidata(hObject, handles)

handles.exitflag = false;
% set(handles.axes1,'visible', 'on');
% set(handles.axes1, ...
%     'Box','off', ...
%     'XDir','normal', ...
%     'YDir','normal', ...
%     'XScale','linear', ...
%     'YScale','log', ...
%     'FontSize',9,...
%     'Color',[200 200 200]./255);
% ylabel(handles.axes1, 'convergence','FontWeight','bold','FontSize',10,'FontAngle','italic');
% xlabel(handles.axes1, 'iteration','FontWeight','bold','FontSize',10,'FontAngle','italic');

dontOpen = false;
nargin   = numel(varargin);
if nargin == 2
  handles.A = varargin{1};
  handles.Y = varargin{2};
end

if dontOpen
   disp('Improper input arguments. Pass a property value pair');
   disp('whose name is "changeme_main" and value is the handle');
   disp('to the changeme_main figure.');
end

% Update handles structure
guidata(hObject, handles);

KLF_1(hObject, handles);
end
end




% UIWAIT makes KLF wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = KLF_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;
varargout{1} = handles.x;
varargout{2} = handles.llh;
delete(hObject);
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.exitflag = true;
set(hObject,'Enable','off');
guidata(hObject, handles);
end


function [x, llh] = KLF_1(hObject, handles)

% Y : (n x m) matrix containing image data with 'n' pixels and 'm' bins
% A : (m x k) matrix containing component data with 'm' bins for 'k' components
% x : (n x k) matrik containing the amplitudes of the components for each pixel

A = handles.A; % Pattern Data
Y = handles.Y; % Image Data with bins

%%
K = A;           % Pattern Data
MNW = maxNumCompThreads; %MNW: maximum number of workers (cores) on the PC
%disp(MNW);
Num_workers = MNW;
numpix = (size(Y,1));
nx = sqrt(size(Y,1));
ny = nx;
X = [];
Conv = [];
Image_data = [];
LLH = [];
update_rate = 1.1;
N_iter = [];
Time = [];
Conv_size = [];
split = (nx*ny)/Num_workers;
split = floor(numpix/Num_workers);
rest = numpix-split * Num_workers;

spmd(Num_workers) % Single Program Multiple Data
    
    switch labindex
       
        case 1
            
        I=Y((labindex-1)*split+1:labindex*split,:);
        [x_out,llh_out,n_iter,conv,t] = Momentum_convergence(I,K,update_rate);
        
        case 2  
            
        I=Y((labindex-1)*split+1:labindex*split,:);
        [x_out,llh_out,n_iter,conv,t] = Momentum_convergence(I,K,update_rate);
        
        case 3  
            
        I=Y((labindex-1)*split+1:labindex*split,:);
        [x_out,llh_out,n_iter,conv,t] = Momentum_convergence(I,K,update_rate);
        
        case 4  
            
        I=Y((labindex-1)*split+1:labindex*split,:);
        [x_out,llh_out,n_iter,conv,t] = Momentum_convergence(I,K,update_rate);
        
        case 5  
            
        I=Y((labindex-1)*split+1:labindex*split,:);
        [x_out,llh_out,n_iter,conv,t] = Momentum_convergence(I,K,update_rate);
        
        case 6  
            
        I=Y((labindex-1)*split+1:labindex*split,:);
        [x_out,llh_out,n_iter,conv,t] = Momentum_convergence(I,K,update_rate);
        
        case 7  
            
        I=Y((labindex-1)*split+1:labindex*split,:);
        [x_out,llh_out,n_iter,conv,t] = Momentum_convergence(I,K,update_rate);
        
        case 8
            
        I=Y((labindex-1)*split+1:labindex*split,:);
        [x_out,llh_out,n_iter,conv,t] = Momentum_convergence(I,K,update_rate);
        
        case 9  
            
        I=Y((labindex-1)*split+1:labindex*split,:);
        [x_out,llh_out,n_iter,conv,t] = Momentum_convergence(I,K,update_rate);
        
        case 10  
            
        I=Y((labindex-1)*split+1:labindex*split,:);
        [x_out,llh_out,n_iter,conv,t] = Momentum_convergence(I,K,update_rate);
        
        case 11  
            
        I=Y((labindex-1)*split+1:labindex*split,:);
        [x_out,llh_out,n_iter,conv,t] = Momentum_convergence(I,K,update_rate);
        
        case 12  
            
        I=Y((labindex-1)*split+1:labindex*split,:);
        [x_out,llh_out,n_iter,conv,t] = Momentum_convergence(I,K,update_rate);
        
        case 13  
            
        I=Y((labindex-1)*split+1:labindex*split,:);
        [x_out,llh_out,n_iter,conv,t] = Momentum_convergence(I,K,update_rate);
        
        %case 14  
            
        %I=Y((labindex-1)*split+1:labindex*split,:);
        %[x_out,llh_out,n_iter,conv,t] = Momentum_convergence(I,K,update_rate);
        
        %case 15 
            
        %I=Y((labindex-1)*split+1:labindex*split,:);
        %[x_out,llh_out,n_iter,conv,t] = Momentum_convergence(I,K,update_rate);
        
        otherwise 
        I=Y((labindex-1)*split+1:labindex*split+rest,:);
        [x_out,llh_out,n_iter,conv,t] = Momentum_convergence(I,K,update_rate);
        
    end
   
end 

%save('FalseResults_workers.mat','X','Conv','Conv_size','Image_data','LLH','N_iter','Time','A');
 for j=1:Num_workers
    X =          [X ; x_out{j}];        %#ok<AGROW>
    Conv =       [Conv  conv{j}];      %#ok<AGROW>
    Conv_size(j)=size(conv{j},2);      %#ok<AGROW>
    Image_data = [Image_data ; I{j}];   %#ok<AGROW>
    LLH =        [LLH ; llh_out{j}];        %#ok<AGROW>
    N_iter =     [N_iter ; n_iter{j}];  %#ok<AGROW>
    Time =       [Time ; t{j}];         %#ok<AGROW>
 end
%save('Results_16workers.mat','X','Conv','Conv_size','Image_data','LLH','N_iter','Time','A');
handles.x=X ;
handles.llh = LLH;
delete(gcp('nocreate'));
msgbox('Linear Unmixing Done','sFLIM Results','warn');
guidata(hObject, handles);
end


