function varargout = LeastSquares(varargin)
% LeastSquares MATLAB code for LeastSquares.fig
%      LeastSquares, by itself, creates a new LeastSquares or raises the existing
%      singleton*.
%
%      H = LeastSquares returns the handle to a new LeastSquares or the handle to
%      the existing singleton*.
%
%      LeastSquares('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LeastSquares.M with the given input arguments.
%
%      LeastSquares('Property','Value',...) creates a new LeastSquares or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LeastSquares_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LeastSquares_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LeastSquares

% Last Modified by GUIDE v2.5 16-Jan-2023 09:43:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LeastSquares_OpeningFcn, ...
                   'gui_OutputFcn',  @LeastSquares_OutputFcn, ...
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

%End initialization code - DO NOT EDIT


% --- Executes just before LeastSquares is made visible.
function LeastSquares_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LeastSquares (see VARARGIN)

% Choose default command line output for LeastSquares
handles.output = hObject;

handles.vParamDist = [];
handles.final   = [];

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

if nargin == 8
  handles.M = varargin{1};
  handles.Y = varargin{2};
  handles.numpix = varargin{3};
  handles.num_vpat = varargin{4};
  handles.num_pat = varargin{5};
  handles.vParam = varargin{6};
  handles.vParamDist_1 = varargin{7};
  handles.final_1 = varargin{8};
end

if dontOpen
   disp('Improper input arguments. Pass a property value pair');
   disp('whose name is "changeme_main" and value is the handle');
   disp('to the changeme_main figure.');
end

% Update handles structure
guidata(hObject, handles);

LeastSquares_1(hObject, handles);
end
end




% UIWAIT makes LeastSquares wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = LeastSquares_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
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


function [vParamDist,final] = LeastSquares_1(hObject, handles)

% Y : (n x m) matrix containing image data with 'n' pixels and 'm' bins
% A : (m x k) matrix containing component data with 'm' bins for 'k' components
% x : (n x k) matrik containing the amplitudes of the components for each pixel

M = handles.M; % Pattern Data
Y = handles.Y; % Image Data with bins
numpix = handles.numpix; % number of pixel
num_vpat = handles.num_vpat; % number of pattern in the pattern family for pattern fit
num_pat = handles.num_pat; % number of standard pattern
vParam = handles.vParam; % 
vParamDist_1 = handles.vParamDist_1; % results part 1
final_1 = handles.final_1; %results part 2

          
%%

MNW = maxNumCompThreads; %MNW: maximum number of workers (cores) on the PC
%disp(MNW);
Num_workers = MNW;

nx = sqrt (numpix); %only for quadratic images with same number of pixel in x and y 
ny = nx;

vParamDist_f = [];
final_f = [];

res = zeros(numpix, num_vpat, num_pat+1);
split = floor(numpix/Num_workers);
rest = numpix-split * Num_workers;


 spmd(Num_workers) % Single Program Multiple Data
   
  switch labindex 
     
    case 1
%         disp case1;
%         disp ((labindex-1)*split+1);
%         disp (labindex*split);
        for pix = (labindex-1)*split+1:labindex*split
            for n = 1:num_vpat

                res(pix,n,1:end-1)  = lsqnonneg(squeeze((M(n,:,:)))', Y(pix,:)');
                
                %res(pix,n,1:end-1)  = squeeze((M(n,:,:)))'\tmp(pix,:)';
                res(pix,n,end)      = sum(Y(pix,:)) - sum(res(pix,n,:));

                dist(pix,n) = res(pix,n,end);

            end
            [m,i] = min(abs(dist(pix,:)));
            vParamDist(pix) = handles.vParam(i);
            final(pix,:) = squeeze(res(pix,i,:));
          end

          

    case 2
%          disp case2;
%          disp ((labindex-1)*split+1);
%          disp (labindex*split);
         for pix = (labindex-1)*split+1:labindex*split
            for n = 1:num_vpat

                res(pix,n,1:end-1)  = lsqnonneg(squeeze((M(n,:,:)))', Y(pix,:)');
                
                %res(pix,n,1:end-1)  = squeeze((M(n,:,:)))'\tmp(pix,:)';
                res(pix,n,end)      = sum(Y(pix,:)) - sum(res(pix,n,:));

                dist(pix,n) = res(pix,n,end);

            end
            [m,i] = min(abs(dist(pix,:)));
            vParamDist(pix) = handles.vParam(i);
            final(pix,:) = squeeze(res(pix,i,:));
        end  
      
    case 3 
        for pix = (labindex-1)*split+1:labindex*split
            for n = 1:num_vpat

                res(pix,n,1:end-1)  = lsqnonneg(squeeze((M(n,:,:)))', Y(pix,:)');
                
                %res(pix,n,1:end-1)  = squeeze((M(n,:,:)))'\tmp(pix,:)';
                res(pix,n,end)      = sum(Y(pix,:)) - sum(res(pix,n,:));

                dist(pix,n) = res(pix,n,end);

            end
            [m,i] = min(abs(dist(pix,:)));
            vParamDist(pix) = handles.vParam(i);
            final(pix,:) = squeeze(res(pix,i,:));
        end 
    case 4
        for pix = (labindex-1)*split+1:labindex*split
            for n = 1:num_vpat

                res(pix,n,1:end-1)  = lsqnonneg(squeeze((M(n,:,:)))', Y(pix,:)');
                
                %res(pix,n,1:end-1)  = squeeze((M(n,:,:)))'\tmp(pix,:)';
                res(pix,n,end)      = sum(Y(pix,:)) - sum(res(pix,n,:));

                dist(pix,n) = res(pix,n,end);

            end
            [m,i] = min(abs(dist(pix,:)));
            vParamDist(pix) = handles.vParam(i);
            final(pix,:) = squeeze(res(pix,i,:));
        end 
    case 5 
        for pix = (labindex-1)*split+1:labindex*split
            for n = 1:num_vpat

                res(pix,n,1:end-1)  = lsqnonneg(squeeze((M(n,:,:)))', Y(pix,:)');
                
                %res(pix,n,1:end-1)  = squeeze((M(n,:,:)))'\tmp(pix,:)';
                res(pix,n,end)      = sum(Y(pix,:)) - sum(res(pix,n,:));

                dist(pix,n) = res(pix,n,end);

            end
            [m,i] = min(abs(dist(pix,:)));
            vParamDist(pix) = handles.vParam(i);
            final(pix,:) = squeeze(res(pix,i,:));
        end 
    case 6
        for pix = (labindex-1)*split+1:labindex*split
            for n = 1:num_vpat

                res(pix,n,1:end-1)  = lsqnonneg(squeeze((M(n,:,:)))', Y(pix,:)');
                
                %res(pix,n,1:end-1)  = squeeze((M(n,:,:)))'\tmp(pix,:)';
                res(pix,n,end)      = sum(Y(pix,:)) - sum(res(pix,n,:));

                dist(pix,n) = res(pix,n,end);

            end
            [m,i] = min(abs(dist(pix,:)));
            vParamDist(pix) = handles.vParam(i);
            final(pix,:) = squeeze(res(pix,i,:));
        end 
    case 7
        for pix = (labindex-1)*split+1:labindex*split
            for n = 1:num_vpat

                res(pix,n,1:end-1)  = lsqnonneg(squeeze((M(n,:,:)))', Y(pix,:)');
                
                %res(pix,n,1:end-1)  = squeeze((M(n,:,:)))'\tmp(pix,:)';
                res(pix,n,end)      = sum(Y(pix,:)) - sum(res(pix,n,:));

                dist(pix,n) = res(pix,n,end);

            end
            [m,i] = min(abs(dist(pix,:)));
            vParamDist(pix) = handles.vParam(i);
            final(pix,:) = squeeze(res(pix,i,:));
        end 
    case 8
        for pix = (labindex-1)*split+1:labindex*split
            for n = 1:num_vpat

                res(pix,n,1:end-1)  = lsqnonneg(squeeze((M(n,:,:)))', Y(pix,:)');
                
                %res(pix,n,1:end-1)  = squeeze((M(n,:,:)))'\tmp(pix,:)';
                res(pix,n,end)      = sum(Y(pix,:)) - sum(res(pix,n,:));

                dist(pix,n) = res(pix,n,end);

            end
            [m,i] = min(abs(dist(pix,:)));
            vParamDist(pix) = handles.vParam(i);
            final(pix,:) = squeeze(res(pix,i,:));
        end 
    case 9
        for pix = (labindex-1)*split+1:labindex*split
            for n = 1:num_vpat

                res(pix,n,1:end-1)  = lsqnonneg(squeeze((M(n,:,:)))', Y(pix,:)');
                
                %res(pix,n,1:end-1)  = squeeze((M(n,:,:)))'\tmp(pix,:)';
                res(pix,n,end)      = sum(Y(pix,:)) - sum(res(pix,n,:));

                dist(pix,n) = res(pix,n,end);

            end
            [m,i] = min(abs(dist(pix,:)));
            vParamDist(pix) = handles.vParam(i);
            final(pix,:) = squeeze(res(pix,i,:));
        end 
    case 10
        for pix = (labindex-1)*split+1:labindex*split
            for n = 1:num_vpat

                res(pix,n,1:end-1)  = lsqnonneg(squeeze((M(n,:,:)))', Y(pix,:)');
                
                %res(pix,n,1:end-1)  = squeeze((M(n,:,:)))'\tmp(pix,:)';
                res(pix,n,end)      = sum(Y(pix,:)) - sum(res(pix,n,:));

                dist(pix,n) = res(pix,n,end);

            end
            [m,i] = min(abs(dist(pix,:)));
            vParamDist(pix) = handles.vParam(i);
            final(pix,:) = squeeze(res(pix,i,:));
        end 
    case 11
        for pix = (labindex-1)*split+1:labindex*split
            for n = 1:num_vpat

                res(pix,n,1:end-1)  = lsqnonneg(squeeze((M(n,:,:)))', Y(pix,:)');
                
                %res(pix,n,1:end-1)  = squeeze((M(n,:,:)))'\tmp(pix,:)';
                res(pix,n,end)      = sum(Y(pix,:)) - sum(res(pix,n,:));

                dist(pix,n) = res(pix,n,end);

            end
            [m,i] = min(abs(dist(pix,:)));
            vParamDist(pix) = handles.vParam(i);
            final(pix,:) = squeeze(res(pix,i,:));
        end 
    case 12
        for pix = (labindex-1)*split+1:labindex*split
            for n = 1:num_vpat

                res(pix,n,1:end-1)  = lsqnonneg(squeeze((M(n,:,:)))', Y(pix,:)');
                
                %res(pix,n,1:end-1)  = squeeze((M(n,:,:)))'\tmp(pix,:)';
                res(pix,n,end)      = sum(Y(pix,:)) - sum(res(pix,n,:));

                dist(pix,n) = res(pix,n,end);

            end
            [m,i] = min(abs(dist(pix,:)));
            vParamDist(pix) = handles.vParam(i);
            final(pix,:) = squeeze(res(pix,i,:));
        end 
    case 13  
        for pix = (labindex-1)*split+1:labindex*split
            for n = 1:num_vpat

                res(pix,n,1:end-1)  = lsqnonneg(squeeze((M(n,:,:)))', Y(pix,:)');
                
                %res(pix,n,1:end-1)  = squeeze((M(n,:,:)))'\tmp(pix,:)';
                res(pix,n,end)      = sum(Y(pix,:)) - sum(res(pix,n,:));

                dist(pix,n) = res(pix,n,end);

            end
            [m,i] = min(abs(dist(pix,:)));
            vParamDist(pix) = handles.vParam(i);
            final(pix,:) = squeeze(res(pix,i,:));
        end 
        
    otherwise %the last worker works in addition on the 'rest' of pixel
        for pix = (labindex-1)*split+1:labindex*split+rest
            for n = 1:num_vpat

                res(pix,n,1:end-1)  = lsqnonneg(squeeze((M(n,:,:)))', Y(pix,:)');
                
                %res(pix,n,1:end-1)  = squeeze((M(n,:,:)))'\tmp(pix,:)';
                res(pix,n,end)      = sum(Y(pix,:)) - sum(res(pix,n,:));

                dist(pix,n) = res(pix,n,end);

            end
            [m,i] = min(abs(dist(pix,:)));
            vParamDist(pix) = handles.vParam(i);
            final(pix,:) = squeeze(res(pix,i,:));
        end 
  end    

end  
        

% put the results of all workers together
for j=1:Num_workers-1 %for all workers but the last

%   disp (j);
%   disp vor;
%   disp(size(vParamDist_f));
%   disp(size(vParamDist{j}.'));
%   disp(size(final_f));
%   disp(size(final{j}));
  vParamDist_2 = vParamDist{j}.';
  vParamDist_f = cat(1,vParamDist_f,vParamDist_2((j-1)*split+1:j*split,:));
  final_2 = final{j};
  final_f = cat(1,final_f,final_2((j-1)*split+1:j*split,:));
%   disp nach;
%   disp(size(vParamDist_f));
%   disp(size(vParamDist{j}.'));
%   disp(size(final_f));
%   disp(size(final{j}));
  
end
%for the rest of pixel with the last worker
  vParamDist_2 = vParamDist{Num_workers}.';
  vParamDist_f = cat(1,vParamDist_f,vParamDist_2((labindex-1)*split+1:labindex*split+rest,:));
  final_2 = final{Num_workers};
  final_f = cat(1,final_f,final_2((labindex-1)*split+1:labindex*split+rest,:));

% transmit the result to the Results_4.m file
setappdata(0,'vParamDist',vParamDist_f);
setappdata(0,'final',final_f);
% disp AusgabeLS;
% disp(size(vParamDist_f));
% disp(size(final_f));

delete(gcp('nocreate'));
msgbox('Linear Unmixing Done','sFLIM Results','warn');
guidata(hObject, handles);
end




