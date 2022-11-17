function varargout = MakeComposite(varargin)
% MAKECOMPOSITE M-file for MakeComposite.fig
%      MAKECOMPOSITE, by itself, creates a new MAKECOMPOSITE or raises the existing
%      singleton*.
%
%      H = MAKECOMPOSITE returns the handle to a new MAKECOMPOSITE or the handle to
%      the existing singleton*.
%
%      MAKECOMPOSITE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAKECOMPOSITE.M with the given input arguments.
%
%      MAKECOMPOSITE('Property','Value',...) creates a new MAKECOMPOSITE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MakeComposite_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MakeComposite_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MakeComposite

% Last Modified by GUIDE v2.5 10-Oct-2013 11:26:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MakeComposite_OpeningFcn, ...
                   'gui_OutputFcn',  @MakeComposite_OutputFcn, ...
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


% --- Executes just before MakeComposite is made visible.
function MakeComposite_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MakeComposite (see VARARGIN)

% Choose default command line output for MakeComposite

data = [];
list = [];
nargin = length(varargin);

if nargin < 2
    errordlg('No data given','Error')
else
    ind  = 1;
    while ind < nargin        
        if strcmpi(varargin{ind},'data') 
           data  = cell2mat(varargin(ind+1));
        end
        if strcmpi(varargin{ind},'list')  
           list  = varargin(ind+1);
           list  = list{:};
        end
        ind = ind + 2;
    end
end


if ~isempty(data)

    load Get_Params;
    
    handles.nx       = size(data,1);
    handles.ny       = size(data,2);
    handles.num_pat  = size(data,3);
    handles.Sat      = zeros(handles.num_pat,1);
    handles.selected = 1;
    
    data    = (reshape(data, [handles.nx*handles.ny handles.num_pat]));
    
    handles.data0 = data;
        
    tmp  = sort(reshape(data,[numel(data) 1]));
    dmax = round(tmp(floor(0.9995.*numel(tmp))));

    handles.data = data./(ones(size(data,1),handles.num_pat).*dmax);
    handles.data(handles.data>1) = 1;
    
    s = {};
    for ind = 1:handles.num_pat
        handles.ROIcolor(ind,:) = ROIcolor(mod(ind,size(ROIcolor,1)),:); %#ok<COLND>
        handles.spectrum(ind,:,:) = (0:0.002:1)'*handles.ROIcolor(ind,:);
        s(ind) = {sprintf('data %d',ind)};
    end
    if ~isempty(list)
        ind = min([handles.num_pat numel(list)]);  %        ind = min([handles.num_pat numel(list{:})]);
        s(1:ind) = list(1:ind);
    end
    set(handles.ROI_list,'Value',1);
    set(handles.ROI_list,'String',s);

    guidata(hObject, handles);
    update_plot(handles);
end

handles.SelectedColor     = [1,0,0];  % Currently selected color in the palette
handles.PalettePanel      = uibuttongroup('Parent',handles.Color_panel, ...
                                          'Units', 'Pixels',...
                                          'Position',[0 0 1 1],...
                                          'Title',{''},...
                                          'BorderType', 'none',...
                                          'SelectionChangeFcn', @(hObject, eventdata)MakeComposite('PalettePanelSelectionChanged',hObject, eventdata, guidata(hObject)));
                        
handles.SelectedColorText = uicontrol('Parent',handles.Color_panel,...
                                      'Units', 'Pixels',...
                                      'Style', 'text');
handles.MoreColorButton   = uicontrol('Parent',handles.Color_panel,...
                                      'Units', 'Pixels',...
                                      'String', 'More Colors ...',...
                                      'Callback',@(hObject, eventdata)MakeComposite('MoreColorButton_Callback',hObject, eventdata, guidata(hObject)));
guidata(hObject, handles);

% Dynamically create the color cells and palette tools and layout component
layoutComponent(handles);

% initalized the displayed color information
localUpdateColor(hObject, handles);

% Return user defined output if it is requested
% mOutputArgs{1} =@getImage; % (hObject)MakeComposite('getImage', guidata(hObject));
% 
% if nargout>0
%      [varargout{1:nargout}] = mOutputArgs{:};
% end
    
handles.output = [];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MakeComposite wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MakeComposite_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = update_plot(handles);
delete(hObject);


% --- Executes on button press in OK_button.
function OK_button_Callback(hObject, eventdata, handles)
% hObject    handle to OK_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume(handles.figure1);


% --- Executes on selection change in ROI_list.
function ROI_list_Callback(hObject, eventdata, handles)
% hObject    handle to ROI_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns ROI_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ROI_list

tmp = get(hObject, 'Value');
if numel(tmp)>0
   handles.selected = tmp(1);
   set(hObject, 'Value', tmp(1));
else
   handles.selected = 1;
   set(hObject, 'Value', 1);    
end
set(handles.slider1,'Value',handles.Sat(handles.selected));
handles.SelectedColor = handles.ROIcolor(handles.selected,:) ; %#ok<COLND>
set(handles.SelectedColorText, 'BackgroundColor', handles.SelectedColor);
set(get(handles.PalettePanel,'SelectedObject'),'Value',0);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ROI_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ROI_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
tmp            = get(hObject, 'Value');
handles.Sat(handles.selected) = tmp;
guidata(hObject, handles);
update_plot(handles);
    

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject,'Min',-4,'Max',4,'Value',0);

% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

sv = 1-get(hObject, 'Value')/100;
data = handles.data0;

tmp  = sort(reshape(data,[numel(data) 1]));
dmax = round(tmp(floor(sv.*numel(tmp))));

handles.data = data./(ones(size(data,1),handles.num_pat).*dmax);
handles.data(handles.data>1) = 1;
guidata(hObject, handles);
update_plot(handles);

% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject,'Min',0,'Max',0.5,'Value',0.05);


function im = update_plot(handles)

tmp = floor(500.*exp(repmat(handles.Sat',[handles.nx*handles.ny 1])).*handles.data)+1;
tmp(tmp>501) = 501;
tmp(tmp<1)   = 1;

ind = 1;
im  = handles.spectrum(ind, tmp(:,ind), :);
while ind < handles.num_pat
    ind = ind + 1;
    im  = im + handles.spectrum(ind, tmp(:,ind), :);
end

im(im>1)=1;

axes(handles.axes1)
hold on;
cla

image(1:handles.nx,1:handles.ny,reshape(im,[handles.nx handles.ny 3]))

set(handles.axes1,'DataAspectRatio', [1,1,1], ...
    'PlotBoxAspectRatio',[1 1 1], ...
    'XDir','normal', ...
    'YDir','reverse', ...
    'FontSize',9,...
    'color', [0 0 0]);

xlabel('x','FontWeight','bold','FontSize',10,'FontAngle','italic');
ylabel('y','FontWeight','bold','FontSize',10,'FontAngle','italic');

hold off


% --- Executes on button press in Reset_button.
function Reset_button_Callback(hObject, eventdata, handles)
% hObject    handle to Reset_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.slider1,'Value',0);
handles.Sat = zeros(handles.num_pat,1);
guidata(hObject, handles);
update_plot(handles);


function MoreColorButton_Callback(hObject, eventdata, handles)
% Callback called when the more color button is pressed.

color = handles.SelectedColor;
if isnan(color)
    color =[0 0 0];
end
color = uisetcolor(color);
if ~isequal(color, handles.SelectedColor)
    handles.SelectedColor =color;
    localUpdateColor(hObject, handles);
end


function PalettePanelSelectionChanged(hObject, eventdata, handles)
% Callback called when the selected color is changed in the colorPlatte

selected = get(hObject,'SelectedObject');
def = get(selected, 'UserData');
if ~isempty(def) && isfield(def,'Callback')
    def.Callback(selected, eventdata, handles);
end


function colorCell_Callback(hObject, eventdata, handles)
% Callback called when any color cell button is pressed

handles.SelectedColor = get(hObject, 'BackgroundColor');
localUpdateColor(hObject, handles)


function localUpdateColor(hObject, handles)
% function that updates the preview of the selected color
set(handles.SelectedColorText, 'BackgroundColor', handles.SelectedColor);
handles.ROIcolor(handles.selected,:) = handles.SelectedColor; %#ok<COLND>
handles.spectrum(handles.selected,:,:) = (0:0.002:1)'*handles.SelectedColor;
guidata(hObject, handles);
update_plot(handles);


function layoutComponent(handles)
% helper function that dynamically creats all the tools and color cells
% in the palette. It also positions all other UI objects properly. 
% get the definision of the layout
    
mColorEntries= {struct('Color','red',         'Callback',@colorCell_Callback),...
                struct('Color','green',       'Callback',@colorCell_Callback),...
                struct('Color','blue',        'Callback',@colorCell_Callback),...
                struct('Color','yellow',      'Callback',@colorCell_Callback),...
                struct('Color','magenta',     'Callback',@colorCell_Callback),...
                struct('Color','cyan',        'Callback',@colorCell_Callback),...
                struct('Color','white',       'Callback',@colorCell_Callback),...
                struct('Color',[0.5,0.5,0.5], 'Callback',@colorCell_Callback)};

mLayout.hgap = 5;
mLayout.vgap = 5;
mLayout.cellSize = 21;
mLayout.cellPerRow = 4;
mLayout.colorSampleSize = 21;
mLayout.moreColorButtonHeight = 30;

% calculate the preferred width and height
width  =  max([2*mLayout.colorSampleSize,(mLayout.cellSize+mLayout.hgap)*mLayout.cellPerRow]);
mLayout.cellRowNumber = ceil(length(mColorEntries)/ceil(width/(mLayout.cellSize+mLayout.vgap)));
height =  mLayout.cellRowNumber*(mLayout.cellSize+mLayout.vgap) ...
        + mLayout.colorSampleSize + mLayout.moreColorButtonHeight;
mLayout.preferredWidth = mLayout.hgap+width;
mLayout.preferredHeight = 2*mLayout.vgap+height;

% change the size of the color palette to the desired size, place
% the components, and then change size back.
setpixelposition(handles.PalettePanel, [0, 0, mLayout.preferredWidth, mLayout.preferredHeight]);

startY = mLayout.preferredHeight - mLayout.colorSampleSize;
startX = mLayout.hgap;

% place color sample
setpixelposition(handles.SelectedColorText, [startX, startY, ...
    mLayout.preferredWidth-2*mLayout.hgap, ...
    mLayout.colorSampleSize]);
% create color cells
for i=1:mLayout.cellRowNumber
    startY = startY - (mLayout.cellSize+mLayout.hgap);
    for j=1:mLayout.cellPerRow
        if ((i-1)*mLayout.cellPerRow + j)>length(mColorEntries)
            break;
        end
        color = mColorEntries{(i-1)*mLayout.cellPerRow + j};
        tooltip = mat2str(color.Color);
        if isfield(color,'Name')
            tooltip = color.Name;
        end
        control = uicontrol('Style','ToggleButton',...
            'TooltipString', tooltip,...
            'BackgroundColor',color.Color,...
            'Parent',handles.PalettePanel, ...
            'Units','Pixels',...
            'UserData',color,...
            'Position',[startX+(j-1)*(mLayout.cellSize+mLayout.hgap),...
            startY, ...
            mLayout.cellSize, mLayout.cellSize]);
    end
end

% place more color button
startY = startY - mLayout.moreColorButtonHeight - mLayout.vgap;
setpixelposition(handles.MoreColorButton,[startX, startY, ...
    mLayout.preferredWidth - 2*mLayout.hgap,mLayout.moreColorButtonHeight]);

% restore palette to the full size
set(handles.PalettePanel, 'units', 'normalized', 'Position', [0 0 1 1]);
