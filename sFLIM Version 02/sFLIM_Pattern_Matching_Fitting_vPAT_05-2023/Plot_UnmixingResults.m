clear variables;
clc;
i=1;
binning = 1; 
list_string = string([]);
list = {};

[filename, pathname] = uigetfile({'*results.mat','Results File (*results.mat)'}, 'Pick a Results file');
A = load([pathname filename]);

%%
num_pattern  = size(A.results.amp,3) - 1;
num_px = size(A.results.amp,1);
num_py = size(A.results.amp,2);
if  myIsField(A,'match')
   handles.match = A.results.match;
 else
   handles.match = 0;
end
if  myIsField(A,'figures')
   handles.figPA = A.results.figures;
 else
   handles.figPA = 1;
end
if  myIsField(A,'vParam')
   handles.vParam = A.results.vParam;
 else
   handles.vParam = [];
end
if  myIsField(A,'vParamDist')
   handles.vParamDist = A.results.vParamDist;
 else
   handles.vParamDist = [];
end

for i=1:num_pattern
  i_s=string(i);
  list_s=append("Pattern ", i_s);
  list_string(i)=string(list_s);
  i=i+1;
end
list_string=list_string';
list=cellstr(list_string);

if  myIsField(A,'list') %with older versions the Pattern List has not been stored
   list = A.results.list;
 else
 
 end

%%Plot Graphs
%%

nx      = num_px;
ny      = num_py;
num_pat = num_pattern;

scrsz = get(0,'ScreenSize');

fwidth = 800;
fheight = round(2*fwidth/3);
fx = 0.5*(scrsz(3)-fwidth);
fy = 0.5*(scrsz(4)-fheight);

py1 = 0.10;
px1 = 0.06;
wy  = 0.88;
wx  = 0.88;

% x = handles.head.ImgHdr.X0+(1:handles.par_binning:handles.head.ImgHdr.PixX)*handles.head.ImgHdr.PixelSize;
% y = handles.head.ImgHdr.Y0+(1:handles.par_binning:handles.head.ImgHdr.PixY)*handles.head.ImgHdr.PixelSize;

%res = reshape(res,[nx*ny num_pat+1]);
% lims = [0 max(max(A.results.amp))];
% if lims(2) <= lims(1)
%    lims(2) = lims(1)+1;
% end

%map  = [0 0 0; 0 0 1; 0 1 1; 0 1 0; 1 1 0; 1 0 0; 0.5 0 0];
%map1 = interp1(linspace(0, 1, 7), map, linspace(0, 1, 256));

map1 = gray;

%res = reshape(res,[nx ny num_pat+1]);

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
    imagesc(A.results.amp(:,:,n));
%    imagesc(x,y,res(:,:,n).*lims(2)./max(max(res(:,:,n))));

    set(R1,'DataAspectRatio', [1,1,1], ...
        'PlotBoxAspectRatio',[1 1 1], ...
        'XDir','normal', ...
        'YDir','reverse', ...
        'FontSize',9,...
        'CLim', [0 max(max(A.results.amp(:,:,n)))], ... % Change
        'color', [0 0 0]);

    xlabel('x  ','FontWeight','bold','FontSize',10,'FontAngle','italic');
    ylabel('y  ','FontWeight','bold','FontSize',10,'FontAngle','italic');
    m = n;

    colorbar('YColor',[0 0 0])
end

if handles.match == 1 || handles.match == 2 || handles.match == 3 || handles.match == 4
    if handles.match == 1
        name = sprintf('%s: lifetime fit / ns', cell2mat(list(k,:)));            
    elseif handles.match == 2 || handles.match ==3 || handles.match == 4
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
    intens = squeeze(squeeze(A.results.amp(:,:,1)));
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

    %image(x,y,im)
    image(im)
    colormap(spectrum)

    set(R1,'DataAspectRatio', [1,1,1], ...
        'PlotBoxAspectRatio',[1 1 1], ...
        'XDir','normal', ...
        'YDir','reverse', ...
        'FontSize',9,...
        'CLim', tmp_lims, ...
        'color', [0 0 0]);

    xlabel('x  ','FontWeight','bold','FontSize',10,'FontAngle','italic');
    ylabel('y  ','FontWeight','bold','FontSize',10,'FontAngle','italic');

    if handles.match == 1
        tst = 2 + round(((tmp_lims(1):.5:tmp_lims(2))-tmp_lims(1))./(tmp_lims(2)-tmp_lims(1)).*(size(spectrum,1)-3));
        n= 1;
        for ttst = tmp_lims(1):0.5:tmp_lims(2)
            stst(n) = cellstr(sprintf('%4.1f',ttst));
            n = n+1;
        end
        %disp(tmp_lims);
        colorbar;
        %colorbar('CLim', tmp_lims, 'YTickMode','manual','XColor',[0 0 0],'YColor',[0 0 0],'Box','off','YTick', tst, ...
            %'FontSize', 9,'TickDir','out', 'YTickLabel',stst);
          
    elseif handles.match == 2 || handles.match == 3 || handles.match == 4
        tst = 2 + round(((tmp_lims(1):5:tmp_lims(2))-tmp_lims(1))./(tmp_lims(2)-tmp_lims(1)).*(size(spectrum,1)-3));
        n= 1;
        for ttst = tmp_lims(1):5:tmp_lims(2)
            stst(n) = cellstr(sprintf('%4.1f',ttst));
            n = n+1;
        end
        %colorbar('CLim', tmp_lims, 'YTickMode','manual','XColor',[0 0 0],'YColor',[0 0 0],'Box','off','YTick', tst, ...
            %'FontSize', 9,'TickDir','out', 'YTickLabel',stst);
         colorbar;   
    end
end


if handles.match == 2 || handles.match == 3 || handles.match == 4

    Binding = 100.*A.results.amp(:,:,1)./(A.results.amp(:,:,1)+A.results.amp(:,:,2));

    %name = sprintf('Binding %d', cell2mat(list(k,:))); 
    name = sprintf('Binding');
    
    figure(handles.figPA+m+2);

    set(gcf,'Position',[fx fy fwidth fheight]);
    set(gcf,'Name',name,'NumberTitle','off');
    R1 = subplot('Position', [px1 py1 wx wy]);

    xl = [399 460 550 640 699];% 699];
    yl = [[0 0 1]; [0 1 1]; [0 1 0]; [1 1 0]; [1 0 0]];%; [1 0 1]];
    lambda   = 399:3:699;
    spectrum = interp1(xl, yl, lambda);

    %intens = squeeze(sqrt(sum(handles.tag,3)));
    intens = squeeze(A.results.amp(:,:,1)+A.results.amp(:,:,2));
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

    %image(x,y,im)
    image(im)

    colormap(spectrum)

    set(R1,'DataAspectRatio', [1,1,1], ...
        'PlotBoxAspectRatio',[1 1 1], ...
        'XDir','normal', ...
        'YDir','reverse', ...
        'FontSize',9,...
        'CLim', tmp_lims, ...
        'color', [0 0 0]);

    xlabel('x','FontWeight','bold','FontSize',10,'FontAngle','italic');
    ylabel('y','FontWeight','bold','FontSize',10,'FontAngle','italic');


    tst = 2 + round(((tmp_lims(1):10:tmp_lims(2))-tmp_lims(1))./(tmp_lims(2)-tmp_lims(1)).*(size(spectrum,1)-3));
    n= 1;
    for ttst = tmp_lims(1):10:tmp_lims(2)
        stst(n) = cellstr(sprintf('%4.1f',ttst));
        n = n+1;
    end
    %colorbar('CLim', tmp_lims, 'YTickMode','manual','XColor',[0 0 0],'YColor',[0 0 0],'Box','off','YTick', tst, ...
        %'FontSize', 9,'TickDir','out', 'YTickLabel',stst);
     colorbar;   
end



name = sprintf('Composite');

im = MakeComposite('data',A.results.amp(:,:,1:num_pat),'list',list);

figure(handles.figPA+num_pat+4);

set(gcf,'Position',[fx fy fwidth fheight]);
set(gcf,'Name',name,'NumberTitle','off');
R1 = subplot('Position', [px1 py1 wx wy]);

hold on;

cla

%image(x,y,reshape(im,[nx ny 3]))
image(num_px, num_py,reshape(im,[nx ny 3]))

set(R1,'DataAspectRatio', [1,1,1], ...
    'PlotBoxAspectRatio',[1 1 1], ...
    'XDir','normal', ...
    'YDir','reverse', ...
    'FontSize',9,...
    'color', [0 0 0]);

  %'CLim', lims,
  
xlabel('x  ','FontWeight','bold','FontSize',10,'FontAngle','italic');
ylabel('y  ','FontWeight','bold','FontSize',10,'FontAngle','italic');

hold off

s = {};
% tst = 1;
% for n = 1:handles.maxROI
%     if handles.ROIlist(n, handles.pulse)==1
%         s(tst) = {sprintf('ROI %d_%d',n,handles.pulse)};
%         tst = tst+1;
%     end
% end
% set(handles.listbox2,'Value',1);
% set(handles.listbox2,'String',s);
% listbox2_Callback(handles.listbox2, [], handles);

%plot_graphs



%%function myisfield

function isFieldResult = myIsField (inStruct, fieldName)
% inStruct is the name of the structure or an array of structures to search
% fieldName is the name of the field for which the function searches
isFieldResult = 0;
f = fieldnames(inStruct(1));
for i=1:length(f)
if(strcmp(f{i},strtrim(fieldName)))
isFieldResult = 1;
return;
elseif isstruct(inStruct(1).(f{i}))
isFieldResult = myIsField(inStruct(1).(f{i}), fieldName);
if isFieldResult
return;
end
end
end
end
