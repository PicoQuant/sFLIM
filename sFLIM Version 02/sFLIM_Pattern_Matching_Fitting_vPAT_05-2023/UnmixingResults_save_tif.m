%%
clear variables;
clc;

load('List.mat');
[filename, pathname] = uigetfile({'*results.mat','Results File (*results.mat)'}, 'Pick a Results file');
A = load([pathname filename]);

filenametif1 = filename(1:end-4); 
% pathname
% filenametif
%% Size
num_pattern  = size(A.results.amp,3);
num_px = size(A.results.amp,1);
num_py = size(A.results.amp,2);
map1  = gray;
%% Save Results of different pattern to tif in the same folder as the result.mat file
%num_pattern = 3;
for n = 1:num_pattern 
    image = A.results.amp(:,:,n);
    imageU16 = uint16(image);
    filenametif = append(filenametif1,'_');
    filenametif = append(filenametif,num2str(n));
    filenametif = append(filenametif,'.tif');
    filenametif = append(pathname,filenametif);
  
    t = Tiff(filenametif,'w');
    tagstruct.ImageLength = num_px;
    tagstruct.ImageWidth = num_py;
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 16;
    tagstruct.SamplesPerPixel = 1;
    %tagstruct.RowsPerStrip = 16;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Software = 'MATLAB';
    %tagstruct % display tagstruct
    setTag(t,tagstruct);
    
    write(t,imageU16); %write the image as tiff on disk
%     figure;
%     colormap(map1);
%     cla
%     imshow(imageU16,[]);
%     colorbar;
end
close(t);

%% Composite Image with tools vary

% im = MakeComposite('data',A.results.amp(:,:,1:num_pattern),'list',list{1:num_pattern});
% name = sprintf('Composite Image');
% figure_handle = figure;
% set(figure_handle,'Name',name,'NumberTitle','off');
% hold on;
% cla;
% 
% image(num_px , num_py,reshape(im,[num_px num_py 3]));
% title(name);
% set(gca,'DataAspectRatio', [1,1,1],'PlotBoxAspectRatio',[1 1 1],'XDir','normal','YDir','reverse');
% xlabel('x (µm)','FontWeight','bold','FontSize',10,'FontAngle','italic');
% ylabel('y (µm)','FontWeight','bold','FontSize',10,'FontAngle','italic');
% hold off