function [im] = Plot_UnmixingResults()
%%
% The function allows user to plot the results already availbale

load('List.mat','list');
[filename, pathname] = uigetfile({'*results.mat','Results File (*results.mat)'}, 'Pick a Results file');
A = load([pathname filename]);
%%
num_pattern  = size(A.results.amp,3) - 1;
num_px = size(A.results.amp,1);
num_py = size(A.results.amp,2);
%% Composite image with tools vary

im = MakeComposite('data',A.results.amp(:,:,1:num_pattern));
name = sprintf('Composite Image');
figure_handle = figure;
set(figure_handle,'Name',name,'NumberTitle','off');
hold on;
cla;

image(num_px , num_py,reshape(im,[num_px num_py 3]));
title(name);
set(gca,'DataAspectRatio', [1,1,1],'PlotBoxAspectRatio',[1 1 1],'XDir','normal','YDir','reverse');
xlabel('x (µm)','FontWeight','bold','FontSize',10,'FontAngle','italic');
ylabel('y (µm)','FontWeight','bold','FontSize',10,'FontAngle','italic');
hold off


end