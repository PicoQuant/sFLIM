function [x] = KLF(pattern_data,image_data)

% Y : (n x m) matrix containing image data with 'n' pixels and 'm' bins
% A : (m x k) matrix containing component data with 'm' bins for 'k' components
% x : (n x k) matrik containing the amplitudes of the components for each pixel

A = pattern_data; % Pattern Data
Y = image_data; % Image Data with bins

x     = [];
llh   = [];

nx          = sqrt(size(Y,1));
ny          = nx;



[Num_workers, ~]  = getWorkersAvailable();
index_split       = floor(linspace(0,nx*ny,Num_workers+1))';
amplitude_out     = [];


parfor labindex  = 1:Num_workers
    
    I = Y; % temporary variables
    split_index = index_split
    [x_out,~,~,~,~] = Convergence(I(index_split(labindex)+1:split_index(labindex+1),:),A)
    amplitude_out =  [amplitude_out;x_out];    
end

x = amplitude_out;

delete(gcp('nocreate'));
msgbox('Linear Unmixing Done','sFLIM Results','warn');

end


