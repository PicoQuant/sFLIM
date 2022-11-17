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

%
% save('R:\Projekte\882_MT200_sFLIM_PhD_Rohilla\WP\08_Publications\02_Publication 2 - HuLu + Cross-linking\ROH_Mac\Cross-labelling_and_Opal_Multiplex\14A_930Z_02_FRET_Slide_ROI_02\14A_930Z_02_FRET_Slide_ROI_02.mat','A','Y','-v7.3');
% keyboard;
% Two options for update rate + parallel processing using CPU
% 1. Start with a higher update_rate and then quench it after few epochs
% 2. Nesterov Accelarated momentum



% Reserved for spmd use
% split = floor((nx*ny)/Num_workers);

[Num_workers, ~]  = getWorkersAvailable();
index_split       = floor(linspace(0,nx*ny,Num_workers+1))';
amplitude_out     = [];


parfor labindex  = 1:Num_workers
    
    I = Y; % temporary variables
    split_index = index_split
    [x_out,~,~,~,~] = Momentum_convergence(I(index_split(labindex)+1:split_index(labindex+1),:),A);
    %     [x_out,~,~,~,~] = Convergence(I(index_split(labindex)+1:split_index(labindex+1),:),A)
    %[x_out,~,~,~,~] = Momentum_convergence_2(I(index_split(labindex)+1:split_index(labindex+1),:),A);
    amplitude_out =  [amplitude_out;x_out];
    
    %     x_out = I(index_split(labindex)+1:split_index(labindex+1),:);
    %     amplitude_out =  [amplitude_out;x_out];
    
end

x = amplitude_out;

delete(gcp('nocreate'));
msgbox('Linear Unmixing Done','sFLIM Results','warn');

%save('FalseResults_workers.mat','X','Conv','Conv_size','Image_data','LLH','N_iter','Time','A');
%  for j=1:Num_workers
%     X =          [X ; x_out{j}];        %#ok<AGROW>
%     Conv =       [Conv  conv{j}];      %#ok<AGROW>
%     Conv_size(j)=size(conv{j},2);      %#ok<AGROW>
%     Image_data = [Image_data ; I{j}];   %#ok<AGROW>
%     LLH =        [LLH ; llh_out{j}];        %#ok<AGROW>
%     N_iter =     [N_iter ; n_iter{j}];  %#ok<AGROW>
%     Time =       [Time ; t{j}];         %#ok<AGROW>
%  end


end


