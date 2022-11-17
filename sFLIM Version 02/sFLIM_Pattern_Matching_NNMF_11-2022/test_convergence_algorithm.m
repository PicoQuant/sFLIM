clear variables;
%%
load('14A_862Z_01_FRET_Slide_01_04_before_unmixing');
nx          = sqrt(size(Y,1));
ny          = nx;

%%
% save('14A_862Z_01_FRET_Slide_01_04_before_unmixing.mat','A','Y','-v7.3');
% gg;
% Three options for update rate + parallel processing using CPU
% 1. Start with a higher update_rate and then quench it after few epochs
% 2. Nesterov Accelarated momentum
% 3. Reserved for: A Smart-ass PhD Student finds another way


% Reserved for spmd use
% split = floor((nx*ny)/Num_workers);



[Num_workers, ~]  = getWorkersAvailable();
index_split       = floor(linspace(0,nx*ny,Num_workers+1))';


%%
tic;
amplitude_out     = [];
parfor labindex  = 1:Num_workers
    
    I               = Y; % temporary variables
    split_index     = index_split
    [x_out,~,~,~,~] = Momentum_convergence(I(index_split(labindex)+1:split_index(labindex+1),:),A)
    amplitude_out   = [amplitude_out;x_out];
    
end

x_convergenceNNMF = amplitude_out;

delete(gcp('nocreate'));
fprintf('\n First algorithm Linear Unmixing Done \n');
toc
%% 
tic;
amplitude_out = [];
parfor labindex  = 1:Num_workers
    
    I               = Y; % temporary variables
    split_index     = index_split
    [x_out,~,~,~,~] = Convergence(I(index_split(labindex)+1:split_index(labindex+1),:),A)
    amplitude_out   = [amplitude_out;x_out];
    
end

x_convergence = amplitude_out;

delete(gcp('nocreate'));
fprintf('Second algorithm Linear Unmixing Done');
toc

%% plot the differences


num_pat = size(x_convergenceNNMF,2);
map1    = gray;



res1 = reshape(x_convergenceNNMF,[nx ny num_pat]);
res2 = reshape(x_convergence,[nx ny num_pat]);


for n = 1:num_pat
    figure;
    colormap(map1);
    hold on;
    cla
    imshow(res2(:,:,n),[]);
%     imshow(res2(:,:,n) - res1(:,:,n),[min(min(res2(:,:,n) - res1(:,:,n)))*0 max(max(res2(:,:,n) - res1(:,:,n)))]);
    colorbar;
    hold off;
end

