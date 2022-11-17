function MomentAnalysis(file)

tmp = strfind(file,'\');
pathname = file(1:tmp(end));
filename = file(tmp(end)+1:end);

load Get_Params;

[head, tag, tcspc, IRF, timname] = load_data([pathname filename]);

nx = floor(head.ImgHdr.PixX);
ny = floor(head.ImgHdr.PixY);
tx = floor(nx/binning);
ty = floor(ny/binning);

nch     = size(tcspc,1);
nbin    = size(tcspc,2);
num_PIE = size(tcspc,3);
tau     = head.tau;

load([pathname timname], 'stim'); 

if binning > 1
    tmp = stim(1:binning*tx,:,:,:,:);
    tmp = reshape(tmp, binning, tx, ny, nch, nbin, num_PIE);
    tmp = shiftdim(sum(tmp,1),1);
    tmp = tmp(:,1:binning*ty,:,:);
    tmp = reshape(tmp, tx, binning, ty, nch, nbin, num_PIE);
    stim = permute(sum(tmp,2), [1 3 4 5 6 2]);
end

tav = zeros(ty, tx, nch, num_PIE, 2);
    
if numel(stim)<2^26
    
    T   = permute(repmat(tau(:), [1 ty tx nch num_PIE]),[2 3 4 1 5]);
    
    H(:,:,1) = permute(sum(repmat(ones(nch,1)*tau   ,[1 1 num_PIE]).*IRF,2),[1 3 2])./permute(sum(IRF,2),[1 3 2]);
    H(:,:,2) = permute(sum(repmat(ones(nch,1)*tau.^2,[1 1 num_PIE]).*IRF,2),[1 3 2])./permute(sum(IRF,2),[1 3 2]);
    
    H = shiftdim(repmat(H,[1 1 1 ty tx]), 3);
    
    F(:,:,:,:,1) = permute(sum(T   .*stim,4),[1 2 3 5 4])./permute(sum(stim,4),[1 2 3 5 4]);
    F(:,:,:,:,2) = permute(sum(T.^2.*stim,4),[1 2 3 5 4])./permute(sum(stim,4),[1 2 3 5 4]);

    F(isnan(F)) = 0;

    tav(:,:,:,:,1) = F(:,:,:,:,1)     -  H(:,:,:,:,1);
    tav(:,:,:,:,2) = F(:,:,:,:,2)./2  - (H(:,:,:,:,2)/2  + H(:,:,:,:,1).*tav(:,:,:,:,1));
    
else
    
    T = permute(repmat(tau(:), [1 tx nch]),[4 2 3 1]);

    for P = 1:num_PIE
        H = zeros(nch,2);
        H(:,1) = permute(sum((ones(nch,1)*tau   ).*IRF(:,:,P),2),[1 3 2])./permute(sum(IRF(:,:,P),2),[1 3 2]);
        H(:,2) = permute(sum((ones(nch,1)*tau.^2).*IRF(:,:,P),2),[1 3 2])./permute(sum(IRF(:,:,P),2),[1 3 2]);

        H = shiftdim(repmat(H,[1 1 tx]), 2);
        if (size(H,2)==1)
            H = permute(H,[2 1 3]);
        end


        for yi = 1:ty
            F(:,:,1) = permute(sum(T   .*stim(yi,:,:,:,P),4),[1 2 3 5 4])./permute(sum(stim(yi,:,:,:,P),4),[1 2 3 5 4]);
            F(:,:,2) = permute(sum(T.^2.*stim(yi,:,:,:,P),4),[1 2 3 5 4])./permute(sum(stim(yi,:,:,:,P),4),[1 2 3 5 4]);
            
            F(isnan(F)) = 0;
            
            tav(yi,:,:,P,1) = F(:,:,1)     -  H(:,:,1);
            tav(yi,:,:,P,2) = F(:,:,2)./2  - (H(:,:,2)/2  + H(:,:,1).*squeeze(tav(yi,:,:,P,1)));
        end
    end

end

% tav(tav<0) = 0;
tav(:,:,:,:,2) = tav(:,:,:,:,2) - tav(:,:,:,:,1).^2;
tav(tav<0) = 0;
tav(:,:,:,:,2) = sqrt(tav(:,:,:,:,2));
clear stim;
save([pathname filename(1:end-4) '_Moments.mat'], 'tav');