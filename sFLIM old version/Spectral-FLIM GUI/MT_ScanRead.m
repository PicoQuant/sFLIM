function [head, tag, tcspc, tcspc_x, tcspc_y, DC, AP] = MT_ScanRead(name, num_PIE, bg, timegate)

if (nargin<4)||isempty(timegate)
    timegate = [];
end

if (nargin<3)||isempty(bg)
    BG_AP = 0;
else
    BG_AP = bg;
end

if (nargin<2)||isempty(num_PIE)
    num_PIE = 1;
end

fname =  [name(1:end-4) '_FLIM.dat'];

if strcmp(name(end-2:end),'ht3')
    
    [head, im_tcspc, im_chan, im_line, im_col] = HT3_ScanRead(name);
    
elseif strcmp(name(end-2:end),'pt3')
    
    [head, im_tcspc, im_chan, im_line, im_col] = PT3_ScanRead(name);
    head.ImgHdr.PixX = head.ImgHdr.ScanWidthX;
    head.ImgHdr.PixY = head.ImgHdr.ScanWidthY;
    
elseif strcmp(name(end-2:end),'ptu')
    
    [head, im_tcspc, im_chan, im_line, im_col] = PTU_ScanRead(name);
    head.ImgHdr.PixX = head.ImgHdr_PixX;
    head.ImgHdr.PixY = head.ImgHdr_PixY;
    head.ImgHdr.X0   = head.ImgHdr_X0;
    head.ImgHdr.Y0   = head.ImgHdr_Y0;
    head.ImgHdr.PixelSize = head.ImgHdr_PixResol;
    head.ImgHdr.PixelTime = head.ImgHdr_PixelTime;
    head.Resolution  = 1e9*head.MeasDesc_Resolution;
    head.SyncRate    = 1./head.MeasDesc_GlobalResolution;
end


if ~isempty(head)
    
    nx = head.ImgHdr.PixX;
    ny = head.ImgHdr.PixY;
    
    maxch = max(im_chan);
    
    maxres     = max([head.Resolution]);
    Resolution = max([maxres 0.032]);
    chDiv      = Resolution/maxres;
    
    im_tcspc = ceil(im_tcspc./chDiv);
    Ngate    = double(max(im_tcspc));
    
    % to avoid late arriving markers in hydraharp
    %Ngate = double(1563);
    
    tcspc   = zeros(maxch, Ngate);
    tcspc_x = zeros(nx, maxch, Ngate);
    tcspc_y = zeros(ny, maxch, Ngate);
    
    for ch = 1:maxch
        tcspc(ch,:) = mHist(double(im_tcspc(im_chan == ch)),1:Ngate);
    end
    
    for y = 1:ny
        ind = (im_line == y);
        tmp1 = im_tcspc(ind);
        tmp2 = im_chan(ind);
        for ch = 1:maxch
            tcspc_y(y,ch,:) = mHist(double(tmp1(tmp2 == ch)),1:Ngate);
        end
    end
    
    for x = 1:nx
        ind = (im_col == x);
        tmp1 = im_tcspc(ind);
        tmp2 = im_chan(ind);
        for ch = 1:maxch
            tcspc_x(x,ch,:) = mHist(double(tmp1(tmp2 == ch)),1:Ngate);
        end
    end
    
    nch                = 1:maxch;
    chind              = sum(tcspc,2)<=100;
    nch(chind)         = [];
    tcspc(chind,:)     = [];
    tcspc_x(:,chind,:) = [];
    tcspc_y(:,chind,:) = [];
    
    bin              = 1:size(tcspc,2);
    %     ind              = sum(tcspc,1)<=10;
    %     bin(ind)         = [];
    %     tcspc(:,ind)     = [];
    %     tcspc_x(:,:,ind) = [];
    %     tcspc_y(:,:,ind) = [];
    
    tau   = ((1:size(tcspc,2))-0.5).*Resolution;
    nbin  = numel(tau);
    anzch = numel(nch);
    
    % Estimate background in measurement
    
    c_t   = [tcspc_x./ny; tcspc_y./nx];
    c_ts  = sum(c_t,3);
    nc    = size(c_ts,1);
    ind   = (1:floor(nc/100):nc);
    AP    = zeros(anzch,1);
    DC    = zeros(anzch,1);
    
    for a = 1:anzch
        [tmp, ord] = sort(c_ts(:,a),'descend');
        cnt = tmp(ind);
        tmp = squeeze(c_t(ord(ind),a,:));
        tmp = sort(tmp,2,'ascend');
        b_g = mean(tmp(:,2:floor(0.25.*nbin)),2);
        
        if BG_AP == 0
            V = [ones(size(cnt)) cnt];
            p  = lsqnonneg(V, b_g);
            DC(a) = p(1);
            AP(a) = p(2);
        else
            AP(a) = BG_AP;
            DC(a) = mean(b_g-cnt.*BG_AP);
        end
    end
    DC(DC<0) = 0;
    AP(AP<0) = 0;
    
    % Determine time-gates for PIE
    
    if isempty(timegate)
        [timegate, Ngate] = DetectTimeGates(tcspc', num_PIE, Resolution);
    else
        Ngate = 1 + timegate(1,2)+timegate(1,4)-timegate(1,1);
    end
    
    timegate(timegate>0) = bin(timegate(timegate>0));
    
    head.timegate = timegate;
    
    tcspc   = zeros(anzch, Ngate, num_PIE);
    tcspc_x = zeros(nx, anzch, Ngate, num_PIE);
    tcspc_y = zeros(ny, anzch, Ngate, num_PIE);
    tag     = zeros(ny, nx,anzch, num_PIE);
    
%     M_avail = memory;
%     
%     if M_avail.MaxPossibleArrayBytes > nx*ny*anzch*Ngate*num_PIE*8
%         
%         tim = zeros(ny, nx, anzch, Ngate, num_PIE);
%         
%         for a = 1:ny                     
%             ind = (im_line == a);
%             tmp1 = im_tcspc(ind);
%             tmp2 = im_chan(ind);
%             tmp3 = im_col(ind);
%             
%             for x = 1:nx
%                 ind = (tmp3 == x);
%                 tmp4 = double(tmp1(ind));
%                 tmp5 = tmp2(ind);
%                 
%                 n = 0;
%                 for pulse = 1:num_PIE
%                     for ch = 1:anzch
%                         n = n+1;
%                         if (timegate(n,3)==0)
%                             tim(a,x,ch,:,pulse) = mHist(tmp4(tmp5==nch(ch)),timegate(n,1):timegate(n,2));
%                         else
%                             tim(a,x,ch,:,pulse) = [mHist(tmp4(tmp5==nch(ch)),timegate(n,1):timegate(n,2)); mHist(tmp4(tmp5==nch(ch)),timegate(n,3):timegate(n,4))];
%                         end
%                     end
%                 end
%             end           
%         end
%         
%         tcspc   = shiftdim(sum(sum(tim,1),2),2);
%         tcspc_y = shiftdim(sum(tim,1),1);
%         tcspc_x = permute(sum(tim,2),[1 3 4 5 2]);
%         tag = permute(sum(tim,4),[1 2 3 5 4]);
%         
%         save([fname(1:end-4) '.mat'],'tim','-v7.3');
%         
%     else
        
        fid = fopen(fname,'w');
        
        for a = 1:ny
            
            tmp = zeros(nx, anzch, Ngate, num_PIE);
            
            ind = (im_line == a);
            tmp1 = im_tcspc(ind);
            tmp2 = im_chan(ind);
            tmp3 = im_col(ind);
            
            for x = 1:nx
                ind = (tmp3 == x);
                tmp4 = double(tmp1(ind));
                tmp5 = tmp2(ind);
                
                n = 0;
                for pulse = 1:num_PIE
                    for ch = 1:anzch
                        n = n+1;
                        if (timegate(n,3)==0)
                            tmp(x,ch,:,pulse) = mHist(tmp4(tmp5==nch(ch)),timegate(n,1):timegate(n,2));
                        else
                            tmp(x,ch,:,pulse) = [mHist(tmp4(tmp5==nch(ch)),timegate(n,1):timegate(n,2)); mHist(tmp4(tmp5==nch(ch)),timegate(n,3):timegate(n,4))];
                        end
                    end
                end
            end
            
            tcspc   = tcspc + shiftdim(sum(tmp,1),1);
            tcspc_x = tcspc_x + tmp;
            tcspc_y(a,:,:,:) = shiftdim(sum(tmp,1),1);
            tag(a,:,:,:) = permute(sum(tmp,3),[1 2 4 3]);
            
            fwrite(fid, uint16(tmp), 'uint16' );
        end
        
        fclose(fid);
        clear tmp im_chan im_col im_line im_tcspc c_t c_ts;
end