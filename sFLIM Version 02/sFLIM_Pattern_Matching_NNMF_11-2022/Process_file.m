function Process_file(name, IRFname)

tmp = strfind(name,'\');
pathname = name(1:tmp(end));

names = dir(name);

for k = 1:numel(names)
    
    fname = names(k).name;
    load('Get_Params.mat');
    
    IRF      = [];
    timegate = [];
    
    if (nargin>1)&&(~isempty(IRFname))
        [IRF, timegate] = IRF_Read(IRFname, num_PIE);
    end
    
    [head, tag, tcspc, tcspc_x, tcspc_y, DC, AP] = MT_ScanRead([pathname fname], num_PIE, BG_AP, timegate);
    
%     save([pathname fname(1:end-4) '_MT_ScanRead_output.mat'], 'tcspc','tag','tcspc_x','tcspc_y','DC','AP','-v7.3');
 
    if (BG_corr == 0) % ignore bg
        
        DC = zeros(size(DC));
        AP = zeros(size(AP));
        
    elseif (BG_corr == 1) % take bg data from specified file
        
        tmp = dir([BG_file(1:end-4) '_DATA.mat']);
        
        if (numel(tmp)==0)
            [dummy, dummy, dummy, dummy, dummy, DC, AP] = MT_ScanRead(BG_file, num_PIE, BG_AP, timegate);
        else
            load([BG_file(1:end-4) '_DATA.mat'],'DC','AP');
        end
        
    end
    
    if ~isempty(head)
        if isempty(IRF)
            IRF = Calc_mIRF(head, tcspc);
        end
        
        maxres = max([head.Resolution]);
        Resolution = max([maxres 0.032]);
        
        nx   = head.ImgHdr.PixX;
        ny   = head.ImgHdr.PixY;
        nch  = size(tcspc,1);
        nbin = min([size(tcspc,2) size(IRF,2)]);
        nPIE = size(tcspc,3);
        
        mp = ceil(1.05./Resolution);      % Goal is to shift this maximum to 1.0 ns
        
        numch = zeros(nch, nPIE);
        
        for ch = 1:nch
            tmp1 = shiftdim(IRF(ch, :, :));
            tmp2 = shiftdim(tcspc(ch, :, :));
            
            for PIE = 1:nPIE
                [m0,m0] = min(abs(cumsum(tmp1(:,PIE)) - sum(tmp1(:,PIE))/2));
                numch(ch, PIE) = mp-m0;
                if numch(ch, PIE) > 0        % we need to add channels to the front
                    if PIE == 1
                        IRF(ch,:,PIE) = [tmp1(end-numch(ch, PIE)+1:end,end  ); tmp1(1:end-numch(ch, PIE),PIE)];
                        tcspc(ch,:,PIE) = [tmp2(end-numch(ch, PIE)+1:end,end  ); tmp2(1:end-numch(ch, PIE),PIE)];
                    else
                        IRF(ch,:,PIE) = [tmp1(end-numch(ch, PIE)+1:end,PIE-1); tmp1(1:end-numch(ch, PIE),PIE)];
                        tcspc(ch,:,PIE) = [tmp2(end-numch(ch, PIE)+1:end,PIE-1); tmp2(1:end-numch(ch, PIE),PIE)];
                    end
                else                         % we need to add channels to the back
                    if PIE == nPIE
                        IRF(ch,:,PIE) = [tmp1(abs(numch(ch, PIE))+1:end,PIE); tmp1(1:abs(numch(ch, PIE)),1)];
                        tcspc(ch,:,PIE) = [tmp2(abs(numch(ch, PIE))+1:end,PIE); tmp2(1:abs(numch(ch, PIE)),1)];
                    else
                        IRF(ch,:,PIE) = [tmp1(abs(numch(ch, PIE))+1:end,PIE); tmp1(1:abs(numch(ch, PIE)),PIE+1)];
                        tcspc(ch,:,PIE) = [tmp2(abs(numch(ch, PIE))+1:end,PIE); tmp2(1:abs(numch(ch, PIE)),PIE+1)];
                    end
                end
                
                for rx = 1:nx
                    tmp3 = shiftdim(tcspc_x(rx,ch,:,:),2);
                    if numch(ch, PIE) > 0        % we need to add channels to the front
                        if PIE == 1
                            tcspc_x(rx,ch,:,PIE) = [tmp3(end-numch(ch, PIE)+1:end,end  ); tmp3(1:end-numch(ch, PIE),PIE)];
                        else
                            tcspc_x(rx,ch,:,PIE) = [tmp3(end-numch(ch, PIE)+1:end,PIE-1); tmp3(1:end-numch(ch, PIE),PIE)];
                        end
                    else                         % we need to add channels to the back
                        if PIE == nPIE
                            tcspc_x(rx,ch,:,PIE) = [tmp3(abs(numch(ch, PIE))+1:end,PIE); tmp3(1:abs(numch(ch, PIE)),1)];
                        else
                            tcspc_x(rx,ch,:,PIE) = [tmp3(abs(numch(ch, PIE))+1:end,PIE); tmp3(1:abs(numch(ch, PIE)),PIE+1)];
                        end
                    end
                end
                for ry = 1:ny
                    tmp3 = shiftdim(tcspc_y(ry,ch,:,:),2);
                    if numch(ch, PIE) > 0        % we need to add channels to the front
                        if PIE == 1
                            tcspc_y(ry,ch,:,PIE) = [tmp3(end-numch(ch, PIE)+1:end,end  ); tmp3(1:end-numch(ch, PIE),PIE)];
                        else
                            tcspc_y(ry,ch,:,PIE) = [tmp3(end-numch(ch, PIE)+1:end,PIE-1); tmp3(1:end-numch(ch, PIE),PIE)];
                        end
                    else                         % we need to add channels to the back
                        if PIE == nPIE
                            tcspc_y(ry,ch,:,PIE) = [tmp3(abs(numch(ch, PIE))+1:end,PIE); tmp3(1:abs(numch(ch, PIE)),1)];
                        else
                            tcspc_y(ry,ch,:,PIE) = [tmp3(abs(numch(ch, PIE))+1:end,PIE); tmp3(1:abs(numch(ch, PIE)),PIE+1)];
                        end
                    end
                end
                
            end
        end
        
        %         for pulse = 1:num_PIE
        %             for ch = 1:nch
        %                 IRF(ch,:,pulse) = IRF(ch,:,pulse) - mean(IRF(ch,1:(mp-10),pulse),2);
        %             end
        %         end
        IRF(IRF<0) = 0;
        IRF = IRF./repmat(sum(IRF,2),[1 size(IRF,2) 1]);
        KRF = IRF;
        
        c_t = shiftdim(sum(permute(tcspc,[2 1 3]),1),1);    % tcspc : [nch nbin num_PIE]
        
        aAP = 1+AP;
        
        BG = c_t - (c_t - repmat(DC(:).*nx.*ny, [1 num_PIE]))./repmat(aAP(:), [1 num_PIE]);
        BG = permute(repmat(BG, [1 1 nbin]), [1 3 2]);
        
        tcspc           = tcspc - BG;
        tcspc(tcspc<=0) = 1e-5;
        
        % rebin the time axis to compress data
        
        t0 = 0.032;
        tbin = round(t0/Resolution);
        if tbin<1
            tbin = 1;
            t0   = Resolution;
        end
        
        m0 = mp;
        m1 = 8+ceil(m0/tbin);
        t1 = t0;
        t2 = (1:m1).*t1;
        
        M0 = permute(IRF, [2 1 3]);
        M1 = M0;
        
        if tbin>1
            ne = floor(nbin/tbin)*tbin;
            M1 = reshape(M0(1:ne,:,:), tbin, ne/tbin, nch, nPIE);
            M1 = shiftdim(mean(M1),1);
        end
        M0(1:m1,:,:,:) = M1(1:m1,:,:,:);
        m2 = m1;
        
        while size(M1,1)>m2+1
            M1 = M1(m2+1:end,:,:,:);
            ne = floor(size(M1,1)/2)*2;
            M1 = reshape(M1(1:ne,:,:,:), 2, ne/2, nch, nPIE);
            M1 = shiftdim(mean(M1,1),1);
            
            m2 = min([8 size(M1,1)]);
            M0(m1+(1:m2),:,:,:) = M1(1:m2,:,:,:);
            m1  = m1+m2;
            t1  = t1*2;
            t2 = [t2 t2(end)+(1:m2)*t1];
        end
        
        IRF = permute(M0(1:m1,:,:), [2 1 3]);
        
        clear M0 M1
        
        M0 = permute(tcspc, [2 1 3]);
        M1 = M0;
        
        m1  = 8+ceil(m0/tbin);
        t1  = t0;
        t2 = (1:m1).*t1;
        
        if tbin>1
            ne = floor(nbin/tbin)*tbin;
            M1 = reshape(M0(1:ne,:,:), tbin, ne/tbin, nch, nPIE);
            M1 = shiftdim(mean(M1),1);
        end
        M0(1:m1,:,:,:) = M1(1:m1,:,:,:);
        m2 = m1;
        
        while size(M1,1)>m2+1
            M1 = M1(m2+1:end,:,:,:);
            ne = floor(size(M1,1)/2)*2;
            M1 = reshape(M1(1:ne,:,:,:), 2, ne/2, nch, nPIE);
            M1 = shiftdim(mean(M1,1),1);
            
            m2 = min([8 size(M1,1)]);
            M0(m1+(1:m2),:,:,:) = M1(1:m2,:,:,:);
            m1  = m1+m2;
            t1  = t1*2;
            t2 = [t2 t2(end)+(1:m2)*t1];
        end
        
        tcspc = permute(M0(1:m1,:,:), [2 1 3]);
        clear M0 M1
        
        head.tauw = diff([0 t2]);
        head.tau  = t2 - diff([0 t2])/2;
        head.m0   = m0;
        
        save([pathname fname(1:end-4) '_DATA.mat'], 'head', 'tag', 'tcspc', 'IRF', 'KRF', 'AP', 'DC');        
                
        if isempty(dir([pathname fname(1:end-4) '_FLIM.dat']))
            
           load([pathname fname(1:end-4) '_FLIM.mat'],'tim') 
           
           c_t = permute(sum(tim,4),[1 2 3 5 4]);    % c_t : [ny nx nch num_PIE]
           
           BG = c_t - (c_t - permute(repmat(DC(:), [1 ny nx num_PIE]),[2 3 1 4]))./permute(repmat(aAP(:), [1 ny nx num_PIE]),[2 3 1 4]);
           BG = permute(repmat(BG, [1 1 1 1 nbin]), [1 2 3 5 4]);
           
           tim        = tim - BG;
           
           clear BG;
           
           tim(tim<0) = 0;
                      
           tim = permute(tim, [4 1 2 3 5]);   % M0 : nbin ny nx nch num_PIE
           
           for ch = 1:nch
               tmp = squeeze(tim(:, :, :, ch, :));
               for PIE = 1:nPIE
                   if numch(ch, PIE) > 0        % we need to add channels to the front
                       if PIE == 1
                           tim(:,:,:,ch,PIE) = [tmp(end-numch(ch, PIE)+1:end,:,:,end  ); tmp(1:end-numch(ch, PIE),:,:,PIE)];
                       else
                           tim(:,:,:,ch,PIE) = [tmp(end-numch(ch, PIE)+1:end,:,:,PIE-1); tmp(1:end-numch(ch, PIE),:,:,PIE)];
                       end
                   else                         % we need to add channels to the back
                       if PIE == nPIE
                           tim(:,:,:,ch,PIE) = [tmp(abs(numch(ch, PIE))+1:end,:,:,PIE); tmp(1:abs(numch(ch, PIE)),:,:,1)];
                       else
                           tim(:,:,:,ch,PIE) = [tmp(abs(numch(ch, PIE))+1:end,:,:,PIE); tmp(1:abs(numch(ch, PIE)),:,:,PIE+1)];
                       end
                   end
               end
           end
           
           m1  = 8+ceil(m0/tbin);
           
           M1 = tim;
           if tbin>1
               ne = floor(nbin/tbin)*tbin;
               M1 = reshape(tim(1:ne,:,:,:), tbin, ne/tbin, nch, nPIE);
               M1 = shiftdim(mean(M1));
           end
           tim(1:m1,:,:,:,:) = M1(1:m1,:,:,:,:);
           m2 = m1;
           
           while size(M1,1)>m2+1
               M1 = M1(m2+1:end,:,:,:,:);
               ne = floor(size(M1,1)/2)*2;
               M1 = reshape(M1(1:ne,:,:,:,:), 2, ne/2, ny, nx, nch, num_PIE);
               M1 = shiftdim(sum(M1),1);
               
               m2 = min([8 size(M1,1)]);
               tim(m1+(1:m2),:,:,:,:) = M1(1:m2,:,:,:,:);
               m1  = m1+m2;
           end
           
           stim = permute(tim(1:m1,:,:,:), [2 3 4 1 5]); %  stim : ny nx nch nbin num_PIE
           
           clear tim M1
           
        else
            num  = nx*nch*nbin*num_PIE;
            
            fid1 = fopen([pathname fname(1:end-4) '_FLIM.dat'],'r');
            
            stim = zeros(ny, nx, nch, size(tcspc,2), num_PIE);
            
            for a = 1:ny
                tim = fread(fid1, num, 'uint16');
                tim = reshape(tim, [nx nch nbin num_PIE]);
                
                c_t = shiftdim(sum(permute(tim,[3 1 2 4]),1));    % c_t : [nx nch num_PIE]
                
                BG = c_t - (c_t - permute(repmat(DC(:), [1 nx num_PIE]),[2 1 3]))./permute(repmat(aAP(:), [1 nx num_PIE]),[2 1 3]);
                BG = permute(repmat(BG, [1 1 1 nbin]), [1 2 4 3]);
                
                tim        = tim - BG;
                tim(tim<0) = 0;
                
                M0 = permute(tim, [3 1 2 4]);   % M0 : nbin nx nch num_PIE
                
                for ch = 1:nch
                    tmp = squeeze(M0(:, :, ch, :));
                    for PIE = 1:nPIE
                        if numch(ch, PIE) > 0        % we need to add channels to the front
                            if PIE == 1
                                M0(:,:,ch,PIE) = [tmp(end-numch(ch, PIE)+1:end,:,end  ); tmp(1:end-numch(ch, PIE),:,PIE)];
                            else
                                M0(:,:,ch,PIE) = [tmp(end-numch(ch, PIE)+1:end,:,PIE-1); tmp(1:end-numch(ch, PIE),:,PIE)];
                            end
                        else                         % we need to add channels to the back
                            if PIE == nPIE
                                M0(:,:,ch,PIE) = [tmp(abs(numch(ch, PIE))+1:end,:,PIE); tmp(1:abs(numch(ch, PIE)),:,1)];
                            else
                                M0(:,:,ch,PIE) = [tmp(abs(numch(ch, PIE))+1:end,:,PIE); tmp(1:abs(numch(ch, PIE)),:,PIE+1)];
                            end
                        end
                    end
                end
                
                m1  = 8+ceil(m0/tbin);
                
                M1 = M0;
                if tbin>1
                    ne = floor(nbin/tbin)*tbin;
                    M1 = reshape(M0(1:ne,:,:), tbin, ne/tbin, nch, nPIE);
                    M1 = shiftdim(mean(M1));
                end
                M0(1:m1,:,:,:) = M1(1:m1,:,:,:);
                m2 = m1;
                
                while size(M1,1)>m2+1
                    M1 = M1(m2+1:end,:,:,:);
                    ne = floor(size(M1,1)/2)*2;
                    M1 = reshape(M1(1:ne,:,:,:), 2, ne/2, nx, nch, num_PIE);
                    M1 = shiftdim(sum(M1),1);
                    
                    m2 = min([8 size(M1,1)]);
                    M0(m1+(1:m2),:,:,:) = M1(1:m2,:,:,:);
                    m1  = m1+m2;
                end
                stim(a,:,:,:,:) = permute(M0(1:m1,:,:,:), [2 3 1 4]); %  stim :  nx nch nbin num_PIE
            end
            
            fclose(fid1);            
            delete([pathname fname(1:end-4) '_FLIM.dat']);
        end
        
        save([pathname fname(1:end-4) '_FLIM.mat'], 'stim', '-v7.3');
        
        clear head tag tcspc IRF stim
        
    end
end

clear all;
