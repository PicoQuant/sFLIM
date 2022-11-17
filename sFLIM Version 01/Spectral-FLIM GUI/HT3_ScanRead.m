function [head, im_tcspc, im_chan, im_line, im_col] = HT3_ScanRead(name, plt)

if (nargin>1)&&~isempty(plt)
    plt = 1;
else
    plt = 0;
end;

if strcmp(name(end-2:end),'ht3')

    head = HT3_Read(name);

    if ~isempty(head)

        nx = head.ImgHdr.PixX;
        ny = head.ImgHdr.PixY;
        dt = 1e-9*nx*head.ImgHdr.PixelTime*head.SyncRate;
                        
        anzch      = 32;
        Resolution = max([head.Resolution]);
        chDiv      = Resolution/head.Resolution;
        Ngate      = ceil(1e9/head.SyncRate./Resolution)+1;
                
        y        = [];
        tmpx     = [];
        chan     = [];
        markers  = [];
        
        im_tcspc = [];
        im_chan  = [];
        im_line  = [];
        im_col   = [];
        Turns    = [];

        cnt      = 0;
        tend     = 0; 
        line     = 1;

        h = waitbar(0,'Please wait ...');

        fprintf('\n\n');

        if head.ImgHdr.Pattern == 0

            [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = HT3_Read(name, [cnt+1 5e6 head.length]);

            while (num>0)

                cnt = cnt + num;                
                if ~isempty(y)
                    tmpy = tmpy+tend;                    
                end;
                                              
                ind = (tmpmarkers>0)|((tmpchan<anzch)&(tmptcspc<Ngate*chDiv));
                
                y       = [y; tmpy(ind)];                         %#ok<AGROW>
                tmpx    = [tmpx; floor(tmptcspc(ind)./chDiv)+1;]; %#ok<AGROW>
                chan    = [chan; tmpchan(ind)+1];                 %#ok<AGROW>
                markers = [markers; tmpmarkers(ind)];             %#ok<AGROW>
                
                Turns   = [Turns; y(markers==1 & chan==5)]; %#ok<AGROW>
                ind     = (markers~=0);
                y(ind)       = [];
                tmpx(ind)    = [];
                chan(ind)    = [];
                markers(ind) = [];
                     
                tend  = y(end)+loc;

                if numel(Turns)>1
                    for j=1:numel(Turns)-1

                        dT = (Turns(2)-Turns(1));
                        t1 = Turns(1)+head.ImgHdr.TStartTo*dT;
                        t2 = Turns(1)+head.ImgHdr.TStopTo*dT;
                        
                        ind = (y<t1);
                        y(ind)       = [];
                        tmpx(ind)    = [];
                        chan(ind)    = [];
                        markers(ind) = [];
                        
                        ind = (y>=t1)&(y<=t2);
                                      
                        im_tcspc  = [im_tcspc; uint16(tmpx(ind))];              %#ok<AGROW>
                        im_chan   = [im_chan; uint8(chan(ind))];                %#ok<AGROW>
                        im_line   = [im_line; uint16(line.*ones(sum(ind),1))];  %#ok<AGROW>
                        im_col    = [im_col;  uint16(1 + floor(nx.*(y(ind)-t1)./(t2-t1)))];  %#ok<AGROW>
                     
                        line = line +1;
                        waitbar(line/ny);
                        drawnow
                        
                        Turns(1) = [];
                        
                    end                    
                end
                
                [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = HT3_Read(name, [cnt+1 5e6 head.length]);

            end;
            
            t1 = Turns(end)+head.ImgHdr.TStartTo*dT;
            t2 = Turns(end)+head.ImgHdr.TStopTo*dT;

            ind          = (y<t1);
            y(ind)       = [];
            tmpx(ind)    = [];
            chan(ind)    = [];

            ind = (y>=t1)&(y<=t2);

            im_tcspc  = [im_tcspc; uint16(tmpx(ind))];              
            im_chan   = [im_chan; uint8(chan(ind))];                
            im_line   = [im_line; uint16(line.*ones(sum(ind),1))];  
            im_col    = [im_col;  uint16(1 + floor(nx.*(y(ind)-t1)./(t2-t1)))];  

            line = line +1;
            waitbar(line/ny);
            drawnow
            
        else  % bidirectional scan

            [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = HT3_Read(name, [cnt+1 5e6 head.length]);

            while (num>0)

                cnt = cnt + num;
                if ~isempty(y)
                    tmpy = tmpy+tend;                    
                end;

                ind = (tmpmarkers>0)|((tmpchan<anzch)&(tmptcspc<=Ngate*chDiv));
                
                y       = [y; tmpy(ind)];                         %#ok<AGROW>
                tmpx    = [tmpx; floor(tmptcspc(ind)./chDiv)+1;]; %#ok<AGROW>
                chan    = [chan; tmpchan(ind)+1];                 %#ok<AGROW>
                markers = [markers; tmpmarkers(ind)];             %#ok<AGROW>
                                
                Turns = [Turns; y(markers==1 & chan==2)]; %#ok<AGROW>
                
                ind = (markers~=0);
                y(ind)       = [];
                tmpx(ind)    = [];
                chan(ind)    = [];
                markers(ind) = [];
                                
                tend = y(end)+loc;

                if numel(Turns)>2
                    for j=1:2:2*floor(numel(Turns)/2-1)

                        dT = (Turns(2)-Turns(1));
                        t1 = Turns(1)+head.ImgHdr.TStartTo*dT;
                        t2 = Turns(1)+head.ImgHdr.TStopTo*dT;

                        ind = (y<t1);
                        y(ind)       = [];
                        tmpx(ind)    = [];
                        chan(ind)    = [];
                        markers(ind) = [];
                        
                        ind = (y>=t1)&(y<=t2);
                                      
                        im_tcspc  = [im_tcspc; uint16(tmpx(ind))];              %#ok<AGROW>
                        im_chan   = [im_chan; uint8(chan(ind))];                %#ok<AGROW>
                        im_line   = [im_line; uint16(line.*ones(sum(ind),1))];  %#ok<AGROW>
                        im_col    = [im_col;  uint16(1 + floor(nx.*(y(ind)-t1)./(t2-t1)))];  %#ok<AGROW>
                                                                                               
                        line = line +1;

                        t1 = Turns(1)+head.ImgHdr.TStartFro*dT;
                        t2 = Turns(1)+head.ImgHdr.TStopFro*dT;

                        ind = (y<t1);
                        y(ind)       = [];
                        tmpx(ind)    = [];
                        chan(ind)    = [];
                        markers(ind) = [];
                        
                        ind = (y>=t1)&(y<=t2);
                                      
                        im_tcspc  = [im_tcspc; uint16(tmpx(ind))];              %#ok<AGROW>
                        im_chan   = [im_chan; uint8(chan(ind))];                %#ok<AGROW>
                        im_line   = [im_line; uint16(line.*ones(sum(ind),1))];  %#ok<AGROW>
                        im_col    = [im_col;  uint16(nx - floor(nx.*(y(ind)-t1)./(t2-t1)))];  %#ok<AGROW>
                        
                        line = line +1;
                        waitbar(line/ny);
                        drawnow

                        Turns(1:2) = [];
                    end
                end
                
                [tmpy, tmptcspc, tmpchan, tmpmarkers, num, loc] = HT3_Read(name, [cnt+1 5e6 head.length]);

            end;
            
            if ~isempty(Turns)
                t1 = Turns(end-1)+head.ImgHdr.TStartTo*dT;
                t2 = Turns(end-1)+head.ImgHdr.TStopTo*dT;
                
                ind = (y<t1);
                y(ind)       = [];
                tmpx(ind)    = [];
                chan(ind)    = [];
                
                ind = (y>=t1)&(y<=t2);
                
                im_tcspc  = [im_tcspc; uint16(tmpx(ind))];
                im_chan   = [im_chan; uint8(chan(ind))];
                im_line   = [im_line; uint16(line.*ones(sum(ind),1))];
                im_col    = [im_col;  uint16(1 + floor(nx.*(y(ind)-t1)./(t2-t1)))];
                
                line = line +1;
                
                t1 = Turns(end-1)+head.ImgHdr.TStartFro*dT;
                t2 = Turns(end-1)+head.ImgHdr.TStopFro*dT;
                
                ind = (y<t1);
                y(ind)       = [];
                tmpx(ind)    = [];
                chan(ind)    = [];
                
                ind = (y>=t1)&(y<=t2);
                
                im_tcspc  = [im_tcspc; uint16(tmpx(ind))];
                im_chan   = [im_chan; uint8(chan(ind))];
                im_line   = [im_line; uint16(line.*ones(sum(ind),1))];
                im_col    = [im_col;  uint16(nx - floor(nx.*(y(ind)-t1)./(t2-t1)))];
                
                line = line +1;
                waitbar(line/ny);
                drawnow
            end
        end;

        close(h);

    end
end


if ~isempty(head)&&(plt==1)
    
    tag = zeros(nx,ny);
    for y = 1:ny
        ind = (im_line == y);
        tmp1 = im_col(ind);
        for x = 1:nx
            tag(y,x) = sum(tmp1 == x);
        end
    end

    x = head.ImgHdr.X0+(1:nx)*head.ImgHdr.PixelSize;
    y = head.ImgHdr.Y0+(1:ny)*head.ImgHdr.PixelSize;
    imagesc(x,y,tag);
    set(gca,'DataAspectRatio', [1,1,1], ...
        'PlotBoxAspectRatio',[1 1 1], ...
        'XDir','normal', ...
        'YDir','reverse');
    xlabel('x / µm');
    ylabel('y / µm');
end