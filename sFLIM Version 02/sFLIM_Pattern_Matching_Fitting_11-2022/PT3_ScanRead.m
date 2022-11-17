function [head, im_tcspc, im_chan, im_line, im_col] = PT3_ScanRead(name, plt)

if (nargin>1)&&~isempty(plt)
    plt = 1;
else
    plt = 0;
end;

if strcmp(name(end-2:end),'pt3')

    head = PT3_Read(name);

    if ~isempty(head)

        nx = head.ImgHdr.ScanWidthX;
        ny = head.ImgHdr.ScanWidthY;
                
        y        = [];
        tmpx     = [];
        chan     = [];
        
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

        if head.ImgHdr.Ident == 1  % E710 Piezo Scanner

            if head.ImgHdr.ScanPattern == 1   % unidirectional scan
                
                [tmpy, tmptcspc, tmpchan, num, loc] = PT3_Read(name, [cnt+1 5e6 head.length]);
                
                while (num>0)
                    
                    cnt = cnt + num;

                    tmpy = tmpy+tend;
                                        
                    y       = [y; tmpy];         %#ok<AGROW>
                    tmpx    = [tmpx; tmptcspc];  %#ok<AGROW>
                    chan    = [chan; tmpchan];   %#ok<AGROW>                    
                    Turns   = [Turns; y(chan==15 & tmpx==2)]; %#ok<AGROW>

                    ind       = (chan==15);
                    y(ind)    = [];
                    tmpx(ind) = [];
                    chan(ind) = [];
                    
                    tend  = y(end)+loc;
                    
                    if numel(Turns)>1
                        for j=1:numel(Turns)-1
                            
                            dT = (Turns(2)-Turns(1));
                            t1 = Turns(1)+head.ImgHdr.TStartTo*dT;
                            t2 = Turns(1)+head.ImgHdr.TStopTo*dT;
                            
                            ind = (y<t1);
                            y(ind)    = [];
                            tmpx(ind) = [];
                            chan(ind) = [];
                            
                            ind = (y>=t1)&(y<=t2);
                            
                            im_tcspc  = [im_tcspc; uint16(tmpx(ind))];                            %#ok<AGROW>
                            im_chan   = [im_chan;  uint8(chan(ind))];                             %#ok<AGROW>
                            im_line   = [im_line;  uint16(line.*ones(sum(ind),1))];               %#ok<AGROW>
                            im_col    = [im_col;   uint16(1 + floor(nx.*(y(ind)-t1)./(t2-t1)))];  %#ok<AGROW>
                            
                            line = line +1;
                            waitbar(line/ny);
                            drawnow
                            
                            Turns(1) = [];                            
                        end
                        
                    end
                    
                    [tmpy, tmptcspc, tmpchan, num, loc] = PT3_Read(name, [cnt+1 5e6 head.length]);
                    
                end;
                
                t1 = Turns(end)+head.ImgHdr.TStartTo*dT;
                t2 = Turns(end)+head.ImgHdr.TStopTo*dT;
                
                ind          = (y<t1);
                y(ind)       = [];
                tmpx(ind)    = [];
                chan(ind)    = [];
                
                ind = (y>=t1)&(y<=t2);
                
                im_tcspc  = [im_tcspc; uint16(tmpx(ind))];
                im_chan   = [im_chan;  uint8(chan(ind))];
                im_line   = [im_line;  uint16(line.*ones(sum(ind),1))];
                im_col    = [im_col;   uint16(1 + floor(nx.*(y(ind)-t1)./(t2-t1)))];
                
                line = line +1;
                waitbar(line/ny);
                drawnow
                
            else % bidirectional scan
                
                [tmpy, tmptcspc, tmpchan, num, loc] = PT3_Read(name, [cnt+1 5e6 head.length]);
                
                while (num>0)
                    
                    cnt = cnt + num;

                    tmpy = tmpy+tend;
                                        
                    y       = [y;    tmpy];        %#ok<AGROW>
                    tmpx    = [tmpx; tmptcspc];    %#ok<AGROW>
                    chan    = [chan; tmpchan];     %#ok<AGROW>
                    
                    Turns   = [Turns; y(chan==15 & tmpx==2)]; %#ok<AGROW>
                    
                    ind       = (chan==15);
                    y(ind)    = [];
                    tmpx(ind) = [];
                    chan(ind) = [];

                    tend = y(end)+loc;
                    
                    if numel(Turns)>2
                        for j=1:2:2*floor(numel(Turns)/2-1)
                            
                            dT = (Turns(2)-Turns(1));
                            t1 = Turns(1)+head.ImgHdr.TStartTo*dT;
                            t2 = Turns(1)+head.ImgHdr.TStopTo*dT;
                            
                            ind = (y<t1);
                            y(ind)    = [];
                            tmpx(ind) = [];
                            chan(ind) = [];
                            
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
                    
                    [tmpy, tmptcspc, tmpchan, num, loc] = PT3_Read(name, [cnt+1 5e6 head.length]);
                    
                end;
                
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
            
        elseif head.ImgHdr.Ident == 3   % LSM
                        
            head.ImgHdr.X0        = 0;
            head.ImgHdr.Y0        = 0;
            head.ImgHdr.PixX      = nx;
            head.ImgHdr.PixY      = ny;
            head.ImgHdr.PixelSize = 1;
            
            LineStart = head.ImgHdr.LineStart;
            LineStop  = head.ImgHdr.LineStop;
            Frame     = head.ImgHdr.Frame;
            
            if Frame < 1
                Frame = -1;
            end
            
            in_frame = false;
            
            if Frame < 1
                in_frame = true;
            end
            
            [tmpy, tmptcspc, tmpchan, num, loc] = PT3_Read(name, [cnt+1 1e6 head.length]);
                                                    
            while (num>0)
                
                t_tcspc  = [];
                t_chan   = [];
                t_line   = [];
                t_col    = [];
                
                cnt = cnt + num;
                
                waitbar(cnt/head.NCounts);
                
                tmpy = tmpy+tend;
                
                y       = [y; tmpy];                       %#ok<AGROW>
                tmpx    = [tmpx; tmptcspc];                %#ok<AGROW>
                chan    = [chan; tmpchan];                 %#ok<AGROW>
                tend    = y(end)+loc;
                
                F  = y(chan==15 & bitand(tmpx,Frame)>0);
                
                while ~isempty(F)
                    
                    if ~in_frame
                        ind = (y<=F(1));
                        y(ind)       = [];
                        tmpx(ind)    = [];
                        chan(ind)    = [];
                        line         = 1;
                        in_frame     = true;
                        F(1)         = [];
                    end
                    
                    if ~isempty(F)
                        ind = y<F(1);
                        
                        f_y  = y(ind);
                        f_x  = tmpx(ind);
                        f_ch = chan(ind);
                        
                        y(ind)    = [];
                        tmpx(ind) = [];
                        chan(ind) = [];
                        F(1)      = [];
                    end
                                        
                    L1 = f_y(f_ch==15 & bitand(f_x,LineStart)>0);
                    L2 = f_y(f_ch==15 & bitand(f_x,LineStop)>0);
                                        
                    ll = line + numel(L2)-1; % this will be the last complete line in the data stack
                    
                    if ll > ny
                        L1 = L1(1:ny-line+1);
                        L2 = L2(1:ny-line+1);
                    end
                    
                    if numel(L1)>1
                        for j=1:numel(L2)
                            
                            ind = (f_y>L1(j))&(f_y<L2(j));
                            
                            t_tcspc  = [t_tcspc; uint16(f_x(ind))];              %#ok<AGROW>
                            t_chan   = [t_chan; uint8(f_ch(ind))];                %#ok<AGROW>
                            t_line   = [t_line; uint16(line.*ones(sum(ind),1))];  %#ok<AGROW>
                            t_col    = [t_col;  uint16(1 + floor(nx.*(f_y(ind)-L1(j))./(L2(j)-L1(j))))];  %#ok<AGROW>
                            
                            line = line +1;
                                                        
                        end
                    end
                    
                    if line>ny
                        in_frame = false;
                    end
                end
                
                im_tcspc  = [im_tcspc; t_tcspc];  %#ok<AGROW>
                im_chan   = [im_chan;  t_chan];   %#ok<AGROW>
                im_line   = [im_line;  t_line];   %#ok<AGROW>
                im_col    = [im_col;   t_col];    %#ok<AGROW>
                
                [tmpy, tmptcspc, tmpchan, num, loc] = PT3_Read(name, [cnt+1 1e6 head.length]);
                
            end;
            
            F  = y(chan==15 & bitand(tmpx,Frame)>0);

            t_tcspc  = []; 
            t_chan   = [];
            t_line   = [];
            t_col    = [];
                                        
            if ~in_frame
                if isempty(F)
                    y       = [];
                    tmpx    = [];
                    chan    = [];
                    line    = 1;
                else
                    ind = (y<=F(1));
                    y(ind)       = [];
                    tmpx(ind)    = [];
                    chan(ind)    = [];
                    line         = 1;
                end
            end
            
            f_y = y;
            f_x = tmpx;
            f_ch = chan;
            
            clear y tmpx chan;
            
            L1 = f_y(f_ch==15 & bitand(f_x,LineStart)>0);
            L2 = f_y(f_ch==15 & bitand(f_x,LineStop)>0);
            
            ll = line + numel(L2)-1; % this will be the last complete line in the data stack
            if ll > ny
                L1 = L1(1:ny-line+1);
                L2 = L2(1:ny-line+1);
            end
            
            if numel(L1)>1
                for j=1:numel(L2)
                    ind = (f_y>L1(j))&(f_y<L2(j));
                    t_tcspc  = [t_tcspc; uint16(f_x(ind))];              %#ok<AGROW>
                    t_chan   = [t_chan; uint8(f_ch(ind))];                %#ok<AGROW>
                    t_line   = [t_line; uint16(line.*ones(sum(ind),1))];  %#ok<AGROW>
                    t_col    = [t_col;  uint16(1 + floor(nx.*(f_y(ind)-L1(j))./(L2(j)-L1(j))))];  %#ok<AGROW>
                    line = line +1;
                end
            end    
            
            im_tcspc  = [im_tcspc; t_tcspc];  
            im_chan   = [im_chan;  t_chan];   
            im_line   = [im_line;  t_line];   
            im_col    = [im_col;   t_col];    

        end
    end;

    close(h);
 
end


if ~isempty(head)&&(plt==1)
    
    tag = zeros(nx,ny);
    for y = 1:ny
        ind = (im_line == y);
        tmp1 = im_col(ind);
        for x = 1:nx
            ind = (tmp1 == x);
            tag(x,y) = sum(ind);
        end
    end
    
    figure;
    x = head.ImgHdr.X0+(1:nx)*head.ImgHdr.PixelSize;
    y = head.ImgHdr.Y0+(1:ny)*head.ImgHdr.PixelSize;
    
    imagesc(x,y,sum(sum(tag,4),3));
    set(gca,'DataAspectRatio', [1,1,1], ...
        'PlotBoxAspectRatio',[1 1 1], ...
        'XDir','normal', ...
        'YDir','reverse');
    xlabel('x / µm');
    ylabel('y / µm');
end