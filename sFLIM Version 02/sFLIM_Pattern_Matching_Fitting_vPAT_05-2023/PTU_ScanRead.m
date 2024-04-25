function [head, im_tcspc, im_chan, im_posy, im_posx] = PTU_ScanRead(name, plt)

if (nargin>1)&&~isempty(plt)
    plt = 1;
else
    plt = 0;
end

if strcmp(name(end-2:end),'ptu')
    head  = PTU_Read_Head(name);

    if ~isempty(head)

        nx     = head.ImgHdr_PixX;
        ny     = head.ImgHdr_PixY;
        if ~isfield(head,'ImgHdr_X0')
            head.ImgHdr_X0 = 0;
            head.ImgHdr_Y0 = 0;
        end

        LStart = 1;
        LStop  = 2;
        Frame  = 0;

        if isfield(head,'ImgHdr_LineStart')
            LStart = 2^(head.ImgHdr_LineStart-1);
        end
        if isfield(head,'ImgHdr_LineStop')
            LStop = 2^(head.ImgHdr_LineStop-1);
        end
        if isfield(head,'ImgHdr_Frame')
            Frame = 2^(head.ImgHdr_Frame-1);
        end

        [im_data, ~, im_param] = PTU_ReadPos(name, [nx ny LStart LStop Frame head.ImgHdr_BiDirect head.length head.TTResultFormat_TTTRRecType(1)]);

        %  im_data contains  detector-no, tcspc-bin, x-, and y-position
        %  in 16 bit blocks

        im_data  = uint64(im_data);
        im_chan  = double(bitand(bitshift(im_data,-48),65535))+1;
        im_tcspc = double(bitand(bitshift(im_data,-32),65535))+1;
        im_posy  = double(bitand(bitshift(im_data,-16),65535));   % this is the line no.
        im_posx  = ceil(double(bitand(im_data,65535))/65536*nx);        % this is the fraction of the current line.

        ind = (im_posx<1)|(im_posy<1)|(im_posx>nx)|(im_posy>ny);
        im_posx(ind)  = [];
        im_posy(ind)  = [];
        im_tcspc(ind) = [];
        im_chan(ind)  = [];

        clear im_data;

        % im_param contains some statistics of the scan

        head.ImgHdr_FrameNum  = im_param(1);  % how many frames were scanned
        head.ImgHdr_LineNum   = im_param(2);  % how many lines were completed
        head.ImgHdr_PixNum    = im_param(3);  % how many pixels wer scanned

        head.ImgHdr_FrameTime = im_param(4)/head.TTResult_SyncRate;  % duration of one frame
        head.ImgHdr_LineTime  = im_param(5)/head.TTResult_SyncRate;  % duration of one line
        head.ImgHdr_PixelTime = im_param(6)/head.TTResult_SyncRate;  % accumulated time per pixel
        head.ImgHdr_DwellTime = im_param(7)/head.TTResult_SyncRate;  % pixel dewll time during one scan
        
    end

end

if ~isempty(head)&&(plt==1)
    
    tag = zeros(nx,ny);
    for y = 1:ny
        ind = (im_posy == y);
        tmp1 = im_posx(ind);
        for x = 1:nx
            tag(y,x) = sum(tmp1 == x);
        end
    end
    
    x = head.ImgHdr_X0+(1:nx)*head.ImgHdr_PixResol;
    y = head.ImgHdr_Y0+(1:ny)*head.ImgHdr_PixResol;
    imagesc(x,y,tag);
    set(gca,'DataAspectRatio', [1,1,1], ...
        'PlotBoxAspectRatio',[1 1 1], ...
        'XDir','normal', ...
        'YDir','reverse');
    xlabel('x / µm');
    ylabel('y / µm');
end