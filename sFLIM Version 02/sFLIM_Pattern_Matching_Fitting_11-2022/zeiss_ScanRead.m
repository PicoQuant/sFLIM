function [head, im_tcspc, im_chan, im_line, im_col] =  zeiss_ScanRead(sync, tcspc, channel, special, head)


tcspc_bin_resolution = 1e9*head.MeasDesc_Resolution; % in Nanoseconds
sync_rate            = ceil(head.MeasDesc_GlobalResolution*1e9); % in Nanoseconds
num_tcspc_channel    = max(tcspc)+1;
%num_tcspc_channel    = floor(sync_rate/tcspc_bin_resolution)+1;
num_pixel_X          = head.ImgHdr_PixX;
num_pixel_Y          = head.ImgHdr_PixY;
num_of_detectors     = max(channel);


% Get Number of Frames
FrameSyncVal       = sync(special  == 4);
num_of_Frames      = size(FrameSyncVal,1);
read_data_range    = find(sync == FrameSyncVal(num_of_Frames));

% Markers necessary to make FLIM image stack
LineStartMarker = 2^(head.ImgHdr_LineStart-1);
LineStopMarker  = 2^(head.ImgHdr_LineStop-1);
FrameMarker     = 2^(head.ImgHdr_Frame-1);

L1  = sync((special == 1));
L2  = sync((special == 2));

% Get pixel dwell time values from header for PicoQuant_FLIMBee or Zeiss_LSM scanner
syncPulsesPerLine = floor(mean(L2(1:100,1)- L1(1:100,1)));

% Initialize Variable
currentLine        = 0;
currentSync        = 0;
syncStart          = 0;
currentPixel       = 0;
countFrame         = -1;
insideLine  = false;
insideFrame = false;
isPhoton    = false;

im_tcspc = [];
im_chan  = [];
im_line  = [];
im_col   = [];

% Read each event separately, and build the image matrix as you go
for event  = 1:read_data_range
    
    currentSync    = sync(event);
    special_event  = special(event);
    currentChannel = channel(event);
    currentTcspc   = tcspc(event); 
    
    if special(event) == 0
        isPhoton = true;
    else
        isPhoton = false;
    end
    
    if ~isPhoton
        
        if (special_event == FrameMarker)
            
            insideFrame  = true;
            countFrame   = countFrame + 1;
            currentLine  = 1;
            %
            %                 % NEW ADDITION
            %                 if countFrame  == num_of_Frames
            %                     insideFrame = false;
            %                 end
        end
        
        if (special_event == LineStartMarker)% && (insideFrame ==  true))
            
            insideLine = true;
            syncStart  = currentSync;
            
        elseif special_event == LineStopMarker
            
            insideLine   = false;
            currentLine  = currentLine + 1;
            syncStart    = 0;
            if (currentLine > num_pixel_Y)
                insideFrame = false;
                currentLine  = 1;
            end
        end
    end
    %update image_data matrix if it's a valid photon
    if (isPhoton && insideLine && insideFrame)
        
        currentPixel = 1 + floor(num_pixel_X*((currentSync - syncStart)/syncPulsesPerLine));
        
        if currentPixel <= num_pixel_X
            im_tcspc  = [im_tcspc; uint16(currentTcspc)];  %#ok<AGROW>
            im_chan   = [im_chan;  uint8(currentChannel)];   %#ok<AGROW>
            im_line   = [im_line;  uint16(currentLine)];   %#ok<AGROW>
            im_col    = [im_col;   uint16(currentPixel)];    %#ok<AGROW> 
        end
    end
    
    
end % end of looping for all events
end