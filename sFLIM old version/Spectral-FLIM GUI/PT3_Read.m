function [sync, tcspc, chan, num, loc, head] = PT3_Read(name, cnts)
%
%  function [sync, tcspc, chan, num, head] = PT3_Read(name, cnts)
%
%  This function reads single-photon data from the file 'name'
%
%  If the argument parameter 'cnts' is missing, or has the value 0 it returns just the 
%  header of the file in the output variable 'sync'.
% 
%  If 'cnts' contains a number larger than 0, the routine reads 'cnts'
%  records the data stream or up the end of the file.
%
%  If 'cnts' contains two numbers [cnts(1) cnts(2)], the routine proceeds
%  to the position cnts(1) before readinf the cnts(2) records of data.
%
%  If 'cnts' contains three numbers [cnts(1) cnts(2) cnts(3)], the routine
%  does not read the header, but jumps over the first cnts(3)-bytes that
%  represent the length of the file-header. After that it procceds to the
%  data position cnts(1) before reading the cnts(2) records of data.
%
%  The output variables contain the followig data:
%  sync    : number of the sync events that preceeded this detection event
%  tcspc   : number of the tcspc-bin of the event
%  chan    : number of the input channel of the event (detector-number)
%  num     : counter of the records that were actually read
%  loc     : number of overcounts after last valid photon
%  head    : the file header 
%

if (nargin<2)||isempty(cnts)
    cnts = [0 0 0];
end;

if numel(cnts)<2
    cnts = [0 cnts 0];
end;

if numel(cnts)<3
    cnts = [cnts 0];
end;


if (nargin<1)||isempty(name)
    fprintf(1,'\n\n      You have to specify a valid file-name. Aborted.\n');
    return;
else    
    fid = fopen(name);
    if fid<1
        fprintf(1,'\n\n      Could not open <%s>. Aborted.\n', name);
        return;
    end
end

%
% The following represents the readable ASCII file header portion 
%

if (nargout==7)||(cnts(3)==0)
    
    head = [];
    
    HeaderLength = fread(fid, 1, 'int32');
    
    if HeaderLength > 0
        fseek(fid,0,'bof');
    end;

    tmp = deblank(char(fread(fid, 16, 'char')'));

    if strcmp(tmp, 'PicoHarp 300')
        
        head.Ident = tmp;    
        head.FormatVersion = fread(fid, 6, '*char')';
        head.CreatorName = fread(fid, 18, '*char')';
        head.CreatorVersion = fread(fid, 12, '*char')';
        head.FileTime = fread(fid, 18, '*char')';
        fread(fid, 2, 'char');
        head.Comment = fread(fid, 256, '*char')';

        head.NCurves = fread(fid, 1, 'int32');
        head.NChannels = 4096;
        head.BitsPerRecord = fread(fid, 1, 'int32');
        head.RoutingChannels = fread(fid, 1, 'int32');
        head.NumberOfBoards = fread(fid, 1, 'int32');
        head.ActiveCurve = fread(fid, 1, 'int32');
        head.MeasMode = fread(fid, 1, 'int32');
        head.SubMode = fread(fid, 1, 'int32');
        head.RangeNo = fread(fid, 1, 'int32');
        head.Offset = fread(fid, 1, 'int32');
        head.TAcq = fread(fid, 1, 'int32');
        head.StopAt = fread(fid, 1, 'int32');
        head.StopOnOvfl = fread(fid, 1, 'int32');
        head.Restart = fread(fid, 1, 'int32');
        head.LinLog = fread(fid, 1, 'int32');
        head.MinAx = fread(fid, 1, 'int32');
        head.MaxAx = fread(fid, 1, 'int32');
        head.MinAxCnt = fread(fid, 1, 'int32');
        head.MaxAxCnt = fread(fid, 1, 'int32');
        head.DispCurves = fread(fid, 16, 'int32');
        head.Params = fread(fid, 9, 'int32');
        head.RepeatMode = fread(fid, 1, 'int32');
        head.RepeatsPerCurve = fread(fid, 1, 'int32');
        head.RepeatTime = fread(fid, 1, 'int32');
        head.RepeatWaitTime = fread(fid, 1, 'int32');
        
        head.ScriptName = fread(fid, 20, '*char')';
        head.anzch = 0;
        for n=(1:head.NumberOfBoards)
            head.BoardIdent(:,n)   = fread(fid, 16, '*char')';
            head.BoardVersion(:,n) = fread(fid, 8, '*char')';
            head.BoardSerial(n)    = fread(fid, 1, 'int32');
            head.SyncDiv(n)        = fread(fid, 1, 'int32');
            head.CFDZeroCross0(n)  = fread(fid, 1, 'int32');
            head.CFDLevel0(n)      = fread(fid, 1, 'int32');
            head.CFDZeroCross1(n)  = fread(fid, 1, 'int32');
            head.CFDLevel1(n)      = fread(fid, 1, 'int32');
            head.Resolution(n)     = fread(fid, 1, 'float');

            head.anzch = head.anzch + 1; 

            if strcmp(head.FormatVersion(1:3),'2.0')
                head.RouterModel(n)        = fread(fid, 1, 'int32');
                head.RouterEnabled(n)      = fread(fid, 1, 'int32');
                head.RtCh1_InputType(n)    = fread(fid, 1, 'int32');
                head.RtCh1_InputLevel(n)   = fread(fid, 1, 'int32');
                head.RtCh1_InputEdge(n)    = fread(fid, 1, 'int32');
                head.RtCh1_CFDPresent(n)   = fread(fid, 1, 'int32');
                head.RtCh1_CFDLevel(n)     = fread(fid, 1, 'int32');
                head.RtCh1_CFDZeroCross(n) = fread(fid, 1, 'int32');
                head.RtCh2_InputType(n)    = fread(fid, 1, 'int32');
                head.RtCh2_InputLevel(n)   = fread(fid, 1, 'int32');
                head.RtCh2_InputEdge(n)    = fread(fid, 1, 'int32');
                head.RtCh2_CFDPresent(n)   = fread(fid, 1, 'int32');
                head.RtCh2_CFDLevel(n)     = fread(fid, 1, 'int32');
                head.RtCh2_CFDZeroCross(n) = fread(fid, 1, 'int32');
                head.RtCh3_InputType(n)    = fread(fid, 1, 'int32');
                head.RtCh3_InputLevel(n)   = fread(fid, 1, 'int32');
                head.RtCh3_InputEdge(n)    = fread(fid, 1, 'int32');
                head.RtCh3_CFDPresent(n)   = fread(fid, 1, 'int32');
                head.RtCh3_CFDLevel(n)     = fread(fid, 1, 'int32');
                head.RtCh3_CFDZeroCross(n) = fread(fid, 1, 'int32');
                head.RtCh4_InputType(n)    = fread(fid, 1, 'int32');
                head.RtCh4_InputLevel(n)   = fread(fid, 1, 'int32');
                head.RtCh4_InputEdge(n)    = fread(fid, 1, 'int32');
                head.RtCh4_CFDPresent(n)   = fread(fid, 1, 'int32');
                head.RtCh4_CFDLevel(n)     = fread(fid, 1, 'int32');
                head.RtCh4_CFDZeroCross(n) = fread(fid, 1, 'int32');
                if head.RouterEnabled(n) == 1
                    head.anzch = head.anzch + 3;
                else 
                    head.anzch = head.RoutingChannels;
                end
            else
                head.anzch      = head.RoutingChannels;
                head.Resolution = head.Resolution.*ones(4,1);
            end
        end;
        
        head.ExternalDev = fread(fid, 1, 'int32');
        head.Reserved1 = fread(fid, 1, 'int32');
        head.Reserved2 = fread(fid, 1, 'int32');
        head.SyncRate = fread(fid, 1, 'int32');
        head.CntRate1 = fread(fid, 1, 'int32');
        head.StopAfter = fread(fid, 1, 'int32');
        head.StopReason = fread(fid, 1, 'int32');
        head.NCounts = fread(fid, 1, 'uint32');

        head.SpecHeaderLength = fread(fid, 1, 'int32');

        tmp = head.SpecHeaderLength/4;
        
        if head.SpecHeaderLength>0
            head.ImgHdr.Dimensions = fread(fid, 1, 'int32');
            head.ImgHdr.Ident = fread(fid, 1, 'int32');
            if head.ImgHdr.Ident==1 % PI E710 Scan Controller
                head.ImgHdr.ScanTimePerPix = fread(fid, 1, 'int32');
                fread(fid, 1, 'int32');
                head.ImgHdr.ScanPattern = fread(fid, 1, 'int32');
                fread(fid, 1, 'int32');
                head.ImgHdr.ScanStartX = fread(fid, 1, 'float');
                head.ImgHdr.ScanStartY = fread(fid, 1, 'float');
                head.ImgHdr.ScanWidthX = fread(fid, 1, 'int32');
                head.ImgHdr.ScanWidthY = fread(fid, 1, 'int32');
                head.ImgHdr.ScanResolution = fread(fid, 1, 'float');
                head.ImgHdr.ScanTStartTo = fread(fid, 1, 'float');
                head.ImgHdr.ScanTStopTo = fread(fid, 1, 'float');
                head.ImgHdr.ScanTStartFro = fread(fid, 1, 'float');
                head.ImgHdr.ScanTStopFro = fread(fid, 1, 'float');
            end
            
            if head.ImgHdr.Ident==3 % LSM
                head.ImgHdr.Frame      = 2^(fread(fid, 1, 'int32')-1);
                head.ImgHdr.LineStart  = 2^(fread(fid, 1, 'int32')-1);
                head.ImgHdr.LineStop   = 2^(fread(fid, 1, 'int32')-1);
                head.ImgHdr.Pattern    = fread(fid, 1, 'int32');
                head.ImgHdr.ScanWidthX = fread(fid, 1, 'int32');
                head.ImgHdr.ScanWidthY = fread(fid, 1, 'int32');
                for n = 9:tmp
                    tmp = fread(fid, 1, 'uint32');
                end
            end
        end
        head.length = ftell(fid);
    end
end


if cnts(2)>0

    if cnts(3)>0
        fseek(fid, (cnts(3)), 'bof');
    end
    if cnts(1)>0
        fseek(fid, 4*(cnts(1)-1), 'cof');
    end

    WRAPAROUND=65536;

    [T3Record num] = fread(fid, cnts(2), 'ubit32'); % all 32 bits:
    
    chan    = bitand(bitshift(T3Record,-28),15);   
    tcspc   = bitand(bitshift(T3Record,-16),4095); 
    sync    = bitand(T3Record,65535); 

    ind  = (tcspc==0 & chan==15);
    sync = sync + WRAPAROUND*cumsum(ind);    
    
    sync(ind)    = [];
    tcspc(ind)   = [];
    chan(ind)    = [];
    loc = num - find(ind==0,1,'last');    
else
    sync = head;
end

fclose(fid);
