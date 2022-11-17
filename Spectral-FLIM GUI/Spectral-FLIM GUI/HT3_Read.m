function [sync, tcspc, chan, special, num, loc, head] = HT3_Read(name, cnts)
%
%  function [sync, tcspc, chan, special, num, head] = HT3_Read(name, cnts)
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
%  special : indicator of the event-type (0: photon; else : virtual photon)
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
    
    tmp = deblank(char(fread(fid, 16, 'char')'));

    if strcmp(tmp, 'HydraHarp')
        
        head.Ident = tmp;
        head.FormatVersion = deblank(char(fread(fid, 6, 'char')'));

%         if not(strcmp(head.FormatVersion,'1.0'))
%             fprintf(1,'\n\n      Warning: This program is for version 1.0 only. Aborted.');
%             fclose(fid);
%             return;
%         end;

        head.CreatorName = deblank(char(fread(fid, 18, 'char')'));
        head.CreatorVersion = deblank(char(fread(fid, 12, 'char')'));
        head.FileTime = deblank(char(fread(fid, 18, 'char')'));
        fread(fid, 2, 'char');
        head.CommentField = deblank(char(fread(fid, 256, 'char')'));

        %
        % The following is binary file header information
        %

        head.Curves = fread(fid, 1, 'int32');
        head.BitsPerRecord = fread(fid, 1, 'int32');
        head.ActiveCurve = fread(fid, 1, 'int32');
        head.MeasurementMode = fread(fid, 1, 'int32');
        head.SubMode = fread(fid, 1, 'int32');
        head.Binning = fread(fid, 1, 'int32');
        head.Resolution = fread(fid, 1, 'double')/1000;
        head.Offset = fread(fid, 1, 'int32');
        head.AcquisitionTime = fread(fid, 1, 'int32');
        head.StopAt = fread(fid, 1, 'int32');
        head.StopOnOvfl = fread(fid, 1, 'int32');
        head.Restart = fread(fid, 1, 'int32');
        head.DispLinLog = fread(fid, 1, 'int32');
        head.DispTimeFrom = fread(fid, 1, 'int32');
        head.DispTimeTo = fread(fid, 1, 'int32');
        head.DispCountFrom = fread(fid, 1, 'int32');
        head.DispCountTo = fread(fid, 1, 'int32');

        for i = 1:8
            head.DispCurveMapTo(i) = fread(fid, 1, 'int32');
            head.DispCurveShow(i) = fread(fid, 1, 'int32');
        end;

        for i = 1:3
            head.ParamStart(i) = fread(fid, 1, 'float');
            head.ParamStep(i) = fread(fid, 1, 'float');
            head.ParamEnd(i) = fread(fid, 1, 'float');
        end;

        head.RepeatMode = fread(fid, 1, 'int32');
        head.RepeatsPerCurve = fread(fid, 1, 'int32');
        head.RepeatTime = fread(fid, 1, 'int32');
        head.RepeatWait = fread(fid, 1, 'int32');
        head.ScriptName = deblank(char(fread(fid, 20, 'char')'));

        %
        % The next is a board specific header
        %

        head.HardwareIdent = deblank(char(fread(fid, 16, 'char')'));
        head.HardwareVersion = deblank(char(fread(fid, 8, 'char')'));
        head.HardwareSerial = fread(fid, 1, 'int32');
        head.NumModules  = fread(fid, 1, 'int32');

        for i = 1:10
            head.Module(i).Serial = fread(fid, 1, 'int32');
            head.Module(i).Version = fread(fid, 1, 'int32');
        end;

        head.BaseResolution = fread(fid, 1, 'double');
        head.InputsEnabled = fread(fid,1,'int64');
        head.InputsPresent = fread(fid,1,'int32');
        head.RefClock = fread(fid,1,'int32');
        head.ExtDevices = fread(fid,1,'int32');
        head.MarkerSettings = fread(fid,1,'int32');

        head.SyncDivider = fread(fid, 1, 'int32');
        head.SyncCFDLevel = fread(fid, 1, 'int32');
        head.SyncCFDZeroCross = fread(fid, 1, 'int32');
        head.SyncOffset = fread(fid, 1, 'int32');

        for i = 1:head.InputsPresent
            head.InpChannel(i).ModuleIdx = fread(fid, 1, 'int32');
            head.InpChannel(i).CFDLevel = fread(fid, 1, 'int32');
            head.InpChannel(i).CFDZeroCross = fread(fid, 1, 'int32');
            head.InpChannel(i).Offset = fread(fid, 1, 'int32');
        end;

        for i = 1:head.InputsPresent
            head.InpChannel(i).CntRate = fread(fid, 1, 'int32');
        end;

        %
        % The next is a T3 mode specific header
        %

        head.SyncRate = fread(fid, 1, 'int32');
        head.StopAfter = fread(fid, 1, 'int32');
        head.StopReason = fread(fid, 1, 'int32');
        head.ImgHdrSize = fread(fid, 1, 'int32');
        head.Records = fread(fid, 1, 'int64');

        %Special header for imaging

        if head.ImgHdrSize > 0
            head.ImgHdr.Dimensions = fread(fid, 1, 'int32');
            head.ImgHdr.Ident = fread(fid, 1, 'int32');
            if head.ImgHdr.Ident == 1
                head.ImgHdr.PixelTime =  2e5*fread(fid, 1, 'int32');
                head.ImgHdr.Acceleration =  fread(fid, 1, 'int32');
                head.ImgHdr.Pattern =  fread(fid, 1, 'int32');
                fread(fid, 1, 'int32');
                head.ImgHdr.X0 = fread(fid, 1, 'float');
                head.ImgHdr.Y0 = fread(fid, 1, 'float');
                head.ImgHdr.PixX =  fread(fid, 1, 'int32');
                head.ImgHdr.PixY =  fread(fid, 1, 'int32');
                head.ImgHdr.PixelSize =  fread(fid, 1, 'float');
                head.ImgHdr.TStartTo  =  fread(fid, 1, 'float');
                head.ImgHdr.TStopTo   =  fread(fid, 1, 'float');
                head.ImgHdr.TStartFro =  fread(fid, 1, 'float');
                head.ImgHdr.TStopFro  =  fread(fid, 1, 'float');
            elseif head.ImgHdr.Ident == 16
                head.ImgHdr.PixelTime =  2e5*fread(fid, 1, 'int32');
                head.ImgHdr.Acceleration =  fread(fid, 1, 'int32');
                head.ImgHdr.Pattern =  fread(fid, 1, 'int32');
                head.ImgHdr.X0 = fread(fid, 1, 'int32')/1000;
                head.ImgHdr.Y0 = fread(fid, 1, 'int32')/1000;
                head.ImgHdr.Z0 = fread(fid, 1, 'int32')/1000;
                head.ImgHdr.PixX =  fread(fid, 1, 'int32');
                head.ImgHdr.PixY =  fread(fid, 1, 'int32');
                head.ImgHdr.PixelSize =  fread(fid, 1, 'float');
                head.ImgHdr.TStartTo  =  fread(fid, 1, 'float');
                head.ImgHdr.TStopTo   =  fread(fid, 1, 'float');
                head.ImgHdr.TStartFro =  fread(fid, 1, 'float');
                head.ImgHdr.TStopFro  =  fread(fid, 1, 'float');

            else
                head.ImgHdr.data = fread(fid, head.ImgHdrSize-2, 'int32');
            end
        end;
        head.length = ftell(fid);
    else
        head.Resolution = 2e-3;   % in nanoseconds
        head.length = 0;
    end
end


if cnts(2)>0

    if cnts(3)>0
        fseek(fid, (cnts(3)), 'bof');
    end
    if cnts(1)>1
        fseek(fid, 4*(cnts(1)-1), 'cof');
    end

    WRAPAROUND=1024;

    [T3Record num] = fread(fid, cnts(2), 'ubit32'); % all 32 bits:
    special = bitand(T3Record,2147483648)>0;       
    chan = bitand(bitshift(T3Record,-25),63);   
    tcspc = bitand(bitshift(T3Record,-10),32767); 
    sync = bitand(T3Record,1023); 

    ind  = (special==1 & chan==63);
    tmp  = sync(ind==1);
    tmp(tmp==0) = 1;
    sync(ind) = tmp;
    sync = sync + WRAPAROUND*cumsum(ind.*sync);    
    
    sync(ind)    = [];
    tcspc(ind)   = [];
    special(ind) = [];
    chan(ind)    = [];
    loc = num - find(ind==0,1,'last');    
else
    sync = head;
end

fclose(fid);

