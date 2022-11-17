function [timegate, Ngate, tcspc] = DetectTimeGates(tcspcdata, cnum, Resolution, pic)
%. 
% 
%  DetectTimeGates divides the TCSPC histogramm into time-windows of the
%  same length in order to separate the data from different excitation
%  pulses.
%
%  usage: [timegate, Ngate, tcspc] = DetectTimeGates(tcspcdata, cnum, Resolution, pic)
%
%    tcspcdata  : TCSPC histogram of a PIE experiment
%    cnum       : number of excitation pulses (default: 1)
%    Resolution : bin-size of TSCPC data in ns (default: 0.002 ns)
%    pic        : flag to display aligned data
%
%  DetectTimeGates returns the following data:
%
%    timegate   : array of the size [cnum*dnum, 4] containing the
%                 respective beginning and end of each time-window. 
%
%                 The order of the list is 
%
%                 timegate(  1   ,:): time-window for pulse  1   in detection channel 1
%                 timegate(  1   ,:): time-window for pulse  1   in detection channel 2
%            ...  timegate(  n   ,:): time-window for pulse  1   in detection channel n
%                 timegate( n+1  ,:): time-window for pulse  2   in detection channel 1 
%                 timegate( n+2  ,:): time-window for pulse  2   in detection channel 2
%            ...  timegate( 2*n  ,:): time-window for pulse  2   in detection channel n
%            ...  timegate(cnum*n,:): time-window for pulse cnum in detection channel n
%
%                 timegate(:, 1)    : begin of time-window 
%                 timegate(:, 2)    : end of time-window 
%                 timegate(:, 3)    : if > 0: continuation of time-window
%                 timegate(:, 4)    : if > 0: end of time-window 
%
%    Ngate      : Number of bins in each time-window
%    tcscpc     : Sorted and aligned TCSPC histogram in the same order as the timegates
%
% .........................................................................................

if (nargin < 4)||isempty(pic)
    pic = 0;
end
pic = (pic>0);

if (nargin < 3)||isempty(Resolution)
    Resolution = 0.002;
end

if (nargin < 2)||isempty(cnum)
    cnum = 1;
end

NChannels = size(tcspcdata,1);
dnum      = size(tcspcdata,2);
timegate  = zeros(dnum*cnum,4); 

%% Determine maxima of all detectors

tmpdata = zeros(size(tcspcdata));

for k = 1:dnum
    tmpdata(:,k) = smooth(tcspcdata(:,k),ceil(0.1./Resolution));
end

irf_w = floor(1.5/Resolution);
win   = floor(NChannels/cnum);     % Max #channels / pulse;

im = zeros(dnum, cnum);
m  = zeros(dnum, cnum);
bg = zeros(dnum, cnum);

[tpm,i1] = max(tmpdata); %#ok<ASGLU>
off      = floor(mod(i1,win)-irf_w);

for l = 1:cnum
    for k = 1:dnum
        ind  = off(k)+(l-1)*win+(1:win-2*irf_w);
        ind(ind<1)         = ind(ind<1) + NChannels;
        ind(ind>NChannels) = ind(ind>NChannels) - NChannels;
        [m(k,l),im(k,l)] = max(tmpdata(ind,k));
        im(k,l) = ind(im(k,l));
        bg(k,l) = mean(tmpdata(ind(2*irf_w:end),k));
    end
end

ind = (m<3*bg)|(m<10);

for l = 1:cnum
    for k = 1:dnum
        if ind(k,l) == 1
            mi = find(im(k,:) == i1(k));
            im(k,l) = im(k,mi)+floor((l-mi)*NChannels/cnum);
        end
    end
end

%% Determine time-gates for photon assignment
if cnum>1
    Ngate = zeros(1,cnum);
    for k = 1:cnum
        p = k-1;
        if p < 1
            p = cnum;
        end
        Ngate(k) = min([abs(im(:,k)-im(:,p)); NChannels-abs(im(:,p)-im(:,k))]);
    end
    Ngate = min(Ngate);
else
    Ngate = NChannels;
end

%% Define final windows for each pulse and detector

g = zeros(cnum*dnum, 4);

for k = 1:dnum
    for l = 1:cnum
        d = im(k,l)-irf_w;
        if d > 0
            if d <= (NChannels - Ngate)+1
                g(k+(l-1)*dnum,1:2) = [d d+Ngate-1];
            else
                g(k+(l-1)*dnum,:) = [d NChannels 1 d+Ngate-NChannels-1];
            end;
        else
            d = d + NChannels;
            if d > NChannels - Ngate +1
                g(k+(l-1)*dnum,:) = [d NChannels 1 d+Ngate-NChannels-1];
            else
                g(k+(l-1)*dnum,1:2) = [d d+Ngate-1];
            end;
        end;
    end
end

%% Sort windows to according excitation pulse

channels = repmat(1:dnum,[1 cnum]);

for n = 1:cnum*dnum
    mch = mod(mean(g(n,1):(g(n,2)+g(n,4))),NChannels);
    pie = floor(mch/(NChannels/(cnum-0.5)));
    timegate(channels(n)+pie*dnum,:) = g(n,:);
end

tcspc = zeros(Ngate,cnum*dnum);

for n = 1:cnum*dnum
    if (timegate(n,3)==0)
        tmp =  tcspcdata(timegate(n,1):timegate(n,2),channels(n));
    else
        tmp = [tcspcdata(timegate(n,1):timegate(n,2),channels(n)); tcspcdata(timegate(n,3):timegate(n,4),channels(n))];
    end
    tcspc(:,n) = tmp;
end

if pic
    figure;
    semilogy((1:Ngate), tcspc);
end

