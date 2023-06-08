function [tcspc, timegate] = IRF_Read(name, num_PIE)

if (nargin<2)||isempty(num_PIE)
    num_PIE = 1;
end;

N1 = [name(1:end-4) '_DATA.mat'];

tmp = dir(N1);

if (numel(tmp)==0)
    [head, tag, tcspc] = MT_ScanRead(name, num_PIE);
    timegate = head.timegate;
else
    load(N1,'KRF','head');
    tcspc    = KRF;
    timegate = head.timegate;
end

