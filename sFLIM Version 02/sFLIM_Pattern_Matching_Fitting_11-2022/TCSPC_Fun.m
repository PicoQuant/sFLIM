function [err, c, zz, z] = TCSPC_Fun(p, t, y, para)

%  Fits the following model function to TCSPC-Data
%
%  TCSPC(N,{q_1},{tau_i},t_0, w, a, dt, T1, T2, bg, t) = 
%   bg + N/Q Sum[ q_i (sqrt(pi/2) w/tau_i exp(w^2/(2 tau_i^2) - t'/tau_i) Erfc((w^2- t' tau_i)/(sqrt(2)w tau_t)) +
%                 a P H(t')) (tau_i/(tau_i/T12 -1)/(tau_1/T1 -1)/T2 exp(-t''\tau_i) + 
%                                 1/(tau_1/T12 -1) exp(-t''/T12) - 
%                                 1/(tau_i/T1  -1) exp(-t''/T1) ) ]
%   t'  = t - t_0
%   t'' = t - t_0 - dt;
%   Q   = sqrt(2 pi w^2) + a T1(1+T1/T2)^{T2/T1}
%   P   = (1-T2/T1)(1+T1/T2)^{T2/T1}
%
%  (see: Walther, K.A. et al, Mol. BioSyst. (2011) doi:10.1039/c0mb00132e)

p = p(:);

if (nargin>3)&&(~isempty(para))
    para = para(:);    
    n = numel(para);
    if n > 6 
        p = [para; p];
    else
        p = [p(1:7-n); para; p(8-n:end)];
    end
end

nex = numel(p) - 7;
tau = p(8:end);

if length(t)<length(y) 
        
    c  = t;
    t  = y(:);
    t  = t - t(1);
        
    IRF = IRF_Fun(p(1:6), t);
    
    zz = zeros(numel(t),nex+2);
    
    zz(:,1) = ones(size(t));
    zz(:,2) = IRF;

    for k = 1:nex
        tmp = Convol(IRF, exp(-t./tau(k))./tau(k));
        zz(:,2+k) = tmp(1:numel(t));
    end;
    
    zz = zz./(ones(length(t),1)*sum(zz));
    
    for j=1:size(c,2)
        err(:,j) = zz*c(:,j);
    end

else

    [m, n] = size(y);
    if m<n 
        y   = y'; 
        [m, n] = size(y);
    end
    
	t = t(isfinite(sum(y,2)));
    t = t - t(1);
    y = y(isfinite(sum(y,2)),:);
    t = t(y>0);
    y = y(y>0);
        
    IRF = IRF_Fun(p(1:6), t);

    zz = zeros(numel(t),nex+2);
    
    zz(:,1) = ones(size(t));
    zz(:,2) = IRF;

    for k = 1:nex
        tmp = Convol(IRF, exp(-t./tau(k))./tau(k));
        zz(:,2+k) = tmp(1:numel(t));
    end
    
    zz = zz./(ones(length(t),1)*sum(zz));

    for j=1:n
        c(:,j) = lsqnonneg(zz,y(:,j)); %%
        z(:,j) = zz*c(:,j);
    end
  
%     if nargin>4 && ~isempty(pic)
%         if pic==1
%             plot(t, y, 'ob', t, z, 'r'); drawnow
%         else
%             semilogy(t, abs(y), 'o', t, abs(z)); drawnow
%         end
%     end

    err = sum((y-z).^2./abs(y));
end

