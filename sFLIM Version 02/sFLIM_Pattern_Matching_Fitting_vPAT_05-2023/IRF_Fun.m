function [z] = IRF_Fun(p, t, pic)

%  Computes the following model IRF 
%
%  IRF(t_0, w, a, dt, T1, T2, t) = 
%           1/Q [ exp(- t'^2/(2 w^2) +
%                 a  H(t') exp(- t'/T1) +
%                 b (1-H(t')) exp( t'/T2)) ] 
%
%   t'  = t - t_0

t = t(:);
p = p';

t_0 = p(1);
w1  = p(2);
T1  = p(3);
T2  = p(4);
a   = p(5);
b   = p(6);

t1 = t-t_0;

H  = ones(size(t));
H(t<(t_0)) = 0;

ind = [false(size(t1)) H==0 H==1];

IRF =  [exp(- t1.^2/(2*w1)) exp(-t1./T1) exp(t1./T2)];
IRF(ind) = 0;
IRF =  IRF./(ones(size(t1))*sum(IRF,1));
IRF =  (ones(size(t1))*[1 a b]).*IRF;

z   =  sum(IRF,2)./sum(sum(IRF));

if nargin>2 && ~isempty(pic)
    if pic==1
        plot(t, z, 'r'); drawnow
    else
        semilogy(t, z); drawnow
    end
end

