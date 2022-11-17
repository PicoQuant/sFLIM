function [err, A, z] = dlsfit(param, t, irf, y, p, pic)

%	DLSFIT(param, t, irf, y, p, pic) returns the Least-Squares deviation between the data y 
%	and the computed values. 
%	LSFIT assumes a function of the form:
%
%	  y =  yoffset + A(1)*convol(irf,exp(-t/tau(1)/(1-exp(-p/tau(1)))) + ...
%
%	param(1) is the color shift value between irf and y.
%	param(2) is the irf offset.
%	param(3:...) are the decay times.
%	irf is the measured Instrumental Response Function.
%	y is the measured fluorescence decay curve.
%	p is the time between to laser excitations (in number of TCSPC channels).

if (nargin==6)&&~isempty(pic)
    pic = true;
else
    pic = false;
end

c   = param(1);
tau = param(2:length(param)); 
tau = tau(:)';

n = length(irf);
t = t(:);
tp = (0:p-1)';

ind = 1+intersect(round(t),round(tp));

irf1 = interp1(t,irf,tp,'pchip');

x = exp(-tp*(1./tau))*diag(1./(1-exp(-p./tau)));

irs = (1-c+floor(c))*irf1(rem(rem(tp-floor(c), p)+p,p)+1) + ...
      (  c-floor(c))*irf1(rem(rem(tp-ceil(c), p)+p,p)+1);

z = convol(irs, x);
z = z(ind,:);

z = [ones(size(z,1),1) z];
A = z\y;  %A = lsqnonneg(z,y);
z = z*A;

if pic
    semilogy(t, irf/max(irf)*max(y), t, y, 'bo', t, z);
	drawnow
end
err = sum((z-y).^2./abs(z))/(n-length(tau));

