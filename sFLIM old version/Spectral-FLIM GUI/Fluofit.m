function [c, offset, A, tau, dc, dtau, irs, zz, t, chi] = Fluofit(irf, y, p, dt, tau, lim, init)
% The function FLUOFIT performs a fit of a multi-exponential decay curve.
% It is called by: 
% [c, offset, A, tau, dc, doffset, dtau, irs, z, t, chi] = fluofit(irf, y, p, dt, tau, limits, init).
% The function arguments are:
% irf 	= 	Instrumental Response Function
% y 	= 	Fluorescence decay data
% p 	= 	Time between laser exciation pulses (in nanoseconds)
% dt 	= 	Time width of one TCSPC channel (in nanoseconds)
% tau 	= 	Initial guess times
% lim   = 	limits for the lifetimes guess times
% init	=	Whether to use a initial guess routine or not 
%
% The return parameters are:
% c	=	Color Shift (time shift of the IRF with respect to the fluorescence curve)
% offset	=	Offset
% A	    =   Amplitudes of the different decay components
% tau	=	Decay times of the different decay components
% dc	=	Color shift error
% doffset	= 	Offset error
% dtau	=	Decay times error
% irs	=	IRF, shifted by the value of the colorshift
% zz	    Fitted fluorecence component curves
% t     =   time axis
% chi   =   chi2 value
% 
% The program needs the following m-files: simplex.m, lsfit.m, mlfit.m, and convol.m.
% (c) 1996 Jörg Enderlein


fitfun = 'lsfit';

close all
irf = irf(:);
offset = 0;
y = y(:);
n = length(irf); 
if nargin>6
    if isempty(init)
        init = 1;
    end
elseif nargin>4 
    init = 0;
else
    init = 1;
end

if init>0 
    [cx, tau, c, c] = DistFluofit(irf, y, p, dt, [-3 3]);    
    cx = cx(:)';
    tmp = cx>0;
    t = 1:length(tmp);
    t1 = t(tmp(2:end)>tmp(1:end-1)) + 1;
    t2 = t(tmp(1:end-1)>tmp(2:end));
    if length(t1)==length(t2)+1 
        t1(end)=[]; 
    end
    if length(t2)==length(t1)+1 
        t2(1)=[]; 
    end
    if t1(1)>t2(1)
        t1(end)=[]; 
        t2(1)=[];
    end
    tmp = [];
    for j=1:length(t1)
        tmp = [tmp cx(t1(j):t2(j))*tau(t1(j):t2(j))/sum(cx(t1(j):t2(j)))];
    end
    tau = tmp;
else
    c = 0;
end

if (nargin<6)||isempty(lim)
    lim = [zeros(1,length(tau)) 100.*ones(1,length(tau))];
end;

p = p/dt;
tp = (1:p)';
tau = tau(:)'/dt; 
lim_min = lim(1:numel(tau))./dt;
lim_max = lim(numel(tau)+1:end)./dt;
t = 1:length(y);
m = length(tau);
x = exp(-(tp-1)*(1./tau))*diag(1./(1-exp(-p./tau)));
irs = (1-c+floor(c))*irf(rem(rem(t-floor(c)-1, n)+n,n)+1) + (c-floor(c))*irf(rem(rem(t-ceil(c)-1, n)+n,n)+1);
z = convol(irs, x);
z = [ones(size(z,1),1) z];
%A = z\y;
A = lsqnonneg(z,y);
z = z*A;
close all
if init<2
    disp('Fit =                Parameters =');
    param = [c; tau'];
    % Decay times and Offset are assumed to be positive.
    paramin = [-1/dt lim_min];
    paramax = [ 1/dt lim_max];
    [param, dparam] = Simplex(fitfun, param, paramin, paramax, [], [], irf(:), y(:), p);
    c = param(1);
    dc = dparam(1);
    tau = param(2:length(param))';
    dtau = dparam(2:length(param));
    x = exp(-(tp-1)*(1./tau))*diag(1./(1-exp(-p./tau)));
    irs = (1-c+floor(c))*irf(rem(rem(t-floor(c)-1, n)+n,n)+1) + (c-floor(c))*irf(rem(rem(t-ceil(c)-1, n)+n,n)+1);
    z = convol(irs, x);
    z = [ones(size(z,1),1) z];
    z = z./(ones(n,1)*sum(z));
    %A = z\y;
    A = lsqnonneg(z,y);
    zz = z.*(ones(size(z,1),1)*A');
    z = z*A;
    dtau = dtau;
    dc = dt*dc;
else
    dtau = 0;
    dc = 0;
end
chi = sum((y-z).^2./abs(z))/(n-m);
t = dt*t;
tau = dt*tau';
c = dt*c;
offset = zz(1,1); 
A(1) = [];
if 1
	hold off
    subplot('position',[0.1 0.4 0.8 0.5])
	plot(t,log10(y),t,log10(irs),t,log10(z));
	v = axis;
	v(1) = min(t);
	v(2) = max(t);
	axis(v);
	xlabel('Time / ns');
	ylabel('Log Count');
	s = sprintf('COF = %3.3f   %3.3f', c, offset);
	text(max(t)/2,v(4)-0.05*(v(4)-v(3)),s);
	s = ['AMP = '];
	for i=1:length(A)
		s = [s sprintf('%1.3f',A(i)/sum(A)) '   '];
	end
	text(max(t)/2,v(4)-0.12*(v(4)-v(3)),s);
	s = ['TAU = '];
	for i=1:length(tau)
		s = [s sprintf('%3.3f',tau(i)) '   '];
	end
	text(max(t)/2,v(4)-0.19*(v(4)-v(3)),s);
    subplot('position',[0.1 0.1 0.8 0.2])
	plot(t,(y-z)./sqrt(abs(z)));
	v = axis;
	v(1) = min(t);
	v(2) = max(t);

    axis(v);
	xlabel('Time / ns');
	ylabel('Residue');
	s = sprintf('%3.3f', chi);
	text(max(t)/2,v(4)-0.1*(v(4)-v(3)),['\chi^2 = ' s]);
    set(gcf,'units','normalized','position',[0.01 0.05 0.98 0.83])
end
