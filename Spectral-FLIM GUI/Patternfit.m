function [c, offset, A, tau, dc, dtau, irs, zz, t, chi] = Patternfit(irf, y, p, dt, tau, lim)
% The function PATTERNFIT performs a fit of a multi-exponential decay curve.
% It is called by:
% [c, offset, A, tau, dc, doffset, dtau, irs, z, t, chi] = Patternfit(irf, y, p, dt, tau, limits, init).
% The function arguments are:
% irf 	= 	Instrumental Response Function
% y 	= 	Fluorescence decay data
% p 	= 	Time between laser exciation pulses (in nanoseconds)
% dt 	= 	Time width of one TCSPC channel (in nanoseconds)
% tau 	= 	Initial guess times
% lim   = 	limits for the lifetimes guess times
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

w = round(diff([0 dt])./(dt(2)-dt(1)));
t1 = find(diff(w));
t2 = find(diff(t1)==1);
w(t1(t2)+1)= w(t1(t2)+2);

dt = dt(:)-dt(1);

if (nargin<6)||isempty(lim)
    lim = [zeros(1,length(tau)) 6.*ones(1,length(tau))];
end;

timeunit = min(diff(dt));

p   = round(p/timeunit);
tp  = (dt./timeunit);
tpp = (0:p-1)';
ind = 1+intersect(round(tpp),round(tp));

tau = tau(:)'/timeunit;
lim_min = lim(1:numel(tau))./timeunit;
lim_max = lim(numel(tau)+1:end)./timeunit;

[err, A, z] = dlsfit([0 tau], tp, irf, y, p);

close all

param = [0; tau'];
paramin = [-1/timeunit lim_min];
paramax = [ 1/timeunit lim_max];

tmp = param;
k   = numel(param);
err = 0;

for casc=1:10
    [ts, s] = min(err);
    r0 = tmp(:, s);
    dlsfit(r0, tp, irf, y, p, 1);
    drawnow;
    for sub=1:10
        rf = r0.*[(2.^(1.2*(rand(k,1)-0.5)./casc))];  % randomize start values
        tmp(:,sub) = Simplex('dlsfit',rf,paramin,paramax,[],[], tp, irf, y, p);
        err(sub)   = dlsfit(tmp(:,sub), tp, irf, y, p);
    end
end

param = mean(tmp(:,err==min(err)),2);
[err, A, z] = dlsfit(param, tp, irf, y, p, 1);

[param, dparam] = Simplex('dlsfit', param, paramin, paramax, [], [], tp, irf, y, p, 1);
c = param(1);
dc = dparam(1);
tau = param(2:length(param))';
dtau = dparam(2:length(param));

chi = sum((y1-z).^2./abs(z))/(n-m);
t = timeunit*t;
tau = timeunit*tau';
c = timeunit*c;
offset = zz(1,1);
A(1) = [];

hold off
subplot('position',[0.1 0.4 0.8 0.5])
plot(t,log10(y1),t,log10(irs),t,log10(z));
v = axis;
v(1) = min(t);
v(2) = max(t);
axis(v);
xlabel('Time in ns');
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
plot(t,(y1-z)./sqrt(abs(z)));
v = axis;
v(1) = min(t);
v(2) = max(t);

axis(v);
xlabel('Time in ns');
ylabel('Residue');
s = sprintf('%3.3f', chi);
text(max(t)/2,v(4)-0.1*(v(4)-v(3)),['\chi^2 = ' s]);
set(gcf,'units','normalized','position',[0.01 0.05 0.98 0.83])
