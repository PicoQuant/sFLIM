clear all
close all

load('F:\Messungen\SFLIM\data\FRET_GFP_and_mRFP_DATA.mat')
load('F:\Messungen\SFLIM\data\Donor.pat','-mat')
Donor = pat;

load('F:\Messungen\SFLIM\data\FRET.pat','-mat')
Fret  = pat;

tau  = head.tau;
w    = head.tauw;
binw = floor(w./w(1));
tr = 1e9/head.SyncRate;
Resolution = 0.032;

decay = [];
fret  = [];
irf   = [];
for k = 1:numel(w)
    irf   = [irf   repmat(IRF(k)/w(k),[1 binw(k)])];
    decay = [decay repmat(Donor(k)/w(k),[1 binw(k)])];
    fret  = [fret  repmat(Fret(k)/w(k),[1 binw(k)])];
end

tau_n = (1:numel(decay)).*Resolution;
tau_n = tau_n+(tau(1)-tau_n(1));

semilogy(tau_n,decay,'o',tau_n,irf,'x')

[c, offset, A, tau, dc, dtau, irs, zz, t, chi] = Fluofit(irf, decay, tr, Resolution, [2 3], [0.5 0.5 10 10], 0);

close;

ind = find(A==max(A));
t_D = tau(ind);

t  = 1:numel(fret);
p  = tr/Resolution;
tp = (1:p)';
n  = numel(irf);
irs = (1-c+floor(c))*irf(rem(rem(t-floor(c)-1, n)+n,n)+1) + (c-floor(c))*irf(rem(rem(t-ceil(c)-1, n)+n,n)+1);

E = (0:99)./100;
for k = 1:numel(E)
    
    t_F = t_D.*(1-E(k));
    
    ttau = [t_F tau']/Resolution;
    
    m = numel(ttau);
    
    x = exp(-(tp-1)*(1./ttau))*diag(1./(1-exp(-p./ttau)));
    z = convol(irs(:), x);
    z = [ones(size(z,1),1) z];
    A = lsqnonneg(z,fret(:));
    z = z*A;
    chi(k) = sum((fret(:)-z).^2./abs(fret(:)))/(n-m);
end

[m, ind] = min(chi);

t_F = t_D.*(1-E(ind));

ttau = [t_F tau']/Resolution;

m = numel(tau);

x = exp(-(tp-1)*(1./ttau))*diag(1./(1-exp(-p./ttau)));
z = convol(irs(:), x);
z = [ones(size(z,1),1) z];
A = z\fret(:); % lsqnonneg(z,fret(:));
z = z*A;

semilogy(tau_n,fret,'x',tau_n,z)
