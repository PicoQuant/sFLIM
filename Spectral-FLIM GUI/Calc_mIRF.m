function IRF = Calc_mIRF(head, tcspc)

maxres = max([head.Resolution]);
Resolution = max([maxres 0.032]);
Pulse      = 1e9/head.SyncRate;

tau = Resolution.*((1:size(tcspc,2))-0.5)';
IRF = zeros(size(tcspc));
nex = 2;

t0 = 1.5;
w1 = 0.5^2;  % width of excitation peak
T1 = 0.10;    % time constant of rise term in IRF
T2 = 0.10;    % time constant of decay term in IRF
a  = 0.1;     % rel. amplitude of second peak in IRF
b  = 0.1;     % rel. amplitude of second peak in IRF
dt = 0.0;     % time shift of second peak in IRF

for PIE = 1:size(tcspc,3)

    p   = [t0     w1     T1   T2   a  b   dt    1 3]';
    pl  = [t0-1.0 1e-3 1e-4 1e-4   0  0  -0.3   0.5*ones(1, nex)]';
    pu  = [t0+1.0 1e-1  1    1      1  1   0.5  5*ones(1, nex)]';
    
    tc = squeeze(sum(tcspc(:,:,PIE),2));
    [tmp, ord] = sort(tc,'descend');

    ch = 1;
    ind = ord(ch);
    y = squeeze(tcspc(ind,:,PIE));

    err = zeros(1,10);

    for casc=1:10
        [ts, s] = min(err);
        r0 = p(:, s);
        for sub=1:10
            rf = r0.*[2.^(1.1*(rand(size(r0))-0.5)./casc)];  % randomize start values
            rf = max([rf pl],[],2);
            rf = min([rf pu],[],2);
            p(:,sub) = Simplex('TCSPC_Fun',rf,pl,pu,[],[],tau, y);
            err(sub) = TCSPC_Fun(p(:,sub), tau, y);
        end
    end

    err1 = min(err);
    p1   = mean(p(:,err==err1),2);
    [tmp, c1, bla1, tmp1] = TCSPC_Fun(p1, tau, y);

    IRF(ind,:,PIE) = IRF_Fun(p1(1:7),tau);
    IRF(ind,:,PIE) = IRF(ind,:,PIE)./max(IRF(ind,:,PIE));

    para = p1(2:7);
    p    = [p1(1);  p1(8:end)];
    pl   = [  0.5; 0.5*ones(nex, 1)];
    pu   = [  3.0;  10*ones(nex, 1)];

    for ch = 2:size(tcspc,1)
        ind = ord(ch);
        y = squeeze(tcspc(ind,:,PIE));
        
        err = zeros(1,10);

        for casc=1:10
            [ts, s] = min(err);
            r0 = p(:, s);
            for sub=1:10
                rf = r0.*[2.^(1.05*(rand(size(r0))-0.5)./casc)];  % randomize start values
                rf = max([rf pl],[],2);
                rf = min([rf pu],[],2);
                p(:,sub) = Simplex('TCSPC_Fun',rf,pl,pu,[],[],tau, y, para);
                err(sub) = TCSPC_Fun(p(:,sub), tau, y, para);
            end
        end

        err1 = min(err);
        p1   = mean(p(:,err==err1),2);
        [tmp, c1, bla1, tmp1] = TCSPC_Fun(p1, tau, y, para);
 
        IRF(ind,:,PIE) = IRF_Fun([p1(1); para; p1(2:end)],tau);
        IRF(ind,:,PIE) = IRF(ind,:,PIE)./max(IRF(ind,:,PIE));
    end;    
end

IRF(IRF<1e-4) = 0;
IRF = IRF./repmat(sum(IRF,2),[1 size(IRF,2) 1]);