function FastFLIM(fname)

close all

tmp = strfind(fname,'\');
pathname = fname(1:tmp(end));
filename = fname(tmp(end)+1:end);

[head, tag, tcspc] = load_data([pathname filename]);

load Get_Params;

nx     = head.ImgHdr.PixX;
ny     = head.ImgHdr.PixY;
tx     = floor(nx/binning);
ty     = floor(ny/binning);
x      = head.ImgHdr.X0+(1:binning:nx)*head.ImgHdr.PixelSize;
y      = head.ImgHdr.Y0+(1:binning:ny)*head.ImgHdr.PixelSize;
nch     = size(tcspc,1);
% nbin    = size(tcspc,2);
num_PIE = size(tcspc,3);
tau     = head.tau;

load([pathname filename(1:end-4) '_Moments.mat'],'tav');

t_m = (tav(:,:,:,:,1));

clear tav;

if binning > 1
    tmp = tag(1:binning*tx,:,:,:);
    tmp = reshape(tmp, binning, tx, ny, nch, num_PIE);
    tmp = shiftdim(sum(tmp,1),1);
    tmp = tmp(:,1:binning*ty,:,:);
    tmp = reshape(tmp, tx, binning, ty, nch, num_PIE);
    tag = permute(sum(tmp,2), [1 3 4 5 2]);
end

ind = (sum(tag,3)) < Threshold;
t_m(ind) = 0; 

ind = (tag >= Threshold);
t_m = squeeze(sum(ind.*t_m.*tag,3)./sum(ind.*tag,3));
t_m(t_m<0) = 0;


intens = squeeze(sqrt(sum(ind.*tag,3)));
intens = intens./max(max(max(intens)));

%save('updatefflim', 't_m', 'ntau', 'nbin', 'stim');

clear S1 S2 S3;

xl = [399 440 485 540 580 610 630 699];
yl = [[0 0 0]; [0 0 1]; [0 1 1]; [0 1 0]; [1 1 0]; [1 0.65 0]; [1 0 0]; [1 0 1]];
lambda   = 399:3:699;
spectrum = interp1(xl, yl, lambda);

p1 = 0.06;
p2 = 0.555;
w  = 0.435;

for pulse = 1:num_PIE

    name = sprintf('FastFLIM %d', pulse);

    figure(figFLIM+pulse);
    set(gcf,'Name',name,'NumberTitle','off');

    S1 = subplot('Position', [p1 p2 w w]);

    tmp = repmat(intens(:,:,pulse)./max(max(intens(:,:,pulse))), [1 1 3]);

    lims = [tau_min tau_max];
    
    val = reshape(t_m(:,:,pulse),size(t_m,1)*size(t_m,2),1);

    val(isnan(val))  = 0;
    val(val<lims(1)) = 0;
    val(val>lims(2)) = lims(2);

    k      = 2 + round((val-lims(1))./(lims(2)-lims(1)).*(size(spectrum,1)-2));
    k(k<2) = 1;

    im = spectrum(k,:);
    im = reshape(im,size(tmp,1),size(tmp,2),3);
    im = im.*tmp;

    image(x,y,im)

    colormap(spectrum)

    set(S1,'DataAspectRatio', [1,1,1], ...
        'PlotBoxAspectRatio',[1 1 1], ...
        'XDir','normal', ...
        'YDir','reverse', ...
        'FontSize',9,...
        'CLim', lims, ...
        'color', [0 0 0]);

    xlabel('x (µm)','FontWeight','bold','FontSize',10,'FontAngle','italic');
    ylabel('y (µm)','FontWeight','bold','FontSize',10,'FontAngle','italic');


    tst = 2 + round(((lims(1):.5:lims(2))-lims(1))./(lims(2)-lims(1)).*(size(spectrum,1)-3));
    n= 1;
    for ttst = lims(1):0.5:lims(2)
        stst(n) = cellstr(sprintf('%4.1f',ttst));
        n = n+1;
    end
    cca = colorbar;
    cca.Label.String = 'lifetime / ns';
    %colorbar('CLim', lims, 'YTickMode','manual','XColor',[0 0 0],'YColor',[0 0 0],'Box','off','YTick', tst, ...
    %    'FontSize', 9,'TickDir','out', 'YTickLabel',stst);

%    save for fastflimgraph
%    save([name '_graph1'], 'x', 'y');

    S2 = subplot('Position', [p1 p1 w w]);

    W = tag(:,:,:,pulse);
    W(W<Threshold) = 0;

    K = (1:size(tag,3));
    K = repmat(reshape(K,[1 1 numel(K)]),[size(tag,1) size(tag,2) 1]);

    Z = W.*K;
    Z = lam_start + lam_step.*(sum(Z,3)./sum(W,3)-1.0);

%    lims = [399 699];
    lims = [clam_min clam_max];

    val = Z;
    val = reshape(val,size(val,1)*size(val,2),1);

    val(isnan(val))  = 0;
    val(val<lims(1)) = 0;
    val(val>lims(2)) = lims(2);

    k      = 2 + round((val-lims(1))./(lims(2)-lims(1)).*(size(spectrum,1)-3));
    k(k<2) = 1;

    im = spectrum(k,:);
    im = reshape(im,size(tmp,1),size(tmp,2),3);
    im = im.*tmp;

    image(x,y,im)

    set(S2,'DataAspectRatio', [1,1,1], ...
        'PlotBoxAspectRatio',[1 1 1], ...
        'XDir','normal', ...
        'YDir','reverse', ...
        'FontSize',9,...
        'CLim', lims, ...
        'color', [0 0 0]);

    xlabel('x (µm)','FontWeight','bold','FontSize',10,'FontAngle','italic');
    ylabel('y (µm)','FontWeight','bold','FontSize',10,'FontAngle','italic');

    tst = 2 + round(((440:40:680)-lims(1))./(lims(2)-lims(1)).*(size(spectrum,1)-3));
    ccb = colorbar;
    ccb.Label.String = 'wavelengths / nm';
    %colorbar('CLim', lims, 'XColor',[0 0 0],'YColor',[0 0 0],'Box','off','YTick', tst, ...
    %    'FontSize',9, 'YTickLabel',{'440', '480', '520', '560', '600', '640', '680'},'TickDir','out');

%    save for fastflimgraph
%    save([name '_graph2'], 'x', 'y', 'im', 'spectrum', 'lims', 'tst')

    S3 = subplot('Position', [p2 p1 w w]);

    s   = squeeze(sum(sum(W,2),1));
    
    %lam = (lam_start+lam_step.*(-1:numel(s)));
    lam = (lam_start+lam_step.*(0:(numel(s)-1)));
     
     disp(lam_step);
    ts  = clam_min:clam_max;
    if max(s)>0 
        %s = [0; s(:); 0]./max(s);
        s = [s(:)]./max(s);
        %cs  = interp1(lam, s, ts,'cubic');
        %cs(ts<lam_start-lam_step) = 0;
        %cs(ts>lam_start+lam_step.*nch) = 0;
    else
        cs = zeros(size(ts));
    end
   
    %plot(ts, cs, '-','Color',[1 0 0]);
     
    disp(s);
    disp(lam);
   
    plot(lam, s,'b--o');

    set(S3,'FontSize',9,'color', [204 204 204]./255);

    %    axis([10*floor(lam(1)/10) 10*ceil(lam(end)/10) 0 1.1])
    axis([clam_min clam_max 0 1.1])

    xlabel('wavelength of \lambda_0 (nm)','FontWeight','bold','FontSize',10,'FontAngle','italic');
    ylabel('rel. spectral irradiance','FontWeight','bold','FontSize',10,'FontAngle','italic');

%    save for fastflimgraph
%    save([name '_graph3'],'ts', 'cs');

    S4 = subplot('Position', [p2 p2 w w]);

    if nch == 1
        k_ind = floor(1+(numel(lambda)-1));
    else
        k_ind = floor(1+((1:nch)-1).*(numel(lambda)-1)./(nch-1));
    end

    m = max(max(tcspc(:,:,pulse)));
    tcspc(tcspc<=1) = 1e-4;
    semilogy(tau,tcspc(1,:,pulse)./m,'-','Color',spectrum(k_ind(1),:),'LineWidth',2);
    hold on;
    for n = 2:size(tcspc,1)
        semilogy(tau,tcspc(n,:,pulse)./m,'-','Color',spectrum(k_ind(n),:),'LineWidth',2);
    end

    hold off

    set(S4,'FontSize',9,'color', [204 204 204]./255);

    axis([0 floor(tau(end)) 1e-4 1.1])

    xlabel('time (ns)','FontWeight','bold','FontSize',10,'FontAngle','italic');
    ylabel('detection frequency','FontWeight','bold','FontSize',10,'FontAngle','italic');

%    save for fastflimgraph
%    save([name '_graph4'],'tau', 'tcspc', 'm', 'k_ind')

end