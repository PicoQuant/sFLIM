function cscatter(h,X,Y, varargin)
% CSCATTER creates a scatter plot coloured by density.
%
%   CSCATTER(X,Y) creates a scatterplot of X and Y at the locations
%   specified by the vectors X and Y (which must be the same size), colored
%   by the density of the points.
%
%   CSCATTER(...,'MARKER',M) allows you to set the marker for the
%   scatter plot. Default is 's', square.
%
%   CSCATTER(...,'MSIZE',MS) allows you to set the marker size for the
%   scatter plot. Default is 10.
%
%   CSCATTER(...,'FILLED',false) sets the markers in the scatter plot to be
%   outline. The default is to use filled markers.
%
%   CSCATTER(...,'BINS',[NX,NY]) allows you to set the number of bins used
%   for the 2D histogram used to estimate the density. The default is to
%   use the number of unique values in X and Y up to a maximum of 200.
%
%   CSCATTER(...,'SMOOTHING',LAMBDA) allows you to set the smoothing factor
%   used by the density estimator. The default value is 20 which roughly
%   means that the smoothing is over 20 bins around a given point.
%
%   CSCATTER(...,'LOGY',true) uses a log scale for the yaxis.
%
% Reference:
% Paul H. C. Eilers and Jelle J. Goeman
% Enhancing scatterplots with smoothed densities
% Bioinformatics, Mar 2004; 20: 623 - 628.

lambda = [];
nbins  = [];
msize  = 10;
marker = 'diamond';
logy   = false;
filled = true;

if nargin > 3
    if rem(nargin,2) == 0
        error('Incorrect number of arguments to %s.',mfilename);
    end
    okargs = {'smoothing','bins','logy','marker','msize','filled'};
    for j=1:2:nargin-3
        pname = varargin{j};
        pval = varargin{j+1};
        k = strmatch(lower(pname), okargs); %#ok
        if isempty(k)
            error('Unknown parameter name: %s.',pname);
        elseif length(k)>1
            error('Ambiguous parameter name: %s.',pname);
        else
            switch(k)
                case 1  % smoothing factor
                    if isnumeric(pval)
                        lambda = pval;
                    else
                        error('Invalid smoothing parameter.');
                    end
                case 2
                    if isscalar(pval)
                        nbins = [pval pval];
                    else
                        nbins = pval;
                    end
                case 3
                    logy = pval;
                    Y = log10(Y);
                case 4
                    marker = pval;
                case 5
                    msize = pval;
                case 6
                    filled = pval;
            end
        end
    end
end

ind = ~isnan(X);
X = X(ind);
Y = Y(ind);
ind = ~isnan(Y);
X = X(ind);
Y = Y(ind);

if numel(X)>10
    minx = min(X,[],1);
    maxx = max(X,[],1);
    miny = min(Y,[],1);
    maxy = max(Y,[],1);

    if isempty(nbins)
        nbins = [min(numel(unique(X)),200) ,min(numel(unique(Y)),200) ];
    end

    if isempty(lambda)
        lambda = 20;
    end

    edges1 = linspace(minx, maxx, nbins(1)+1);
    edges1 = [-Inf edges1(2:end-1) Inf];
    edges2 = linspace(miny, maxy, nbins(2)+1);
    edges2 = [-Inf edges2(2:end-1) Inf];

    n = size(X,1);
    bin = zeros(n,2);
    % Reverse the columns to put the first column of X along the horizontal
    % axis, the second along the vertical.
    [dum,bin(:,2)] = histc(X,edges1);
    [dum,bin(:,1)] = histc(Y,edges2);
    H = accumarray(bin,1,nbins([2 1])) ./ n;
    G = smooth1D(H,nbins(2)/lambda);
    F = smooth1D(G',nbins(1)/lambda)';

    if logy
        Y = 10.^Y;
    end

    if numel(F)>0
        F = F./max(F(:));
        ind = sub2ind(size(F),bin(:,1),bin(:,2));
        col = F(ind);
        if filled
            h = scatter(h,X,Y,msize,col,marker,'filled');
        else
            h = scatter(h,X,Y,msize,col,marker);
        end
    else
        cla;
    end
else
    if filled
        h = scatter(h,X,Y,msize,[0.5 0.3 0.5],marker,'filled');
    else
        h = scatter(h,X,Y,msize,[0.5 0.3 0.5],marker);
    end
end


%set(h,'yscale','log');
%set(h,'xscale','log');

% if nargout > 0
%     hAxes = get(h,'parent');
% end
