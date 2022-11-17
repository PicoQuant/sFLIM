function [z, xv] = mHist(x,xv)

if size(x,1)==1 || size(x,2)==1
    x = x(:);
    x(~isfinite(x)) = [];
end

xdim = max([1 size(x,2)]);

if nargin>1 && ~isempty(xv)
    xmin = xv(1);
    xmax = xv(end);
    ind = x>xmax | x<xmin;
    x(ind) = [];
 
    if sum(diff(diff(xv)))==0
        dx = 0;
        if numel(xv)>1
            dx = xv(2)-xv(1);
        end
        if dx ~=0
            x = round((x-xmin)/dx)+1;
            xmax = round((xmax-xmin)/dx)+1;
        else
            x = round(x-xmin)+1;
            xmax = 0;            
        end
    else
        x = round(interp1(xv,1:length(xv),x));
        xmax = round(interp1(xv,1:length(xv),xmax));
    end
end

z = zeros(length(xv)*xdim,1);
if ~isempty(x)
    num = sort(x+xmax*ones(size(x,1),1)*(0:xdim-1));
    num = num(:);
    z(num) = 1;
    tmp = diff(diff([0; num; 0])==0);
    ind = (1:length(num))';
    z(num(tmp==1)) = z(num(tmp==1))-ind(tmp==1)+ind(tmp==-1);
    z = reshape(z,length(xv),xdim);
end;

if nargout==0
    bar(xv,z,1,'b'); 
    clear z xv
end
