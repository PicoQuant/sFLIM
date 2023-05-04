function varargout = EFilter(sze, cutoffM, cutoffm, n, varargin)

% EFILTER constructs an elliptical lowpass filter Butterworth filter
%   E = EFILTER(Size, cutoffM, cutoffm, n) designs an Nth order elliptical
%   lowpass digital Butterworth filter where Size is a two element
%   [rows, cols] vector specifying the size of the filter to construct,
%   cutoffM and cutoffm, the cutoff freqency on the major and minor
%   axes are 0 < cutoff <= 1.
%
%   If E = EFilter(Size, cutoffM, cutoffm, n, alpha), where alpha is an angle
%   in radians, it will return and plot an elliptical filter rotated
%   counter-clockwise through alpha.
%
%   If E = EFilter(Size, cutoffM, cutoffm, n, alpha, xoff, yoff), where xoff
%   and yoff are offsets in the x and y direction, it will return and
%   plot an eliptical filter which is offset by the specified amount.
%   An offset of 0 corresponds to the center and an offset of 1
%   corresponds to the edge of the filter. A positive offset shifts the
%   filter in the positive direction.
%
%   Calling EFilter(...) without assigning the output variable
%   plots the 3D surface described by the function.

% Katie Streit   kstreit@rice.edu
% ELEC 301
% Rice University
%
% December 2001

% Much of this code was based on Peter Kovesi's  (pk@cs.uwa.edu.au)
% Matlab function for a lowpass Butterworth filter.

if nargin == 4
    alpha = 0;
    offx = 0;
    offy = 0;
elseif nargin == 5
    offx = 0;
    offy = 0;
    alpha = varargin{1};
elseif nargin == 7
    alpha = varargin{1};
    offx = varargin{2};
    offy = varargin{3};
else
    error('Invalid number of input arguments');
end

if nargout > 1
    error('Invalid number of output arguments');
end

if cutoffM < 0 || cutoffM > 1
    error('cutoffM frequency must be between 0 and 1');
end

if cutoffm < 0 || cutoffm > 1
    error('cutoffm frequency must be between 0 and 1');
end

if rem(n,1) ~= 0 || n < 1
    error('n must be an integer >= 1');
end

% extract the sizes from sze

rows = sze(1);
cols = sze(2);

% x and y matrices normalized to +/-.5 and an offset of offx or offy

x =  ((((ones(rows,1) * [1:cols])-offx*rows/2)  - (fix(cols/2)+1))/cols);
y =  ((([1:rows]' * ones(1,cols))-offy*rows/2) - (fix(rows/2)+1))/rows;

% apply a linear transformation to rotate through alpha. Note that it 
% uses negative alpha, which is caused by x and y being independent matrices.

x2 = (x*cos(alpha) - y*sin(-alpha));
y2 = (x*sin(-alpha) + y*cos(alpha));

% constructs an elliptical cone (defined by a and b) of height r on each at
% each elliptical ring. (r is effectively the "radius")
% r = sqrt(((x2/a).^2 + (y2/b).^2));

% Designs the filter
% f = 1./(1.0 + (r./cutoff).^(2*n));

a = cutoffM/2;
b = cutoffm/2;

f = 1./(1+((x2/(a)).^2 + (y2/(b)).^2).^n);

if nargout > 0
    varargout{1} = f;
else
    %Plots a normalized (+/- 1), interpolated 3D image of the filter
    surf([-1:2/(cols-1):1],[-1:2/(rows-1):1], f);
    shading interp;
    title('Elliptical Butterworth filter');
    xlabel('x');
    ylabel('y');
    zlabel('intensity');
    grid on;
end