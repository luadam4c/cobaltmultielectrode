function [maxtab, mintab]=peakdet(v, delta, x) 
%PEAKDET find peaks in a time series
% [maxtab, mintab]=peakdet(v, delta, x) finds the local maxiam
% and minima in the vector v.

maxtab = [];
mintab = [];

v = v(:); 

% Just in case this wasn�t a proper vector
if nargin < 3
    x = (1:length(v))';
else
    x = x(:) ;
    if length(v)~= length(x)
    error('Input vectors v and x must have same length');
    end
end

if (length(delta(:)))>1 
    error('Input argument DELTA must be a scalar');
end

if delta <= 0 
    error('Input argument DELTA must be positive');
end

mn = Inf; mx =-Inf; 
mnpos = NaN; mxpos = NaN;

lookformax = 1;

for i=1:length(v)
    this = v(i);
    if this > mx, mx = this; mxpos = x(i); end 
    if this < mn, mn = this; mnpos = x(i); end
    
    if lookformax 
            if this < mx-delta 
                maxtab = [maxtab ; mxpos mx];
                mn = this; mnpos = x(i);
                lookformax = 0;
            end
    else
        if this > mn+delta
            mintab = [mintab ; mnpos mn];
            mx = this; mxpos = x(i);
            lookformax = 1;
        end
    end
end