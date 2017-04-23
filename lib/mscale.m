function [scaled] = mscale(data, nmin, nmax)
    cmax = max(data(:));
    cmin = min(data(:));
    scaled =((data-cmin)*(nmax-nmin))/(cmax-cmin) + nmax;
end