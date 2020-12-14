function [ibin] = circ_bini(x,nbin,pbin)

if nargin < 3
    pbin = 1/nbin;
end
if nargin < 2
    error('Missing input arguments!');
end

if ~isvector(x)
    error('Wrong input size!');
end
quantbeg  = mod(linspace(0-pbin/2,1-1/nbin-pbin/2,nbin),1);
quantend  = mod(quantbeg+pbin,1);
xbinbeg   = quantile(x(:),quantbeg);
xbinend   = quantile(x(:),quantend);

ibin = false(numel(x),nbin);
for i = 1:nbin
    if quantbeg(i)<quantend(i) %no wrap
        ibin(:,i) = x(:) >= xbinbeg(i) & x(:) <= xbinend(i);
    else %wrap
        ibin(:,i) = x(:) >= xbinbeg(i) | x(:) <= xbinend(i);
    end
end
end