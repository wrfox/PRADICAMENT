function [X,I]=prad_fwd(x0,b,KB,I0)
%% PRAD_FWD:  model 1-D proton fluence given magnetic field
% 
% Given the proton mapping x' = x + b/KB,
% where b is the line-integrated magnetic field, and KB a constant.
% This routine obtains the forward proton fluence I using histogramming
%
% [X,I] =  PRAD_FWD (x0, b, KB, I0)
%
% inputs:
% x0: spatial coordinate where B defined
% b: deflecting line-integrated B field
% KB = deflection parameter (constant)
% I0: undisturbed proton fluence, defined at points x0
%
% outputs:
% X = points in plasma plane where I is determined
% I = proton fluence


% if I0 a singleton, assume constant
if (length(I0) == 1)
    I0 = I0 + 0*x0;
end

if any( size(x0) ~= size(b) )
    error('prad_inv: Mismatch in size(x0) and size(b)')
end



% Generate non-uniform launching points to follow I0

scale=10000;                % just a large number
nx0=length(x0)-1;
np=nx0*scale;

x0_samp = linspace(x0(1),x0(end),np);  %initial sampling of launching points
I0_samp = interp1(x0,I0,x0_samp);       % upsample the I0 function

I0_sum = cumsum(I0_samp);         % integrate
I0_sum = I0_sum / I0_sum(end);    % normalize I0_sum(end) = 1

ns = (1:np)/np;
x0s = interp1(I0_sum,x0_samp,ns);   % sample the launching points according to the BG


%forward scattering
bs = interp1(x0,b,x0s); %
xps = x0s + bs / KB;

% set edges according to final scattering points
dx = mean(diff(x0));
edges = [min(xps):dx:max(xps)];

I = histcounts(xps,edges) * trapz(x0,I0) / dx / np;       % histogram onto the original points

size(I);
size(x0);

% set return data points at center of bins
X = edges(1:end-1)+dx/2;


end