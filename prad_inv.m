
function [x,b, I0r]=prad_inv(X0,I,I0, KB, x_bc, b_bc)
%% PRAD_INV: solve prad inverse problem with boundary conditions
% 
% Given the proton mapping x' = x + b/KB,
% where b is the line-integrated magnetic field, and KB a constant,
% the detected fluence I is related to the input fluence
% I0(x) by  the Jacobion transformation,
%
% I(x') = I0(x) / |dx'/dx|
%
% The inverse problem solves for b(x) given I and I0 as data,
% such that the boundary conditions {x_bc, b_bc} are satisfied.
% I0 is renormalized to satisfy boundary conditions.
%
% [x,b,I0r] =  PRAD_INV (X, I, I0, KB, x_bc, b_bc)
%
% inputs:
% X0: spatial coordinate
% I: measured proton fluence at points X0
% I0: undisturbed proton fluence at points X0.
%       If I0 is scalar, indicates uniform through domain.
% nu = deflection parameter (constant)
% x_bc = [x1 x2] - boundary condition coordinates
% b_bc = [b1 b2] - boundary condition b field
%
% outputs:
% x = points in plasma plane where B was determined
% b = line-integrated field from solving inverse problem
% I0r = the renormalized fluence profile used
%
% The routine throws a warning if x' or x are generated which
% are out of the data domain of X.   Simple extrapolation is used 
% to allow moderately out of domain data in this case, but you are warned.
%
% See also prad_inv_I0, prad_fwd

if length(x_bc) ~= 2
    error('Expect x_bc = [x1 x2] two element vector')
end

if length(b_bc) ~= 2
    error('Expect b_bc = [b1 b2] two element vector')
end

% transpose to required orientation
if size(X0,1) == 1
    X0 = X0';
end

if size(I,1) == 1
    I = I';
end

if size(I0,1) == 1
    I0 = I0';
end

% commensurate size checks
if any( size(I) ~= size(X0) )
    error('prad_inv: Mismatch in size(X0) and size(I)')
end

if length(I0)>1 && any( size(I0) ~= size(X0) )
    error('prad_inv: Mismatch in size(X0) and size(I0), or I0 should be a scalar')
end


[x_bc,ind] = sort(x_bc);
b_bc = b_bc(ind);

x1 = x_bc(1);
x2 = x_bc(2);

b1 = b_bc(1);
b2 = b_bc(2);


xp1 = x1 + b1/KB;
xp2 = x2 + b2/KB;


% re-normalize I0 to fit boundary conditions

X0_new = unique(sort( [x1; x2; X0] ));
X0_new = X0_new ( X0_new >= x1 & X0_new <= x2);

% calculate integral I0 dx from x1 to x2
if length(I0)>1
    I0_new = interp1(X0, I0, X0_new, 'linear', 'extrap');
    I0_sum = trapz(X0_new, I0_new);
else
    % scalar I0, no interpolation needed
    I0_sum = abs(x2-x1) * I0;
end

% calculate integral I dx' from x'1 to x'2
Xp_new = unique(sort( [xp1; xp2; X0] ));
Xp_new = Xp_new ( Xp_new >= xp1 & Xp_new <= xp2);

I_new = interp1(X0, I, Xp_new, 'linear', 'extrap');
I_sum = trapz(Xp_new, I_new);

% renormalized I0 to achieve boundary conditions
I0r = I0 * (I_sum / I0_sum);

% now just call prad_inv_I0 as usual, starting from (x1,B1)
[x,b]=prad_inv_I0(X0,I,I0r, KB, x1, b1);
   
% return I0r, either as scalar, or interpolated to x
if length(I0)>1
   I0r = interp1(X0, I0r, x, 'linear' ,'extrap');
end

 
end
