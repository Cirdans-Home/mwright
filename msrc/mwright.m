function [w,N] = mwright(x,t,lambda,mu,varargin)
%MWRIGHT This routine computes the Wright function on the real line:
%:math:`t^{-\lambda} W_{\lambda,\mu}(-|x|t^{-\lambda})`
%
%  :param: x (may be a vector)
%  :param: t has to be a scalar > 0
%  :param: -1 < lambda < 0
%  :param: mu scalar
%  :param: Optional, the number of quadrature points N
%  :Returns: w and number of used quadrature points N

%% Check the inputs
if ~isreal(x)
    error('mwright:incorrectType',...
        'Error. \nInput x must be a real, not a %s.',class(n))
end
if ~isreal(t)
    error('mwright:incorrectType',...
        'Error. \nInput t must be a real, not a %s.',class(n))
end
if lambda <= -1
    error('mwright:incorrectValue',...
        'Error. \nInput lambda must be > -1, not %f.',lambda)
end
if t < 0
    error('mwright:incorrectValue',...
        'Error. \nInput t must be > 0, not %f.',lambda)
end

%% Determine the value of N
% We select it either from the theoretical analysis or let the user input 
% it, maybe the user knows what is doing.
if nargin == 4
    N = round(-3/(2*pi)*log(eps));
else
    N = varargin{1};
end
%% Build the integration contour
gamma = (pi*N)./(12*t);
z  = @(u) gamma*(1i*u + 1).^2;
zp = @(u) 2*1i*gamma*(1i*u + 1);
%% Determine the nodes
h = 3/N;
uk = -N:h:N;
%% Apply the quadrature rule
w = zeros(size(x));
zuk = z(uk);
for k=1:length(w)
    w(k) = h/(2*pi*1i)*sum(exp(zuk.*t).*zuk.^(-mu).*exp(-abs(x(k)).*zuk.^(-lambda)).*zp(uk));
end

%% Check that everything is real
w = real(w);

end

