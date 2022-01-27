function [w,N,bound] = mwright(x,t,lambda,mu,varargin)
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
        'Error. \nInput x must be a real, not a %s.',class(x))
end
if ~isreal(t)
    error('mwright:incorrectType',...
        'Error. \nInput t must be a real, not a %s.',class(t))
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
l = -log(eps);
ltol = -log(1e-15);
if nargin == 4
    if real(mu) < 2
        N = ceil(sqrt(2*l*ltol)/pi);
        h = (4*l)/(pi*N^2);
        gamma = (pi^2*N^2)/(16*t*l);
        bound = exp(-2*pi*1/h);
    elseif real(mu) == 2
        Nfun = @(c,mu) (sqrt(l*ltol)/pi)*...
            sqrt( 1 + (1./c).*(1 + 1/ltol*log(-log( (l-ltol)*(1-c).^2))));
        c  = fminbnd( @(c) Nfun(c,mu), 1-1/(l-ltol)+eps, 1);
        N  = ceil(Nfun(c,mu));
        xi = 2/( 1 + 1/ltol*log(-log( (l-ltol)*(1-c).^2)));
        h  = ((2+xi*c)*l)/(pi*N^2);
        gamma = (pi^2*N^2)/((2+xi*c)^2*t*l);
        bound = log(gamma*t*(1-c)^2)*exp(-2*pi*c/h);
    else
        Nfun = @(c,mu) (sqrt(l*ltol)/pi)*...
            sqrt( 1 + (1./c).*(1 + ((2 - real(mu))./ltol).*log(1-c)));
        c  = fminbnd( @(c) Nfun(c,mu),0,1);
        N  = ceil(Nfun(c,mu));
        xi = 2/(1 + ((2-real(mu))./ltol)*log(1-c));
        h  = ((2+xi*c)*l)/(pi*N^2);
        gamma = (pi^2*N^2)/((2+xi*c)^2*t*l);
        bound = (1-c)^(2-real(mu))*exp(-2*pi*c/h);
    end
else
    N = varargin{1};
    h = (4*l)/(pi*N^2);
    gamma = (pi^2*N^2)/(16*t*l);
    bound = exp(-2*pi*1/h);
end
%% Build the integration contour
z  = @(u) gamma*(1i*u + 1).^2;
zp = @(u) 2*1i*gamma*(1i*u + 1);
%% Determine the nodes
uk = (-N:N)*h;
%% Apply the quadrature rule
w = zeros(size(x));
zuk = z(uk);
for k=1:length(w)
    w(k) = (h/(2*pi*1i))*KahanBabushkaKleinSum(zuk.^(-mu).*exp( (zuk.*t) - (abs(x(k)).*zuk.^(-lambda))).*zp(uk));
end

%% Check that everything is real
if isreal(mu)
    w = real(w);
end


end

function sum = KahanBabushkaKleinSum(input)
%%KAHANBABUSHKAKLEINSUM Implementation of the sum algorithm from
%
% A., Klein (2006). "A generalized Kahan–Babuška-Summation-Algorithm".
% Computing. Springer-Verlag. 76 (3–4): 279–293.
% doi:10.1007/s00607-005-0139-x. S2CID 4561254.
%
% This should avoid cancellation
sum = 0.0;
cs  = 0.0;
ccs = 0.0;

for i = 1:length(input)
    t = sum + input(i);
    if abs(sum) >= abs(input(i))
        c = (sum - t) + input(i);
    else
        c = (input(i) - t) + sum;
    end
    sum = t;
    t = cs + c;
    if abs(cs) >= abs(c)
        cc = (cs - t) + c;
    else
        cc = (c - t) + cs;
    end
    cs = t;
    ccs = ccs + cc;
end

sum = sum + cs + ccs;

end
