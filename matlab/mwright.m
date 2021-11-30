function [w] = mwright(x,t,lambda,mu)
%MWRIGHT This routine computes the Wright function on the real line
%   Detailed explanation goes here

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

%% Build the contour
% These handle define the parabolic contour
gamma = 1;
z  = @(u) gamma*(1i*u + 1).^2;
zp = @(u) 2*1i*gamma*(1i*u + 1);

N = 30;
h = 1/N;
uk = -N:h:N;

w = zeros(size(x));

zuk = z(uk);
for k=1:length(w)
    w(k) = h/(2*pi*1i)*sum(exp(zuk.*t).*zuk.^(-mu).*exp(-abs(x(k)).*zuk.^(-lambda)).*zp(uk));
end

if norm(imag(w),'inf') < 10*eps
    w = real(w);
else
    warning("Non negligible immaginary part %1.2e",norm(imag(w),'inf'));
    w = real(w);
end

end

