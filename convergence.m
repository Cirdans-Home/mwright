%% Convergence 

M = 100;
x = linspace(-5,5,M).';

%% \nu = 1/2
true = @(z) exp(-z.^2/4)./sqrt(pi);
nu = 1/2;
w = zeros(M,1);
for i=1:M
   w(i) = mwright(x(i),1,-nu,1-nu); 
end

figure(2)
hold on
semilogy(x,abs(w-true(x))./true(x),'DisplayName','N = 20','LineWidth',2)
hold off

%% \nu = 1/3
true = @(z) (3^(2/3))*airy(0,abs(z)./(3^(1/3)),0);
nu = 1/3;
w = zeros(M,1);
for i=1:M
   w(i) = mwright(x(i),1,-nu,1-nu); 
end

figure(3)
hold on
semilogy(x,abs(w-true(x))./true(x),'DisplayName','N = 20','LineWidth',2)
hold off

%% \lambda = -1/2 \mu = 0
true = @(z) abs(z).*exp(-abs(z).^2/4)./(2*sqrt(pi)); 
nu = 1/3;
w = zeros(M,1);
for i=1:M
   w(i) = mwright(x(i),1,-1/2,0); 
end

figure(4)
hold on
semilogy(x,abs(w-true(x))./true(x),'DisplayName','N = 20','LineWidth',2)
hold off


