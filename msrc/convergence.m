%% Convergence for the Matlab Implementation
% This test verifies the convergence for the MATLAB implementation with
% respect to closed form solutions for the Mainardi version of the
% function.

clear; clc; close all;

M = 25;
x = linspace(-5,5,M).';

%% \nu = 1/2
true = @(z) exp(-z.^2/4)./sqrt(pi);
nu = 1/2;
w = zeros(M,1);
for i=1:M
   [w(i),N] = mwright(x(i),1,-nu,1-nu); 
end

Nvalues = unique(round(linspace(2,N,4)));

for N = Nvalues
    w = zeros(M,1);
    for i=1:M
        w(i) = mwright(x(i),1,-nu,1-nu,N);
    end
    
    figure(1)
    hold on
    semilogy(x,abs(w-true(x))./true(x),'DisplayName',sprintf('N = %d',N),'LineWidth',2)
    hold off
    
end

legend;
title('$\lambda = -1/2, \; \mu = 1/2$','Interpreter','latex')
set(gca,'YScale','log');

%% \nu = 1/3
true = @(z) (3^(2/3))*airy(0,abs(z)./(3^(1/3)),0);
nu = 1/3;
w = zeros(M,1);
for i=1:M
   [w(i),N] = mwright(x(i),1,-nu,1-nu); 
end

Nvalues = unique(round(linspace(2,N,4)));

for N = Nvalues
    w = zeros(M,1);
    for i=1:M
        w(i) = mwright(x(i),1,-nu,1-nu,N);
    end
    
    figure(2)
    hold on
    semilogy(x,abs(w-true(x))./true(x),'DisplayName',sprintf('N = %d',N),'LineWidth',2)
    hold off
    
end

legend;
title('$\lambda = -1/3, \; \mu = 2/3$','Interpreter','latex')
set(gca,'YScale','log');

%% \nu = 0
true = @(z) exp(-abs(z));
nu = 0;
w = zeros(M,1);
for i=1:M
   [w(i),N] = mwright(x(i),1,-nu,1-nu); 
end

Nvalues = unique(round(linspace(2,N,4)));

for N = Nvalues
    w = zeros(M,1);
    for i=1:M
        w(i) = mwright(x(i),1,-nu,1-nu,N);
    end
    
    figure(3)
    hold on
    semilogy(x,abs(w-true(x))./true(x),'DisplayName',sprintf('N = %d',N),'LineWidth',2)
    hold off
    
end

legend;
title('$\lambda = 0, \; \mu = 0$','Interpreter','latex')
set(gca,'YScale','log');