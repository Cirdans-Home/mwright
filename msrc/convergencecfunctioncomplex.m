%% Convergence with respect to the values computed with the ARB series
% This test uses the implementation with the ARB library of the series to 
% compute benchmark values for the comparison with the approach given by
% the inversion of the Laplace transform - Complex values of :math:`\mu`

clear; clc; close all;

%% Parameters of the function
lambda = -0.2; 
mu = -0.5 - 2*1i;
n = 1000;
prec = 128;

%% Set of points for the evaluation
m = 100;
x = linspace(-5,0,m).';
fid = fopen('input.inp','w+');
fprintf(fid,"%d",m);
for i=1:m
   fprintf(fid,"%1.16f\n",x(i));
end

%% Execute the ARB code
stringa = sprintf("../build/zwrighttest %1.16f %1.16f %1.16f %d %d input.inp",...
    lambda, real(mu), imag(mu), n, prec);
system(stringa);

system('sed -i "s/- /-/g" wrighttest.out');
system('sed -i "s/+ /+/g" wrighttest.out');

%% Read the result to MATLAB for plotting
wseries = zeros(m,1);
fid = fopen('wrighttest.out','r');
for i=1:m
    line = fgetl(fid);
    res = sscanf(line,"(%e %e)");
    wseries(i) = res(1) + res(2)*1i;
end
fclose(fid);

%% Compute the values with mwright function
t = 1;

[wmatlab,N] = mwright(x,t,lambda,mu);
N = unique(round(linspace(2,N,4)));

figure(1)
subplot(1,2,1);
plot(x,real(wmatlab),'r-',x,real(wseries),'b--','LineWidth',2)
title("Real");
subplot(1,2,2);
plot(x,imag(wmatlab),'r-',x,imag(wseries),'b--','LineWidth',2)
title("Imag");
legend('Matlab','Series')


for Nval = N
    wmatlab = mwright(x,t,lambda,mu,Nval);
    figure(2)
    hold on
    semilogy(x,abs(wmatlab-wseries)./abs(wseries),... %
        'DisplayName',sprintf('N = %d',Nval))
    hold off
end

title(sprintf('$\\lambda =  %1.2f$, $\\mu = %1.2f $',lambda,mu),...
    'Interpreter','latex')
legend
set(gca,'YScale','log')