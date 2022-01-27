%% Comparison of Cauchy and Time Integrator
% This code makes a comparison between solving the time-fractiona
% differential equation with the expression of the solution obtained by
% means of the Wright function and by using the Fractional Trapezoidal
% rule. To run the example you need the FLMM2 code by R. Garrappa.
%
%   [1] C. Lubich, Discretized fractional calculus, SIAM J. Numer. Anal.
%   17(3) (1986), 704719 
%
%   [2] R. Garrappa, Trapezoidal methods for fractional differential
%   equations: Theoretical and computational aspects. Mathematics and
%   Computers in Simulation, DOI: 10.1016/j.matcom.2013.09.012 

clear; clc; close all;

addpath('flmm2/');

%% Spatial Domain
xlr=20; np=256; dx=2*xlr/np;
x=linspace(-xlr,xlr-dx,np);

%% Building Spatial Discretization
e = ones(np,1);
T = spdiags([-e,2*e,-e]./dx.^2,-1:1,np,np);

%% Time Domain
nu = 0.3;
alpha = 2*nu;
t0 = 0;
tfinal = 1;

%% Initial Data
g = @(x) real((abs(x) <= 1));
y0 = g(x).';
h = 0.01;

%% Solution with Fractional Linear Multistep Method
fdefun = @(t,y) -T*y;
Jfdefun = @(t,y) -T;
[t, y] = flmm2(alpha,fdefun,Jfdefun,t0,tfinal,y0,h);


%% Solution with Wright function
tic;
Gc = @(x,t) mwright(x,t,-nu,1-nu)/2;
Gcx = Gc(x,t(end));
fu1 = fft(y0.');
fGcx  = fft(Gcx);
fgt = fGcx.*fu1;
ysol = dx*fftshift( ifft( fgt ) );
timew = toc;

%% Compute timings
h = logspace(-3,-1,4);
npval = [2048,1024,512,256];
ind = 1;
for hval = h
    % Spatial Domain
    xlr=20;
    np = npval(ind); %floor((sqrt(4*xlr^2/hval^(2-alpha))));
    dx=2*xlr/np;
    xh=linspace(-xlr,xlr-dx,np);
    
    % Building Spatial Discretization
    e = ones(np,1);
    T = spdiags([-e,2*e,-e]./dx.^2,-1:1,np,np);
    y0 = g(xh).';
    
    fdefun = @(t,y) -T*y;
    Jfdefun = @(t,y) -T;
    fprintf('Launching FLMM np = %d h = %f.\n',np,hval);
    tic;
    [t, yh] = flmm2(alpha,fdefun,Jfdefun,t0,tfinal,y0,hval);
    time(ind) = toc;
    timeperstep(ind) = time(ind)/length(t);
    
    fprintf('Launching Wright function np = %d.\n',np);
    tic;
    Gc = @(x,t) mwright(xh,t,-nu,1-nu)/2;
    Gcx = Gc(xh,t(end));
    fu1 = fft(y0.');
    fGcx  = fft(Gcx);
    fgt = fGcx.*fu1;
    ysol = dx*fftshift( ifft( fgt ) );
    timew(ind) = toc;
    err(ind) = norm(yh(:,end)-ysol.','inf');
    
    
    figure(2)
    subplot(2,2,ind);
    semilogy(xh,abs(ysol-yh(:,end).'),'k-','LineWidth',2);
    title(sprintf('$n = %d$, $h = %1.3f$',np,hval),'Interpreter','latex');
    
    
    ind = ind + 1;
end

%% Comparison figures
figure(1)
subplot(2,3,[1,4])
plot(x,y(:,end),'r-',x,ysol,'b--','LineWidth',2)
legend({'FLMM','Wright'},'Location','northeast');
subplot(2,3,[2,5])
semilogy(x,abs(y(:,end)-ysol.'),'k','LineWidth',2)
title('Absolute Difference')
subplot(2,3,3)
loglog(time,err,'kx-','LineWidth',2)
xlabel('Time to solution (s)');
ylabel('$\|\cdot\|$-error','Interpreter','LaTeX');
subplot(2,3,6)
bar(1:4,[timew.',time.',timeperstep.'])
set(gca,'YScale','log')
legend({'Time Wright','Time FLMM','Time FLMM/N° Steps'})