%% Plots of the f_{\lambda,mu}(1,x) for Section 2

clear; clc; close all;

M = 100;
x = linspace(-5,5,M);

%% \nu = 1/2
nu = 1/2;
w1 = zeros(M,1);
for i=1:M
   w1(i) = mwright(x(i),1,-nu,1-nu); 
end
figure(1)
subplot(1,4,1);
plot(x,w1,'k-','LineWidth',2);
xlabel('x');
ylabel('$f_{-\nu,1-\nu}(1,x)$','Interpreter','latex');
title('$\nu = \frac{1}{2}$','Interpreter','latex');
axis([-5 5 0 1])
axis square

%% \nu = 3/8
nu = 3/8;
w1 = zeros(M,1);
for i=1:M
   w1(i) = mwright(x(i),1,-nu,1-nu); 
end
figure(1)
subplot(1,4,2);
plot(x,w1,'k-','LineWidth',2);
xlabel('x');
ylabel('$f_{-\nu,1-\nu}(1,x)$','Interpreter','latex');
title('$\nu = \frac{3}{8}$','Interpreter','latex');
axis([-5 5 0 1])
axis square

%% \nu = 1/4
nu = 1/4;
w1 = zeros(M,1);
for i=1:M
   w1(i) = mwright(x(i),1,-nu,1-nu); 
end
figure(1)
subplot(1,4,3);
plot(x,w1,'k-','LineWidth',2);
xlabel('x');
ylabel('$f_{-\nu,1-\nu}(1,x)$','Interpreter','latex');
title('$\nu = \frac{1}{4}$','Interpreter','latex');
axis([-5 5 0 1])
axis square

%% \nu = 1/8
nu = 1/8;
w1 = zeros(M,1);
for i=1:M
   w1(i) = mwright(x(i),1,-nu,1-nu); 
end
figure(1)
subplot(1,4,4);
plot(x,w1,'k-','LineWidth',2);
xlabel('x');
ylabel('$f_{-\nu,1-\nu}(1,x)$','Interpreter','latex');
title('$\nu = \frac{1}{8}$','Interpreter','latex');
axis([-5 5 0 1])
axis square


