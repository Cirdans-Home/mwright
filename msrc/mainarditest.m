%% Plot the results for the Mainardi test function

clear; clc; close all;

m1 = csvread('mainardi1.out');

figure(1)
plot(m1(:,1),abs(m1(:,2)-m1(:,3))./abs(m1(:,3)),'r--');

m2 = csvread('mainardi1.out');

figure(2)
plot(m2(:,1),abs(m2(:,2)-m2(:,3))./abs(m2(:,3)),'r--');
