%% Plot the results for the Mainardi test function in double precision

clear; clc; close all;

result = system('../build/mainarditest');
if (result ~= 0)
    error("Error in running Fortran code");
end

m1 = csvread('mainardi1.out');

figure(1)
plot(m1(:,1),abs(m1(:,2)-m1(:,3))./abs(m1(:,3)),'r--');


