close all; clear; clc;

data1 = load('data1.dat');

figure()
hold on
plot(data1(:,1), data1(:,2))
plot(data1(:,1), data1(:,3))
plot(data1(:,1), data1(:,4))
hold off