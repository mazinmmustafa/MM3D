close all; clear; clc;

RCS_1 = load('data1.dat');
RCS_2 = load('data2.dat');

figure()
plot(RCS_1(:,1), RCS_1(:,2))

figure()
plot(RCS_2(:,1), RCS_2(:,3))