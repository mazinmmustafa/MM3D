close all; clear; clc;
%%
j       =   sqrt(-1);
Zmn     =   load('Z_mn_r.dat')+j*load('Z_mn_i.dat');
%%
figure()
pcolor(flipud(abs(Zmn)))
shading flat
colormap gray
axis equal
colorbar
%%
Ymn     =   inv(Zmn);
%%
figure()
pcolor(flipud(abs(Ymn)))
shading flat
colormap gray
axis equal
colorbar
%%