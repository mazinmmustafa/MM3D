close all; clear; clc;
%%
Data1   =   load('data1.dat');
Data2   =   load('data2.dat');
%%
figure()
plot(Data1(:,1),Data1(:,2),'-k','LineWidth',1)
xlabel('$\theta$ [deg]','Interpret','Latex','FontSize',14)
ylabel('$\sigma_{\theta\theta}$ [dB]','Interpret','Latex','FontSize',14)
set(gca,'TickLabel','Latex','FontSize',14)
xlim([0 180])
ylim([-10 +15])
%%
exportgraphics(gcf,'Figure72.pdf','ContentType','vector')
%%
figure()
plot(Data2(:,1),Data2(:,3),'-k','LineWidth',1)
xlabel('$\theta$ [deg]','Interpret','Latex','FontSize',14)
ylabel('$\sigma_{\varphi\varphi}$ [dB]','Interpret','Latex','FontSize',14)
set(gca,'TickLabel','Latex','FontSize',14)
xlim([0 180])
ylim([-10 +15])
%%
exportgraphics(gcf,'Figure73.pdf','ContentType','vector')
%%