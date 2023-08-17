close all; clear; clc;
%%
Data1   =   load('data1.dat');
Data2   =   load('data2.dat');
%%
figure()
hold on
plot(Data1(:,1),Data1(:,2),'-k','LineWidth',1)
plot(Data1(:,1),Data1(:,3),'--k','LineWidth',1)
plot(Data1(:,1),Data1(:,4),'-.k','LineWidth',1)
hold off
xlabel('$x$ [cm]','Interpret','Latex','FontSize',14)
ylabel('$|E|$ [V/m]','Interpret','Latex','FontSize',14)
set(gca,'TickLabel','Latex','FontSize',14)
legend('$E_{x}$','$E_{y}$','$E_{z}$',...
    'Interpreter','Latex','FontSize',14)
xlim([-100 +100])
%%
exportgraphics(gcf,'Figure1.pdf','ContentType','vector')
%%
figure()
hold on
plot(Data2(:,1),Data2(:,2),'-k','LineWidth',1)
plot(Data2(:,1),Data2(:,3),'--k','LineWidth',1)
plot(Data2(:,1),Data2(:,4),'-.k','LineWidth',1)
hold off
xlabel('$x$ [cm]','Interpret','Latex','FontSize',14)
ylabel('$|H|$ [mA/m]','Interpret','Latex','FontSize',14)
set(gca,'TickLabel','Latex','FontSize',14)
legend('$H_{x}$','$H_{y}$','$H_{z}$',...
    'Interpreter','Latex','FontSize',14)
xlim([-100 +100])
%%
exportgraphics(gcf,'Figure2.pdf','ContentType','vector')
%%