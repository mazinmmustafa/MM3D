close all; clear; clc;
%%
MeshData    =   load('mesh_data.dat');
BasisData   =   load('basis.dat');
%%
[N_tri, ~] 	=   size(MeshData);
[N_basis, ~]=   size(BasisData);
n           =   0.1;
%%
figure()
hold on
for i=1:N_tri
    x1 = MeshData(i,1); x2 = MeshData(i,4); x3 = MeshData(i,7); 
    y1 = MeshData(i,2); y2 = MeshData(i,5); y3 = MeshData(i,8);
    z1 = MeshData(i,3); z2 = MeshData(i,6); z3 = MeshData(i,9);
    fill3([x1 x2 x3],[y1 y2 y3],[z1 z2 z3],[1 1 1]) 
    %% Show numbers
    str = sprintf('%d',i);
    text((x1+x2+x3)/3,(y1+y2+y3)/3,(z1+z2+z3)/3,...
        str,'Interpret','Latex','FontSize',15) 
end
axis equal
% view([45 45])
view([0 90])
xlabel('$x$','Interpret','Latex','FontSize',15)
ylabel('$y$','Interpret','Latex','FontSize',15)
zlabel('$z$','Interpret','Latex','FontSize',15)
set(gca,'TickLabel','Latex','FontSize',15)
for i=1:N_basis
    x1 = BasisData(i,1);
    y1 = BasisData(i,2);
    z1 = BasisData(i,3);
    x2 = (BasisData(i,4)+BasisData(i,10))/2;
    y2 = (BasisData(i,5)+BasisData(i,11))/2;
    z2 = (BasisData(i,6)+BasisData(i,12))/2;
    plot3([x1 x2],[y1 y2],[z1 z2],'-b','LineWidth',1)
    x1 = BasisData(i,7);
    y1 = BasisData(i,8);
    z1 = BasisData(i,9);
    x2 = (BasisData(i,4)+BasisData(i,10))/2;
    y2 = (BasisData(i,5)+BasisData(i,11))/2;
    z2 = (BasisData(i,6)+BasisData(i,12))/2;
    plot3([x1 x2],[y1 y2],[z1 z2],'-r','LineWidth',1)
%     Normal unit vectors
%     xc = (BasisData(i,1)+BasisData(i,4)+BasisData(i,10))/3;
%     yc = (BasisData(i,2)+BasisData(i,5)+BasisData(i,11))/3;
%     zc = (BasisData(i,3)+BasisData(i,6)+BasisData(i,12))/3;
%     xn = BasisData(i,13);
%     yn = BasisData(i,14);
%     zn = BasisData(i,15);
%     plot3([xc xc+n*xn],[yc yc+n*yn],[zc zc+n*zn],'-k','LineWidth',1)
%     xc = (BasisData(i,4)+BasisData(i,7)+BasisData(i,10))/3;
%     yc = (BasisData(i,5)+BasisData(i,8)+BasisData(i,11))/3;
%     zc = (BasisData(i,6)+BasisData(i,9)+BasisData(i,12))/3;
%     xn = BasisData(i,16);
%     yn = BasisData(i,17);
%     zn = BasisData(i,18);
%     plot3([xc xc+n*xn],[yc yc+n*yn],[zc zc+n*zn],'-k','LineWidth',1)
end
hold off
%%



