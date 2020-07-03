clear all
close all
clc

% wdmp = cell2mat(struct2cell(load('wbestDMP1.mat','w')));
% [wdmp,Nt,y0dot,y0,x0,alphax,tau,dt,sigma_or_rho,c_or_alpha,alphay,betay,g,desiredtraj,ShowPlot,OneforDMPtwoforRMP]...
%     = cell2mat(struct2cell(load('wbestDMP1.mat')));
load('wbestDMP1.mat')
wdmp = w;
[d,dcomponents,ftrace1] = Eval(wdmp, Nt, y0dot, y0, x0, alphax, tau, dt, sigma, c, alphay, betay,g, desiredtraj, 1, OneforDMPtwoforRMP);

% wrmp = cell2mat(struct2cell(load('wbestRMP1.mat','w')));
load('wbestRMP1.mat')
wrmp = w;
[d,dcomponents,ftrace2] = Eval(wrmp, Nt, y0dot, y0, x0, alphax, tau, dt, sigma, c, alphay, betay,g, desiredtraj, 1, OneforDMPtwoforRMP);

figure(7)
hold on
plot((1:Nt)*dt,ftrace1(:,1),'ro','MarkerSize',0.7)
plot((1:Nt)*dt,ftrace1(:,2),'bo','MarkerSize',0.7)

plot((1:Nt)*dt,ftrace2(:,1),'Linewidth',1.5,'Color',[0.6 0 0])
plot((1:Nt)*dt,ftrace2(:,2),'Linewidth',1.5,'Color',[0 0 0.6])

figure(8)
hold on
plot(ftrace1(:,1),ftrace1(:,2),'.k')
plot(ftrace2(:,1),ftrace2(:,2),'-k')