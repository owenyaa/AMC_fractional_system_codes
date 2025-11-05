clear;
clc; %close all;
global tf h Tdata

%待识别参数
e=-0.8;
alpha1=1.5;
alpha2=0.9;
k=1;
w0=1;
tf=15;h=0.01;Tdata=(0:h:tf)';
iteration=length(Tdata);


parameter_a=[e,alpha1,alpha2,k,w0];
x_cal=cal_x(parameter_a);
x_cal=x_cal';
x=x_cal(:,1);dx=x_cal(:,2);ddx=x_cal(:,3);


x_cal=[x,dx,ddx];
savefile='simple_fre_datanewa.mat';
save(savefile,'Tdata','x_cal');
figure
plot(Tdata,x_cal(1:end-1,1),'k-','LineWidth',1);

% h=legend('displacement');
% set(h,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','','FontSize',15,'LineWidth',1.5);
% xlabel('time/s','Fontsize',15);
% ylabel('displacement','Fontsize',15);
%相图;
figure
plot(x_cal(fix(0.8*length(Tdata)):end,1),x_cal(fix(0.8*length(Tdata)):end,2),'k-','LineWidth',1);
% xlabel('displacement','Fontsize',15);
% ylabel('velocity','Fontsize',15);
h1=legend('$$Newmark-\beta$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);


