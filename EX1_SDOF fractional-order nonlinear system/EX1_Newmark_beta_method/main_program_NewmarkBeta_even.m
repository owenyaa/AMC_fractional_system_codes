clear;clc;
% close all;
global h Tdata

%待识别参数
e=-0.8;
alpha=0.5;
w0=1;
k=0;
para_constant=1.1;

% mass_matrix=1; nonliear_stiffness_matrix=1;
% epsilon=-0.8; alpha=0.5;
% omega_2=1; parameter_constant=1.1;


tf=20;h=0.01;Tdata=(0:h:tf)';
iteration=length(Tdata);


parameter_a=[e,alpha,w0,k,para_constant];
x_cal=cal_x_even(parameter_a);
x_cal=x_cal';
x=x_cal(:,1);dx=x_cal(:,2);ddx=x_cal(:,3);


x_cal=[x,dx,ddx];
savefile='simple_fre_datanewa.mat';
save(savefile,'Tdata','x_cal');
figure
subplot(2,1,1);
plot(Tdata(1:20:end),x_cal(1:20:end-1,1),'k.','MarkerSize',12);
h1=legend('$$TMRM$$','$$Newmark-\beta$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% h=legend('位移');
% set(h,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% xlabel('时间/s','Fontsize',15);
% ylabel('响应','Fontsize',15);
% %相图;
figure
plot(x_cal(1590:20:end-1,1),x_cal(1590:20:end-1,2),'k.','MarkerSize',12);
% xlabel('displacement','Fontsize',15);
% ylabel('velocity','Fontsize',15);
h1=legend('$$Newmark-\beta$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);


