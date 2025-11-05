% 绘制x2
% clear; clc; 
close all;
uiopen('E:\Dropbox\科研\博后第三年\王海苏\分数阶系统\wanghaisu_new\程序\EX3_fractional_van_der_Pol_TMRM_new\alpha1_bifurcation\untitled.fig',1)
% load x2_num_xcal_1_1.906.mat;
index_bifurcation=471;
index_bifurcation1=1;
for alpha=1.470:0.01:1.532
    alpha_inter(index_bifurcation1)=num_xcal(index_bifurcation).alpha1;
    num=num_xcal(index_bifurcation).x_cal;
    xmax(index_bifurcation1).xmax=getmax(num(48000:end,1));
    index_bifurcation=index_bifurcation+10;
    index_bifurcation1=index_bifurcation1+1;
end


hold on;
for i=1:1:length(1.47:0.01:1.532)
    QQ=alpha_inter(i)*ones(1,length(xmax(i).xmax));
    plot(QQ,xmax(i).xmax,'r.','MarkerSize',10);
    hold on;
end


index_bifurcation=680;
index_bifurcation1=1;
for alpha=1.679:0.01:1.739
    alpha_inter(index_bifurcation1)=num_xcal(index_bifurcation).alpha1;
    num=num_xcal(index_bifurcation).x_cal;
    xmax(index_bifurcation1).xmax=getmax(num(48000:end,1));
    index_bifurcation=index_bifurcation+10;
    index_bifurcation1=index_bifurcation1+1;
end


hold on;
for i=1:1:length(1.679:0.01:1.739)
    QQ=alpha_inter(i)*ones(1,length(xmax(i).xmax));
    plot(QQ,xmax(i).xmax,'r.','MarkerSize',10);
    hold on;
end


% index_bifurcation=727;
% index_bifurcation1=1;
% for alpha=1.726:0.015:1.906
%     alpha_inter(index_bifurcation1)=num_xcal(index_bifurcation).alpha1;
%     num=num_xcal(index_bifurcation).x_cal;
%     xmax(index_bifurcation1).xmax=getmax(num(48000:end,1));
%     index_bifurcation=index_bifurcation+15;
%     index_bifurcation1=index_bifurcation1+1;
% end
% 
% 
% hold on;
% for i=1:1:length(1.726:0.015:1.906)
%     QQ=alpha_inter(i)*ones(1,length(xmax(i).xmax));
%     plot(QQ,xmax(i).xmax,'r.','MarkerSize',10);
%     hold on;
% end


index_bifurcation=534;
index_bifurcation1=1;
for alpha=1.533:0.001:1.678
    alpha_inter(index_bifurcation1)=num_xcal(index_bifurcation).alpha1;
    num=num_xcal(index_bifurcation).x_cal;
    xmax(index_bifurcation1).xmax=getmax(num(30000:end,1));
    index_bifurcation=index_bifurcation+1;
    index_bifurcation1=index_bifurcation1+1;
end


hold on;
for i=1:1:length(1.533:0.001:1.678)
    QQ=alpha_inter(i)*ones(1,length(xmax(i).xmax));
    plot(QQ,xmax(i).xmax,'r.','MarkerSize',2);
    hold on;
end
