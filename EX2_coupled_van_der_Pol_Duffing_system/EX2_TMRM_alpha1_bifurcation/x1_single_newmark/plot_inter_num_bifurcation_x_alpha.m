% 绘制x2
clear; clc; %close all;
% uiopen('E:\Dropbox\科研\博后第三年\王海苏\分数阶系统\wanghaisu_new\程序\EX3_fractional_van_der_Pol_TMRM_new\alpha1_bifurcation\x1_single_newmark\untitled.fig',1);
load x1_num_xcal_1_1.976.mat;


index_bifurcation=1;
index_bifurcation1=1;
for alpha=1:0.001:1.394
    alpha_inter(index_bifurcation1)=num_xcal(index_bifurcation).alpha1;
    num=num_xcal(index_bifurcation).x_cal;
    xmax(index_bifurcation1).xmax=getmax(num(30000:end,1));
    index_bifurcation=index_bifurcation+1;
    index_bifurcation1=index_bifurcation1+1;
end


hold on;
for i=1:1:length(1:0.001:1.394)
    QQ=alpha_inter(i)*ones(1,length(xmax(i).xmax));
    plot(QQ,xmax(i).xmax,'r.','MarkerSize',2);
    hold on;
end

index_bifurcation=397;
index_bifurcation1=1;
for alpha=1.396:0.015:1.577
    alpha_inter(index_bifurcation1)=num_xcal(index_bifurcation).alpha1;
    num=num_xcal(index_bifurcation).x_cal;
    xmax(index_bifurcation1).xmax=getmax(num(40000:end,1));
    index_bifurcation=index_bifurcation+15;
    index_bifurcation1=index_bifurcation1+1;
end


hold on;
for i=1:1:length(1.396:0.015:1.577)
    QQ=alpha_inter(i)*ones(1,length(xmax(i).xmax));
    plot(QQ,xmax(i).xmax,'r.','MarkerSize',10);
    hold on;
end

index_bifurcation=579;
index_bifurcation1=1;
for alpha=1.578:0.001:1.674
    alpha_inter(index_bifurcation1)=num_xcal(index_bifurcation).alpha1;
    num=num_xcal(index_bifurcation).x_cal;
    xmax(index_bifurcation1).xmax=getmax(num(30000:end,1));
    index_bifurcation=index_bifurcation+1;
    index_bifurcation1=index_bifurcation1+1;
end


hold on;
for i=1:1:length(1.578:0.001:1.674)
    QQ=alpha_inter(i)*ones(1,length(xmax(i).xmax));
    plot(QQ,xmax(i).xmax,'r.','MarkerSize',2);
    hold on;
end


index_bifurcation=676;
index_bifurcation1=1;
for alpha=1.675:0.015:1.705
    alpha_inter(index_bifurcation1)=num_xcal(index_bifurcation).alpha1;
    num=num_xcal(index_bifurcation).x_cal;
    xmax(index_bifurcation1).xmax=getmax(num(40000:end,1));
    index_bifurcation=index_bifurcation+15;
    index_bifurcation1=index_bifurcation1+1;
end


hold on;
for i=1:1:length(1.675:0.015:1.705)
    QQ=alpha_inter(i)*ones(1,length(xmax(i).xmax));
    plot(QQ,xmax(i).xmax,'r.','MarkerSize',10);
    hold on;
end

index_bifurcation=706;
index_bifurcation1=1;
for alpha=1.705:0.005:1.72
    alpha_inter(index_bifurcation1)=num_xcal(index_bifurcation).alpha1;
    num=num_xcal(index_bifurcation).x_cal;
    xmax(index_bifurcation1).xmax=getmax(num(40000:end,1));
    index_bifurcation=index_bifurcation+5;
    index_bifurcation1=index_bifurcation1+1;
end


hold on;
for i=1:1:length(1.705:0.005:1.72)
    QQ=alpha_inter(i)*ones(1,length(xmax(i).xmax));
    plot(QQ,xmax(i).xmax,'r.','MarkerSize',10);
    hold on;
end

index_bifurcation=721;
index_bifurcation1=1;
for alpha=1.72:0.015:1.874
    alpha_inter(index_bifurcation1)=num_xcal(index_bifurcation).alpha1;
    num=num_xcal(index_bifurcation).x_cal;
    xmax(index_bifurcation1).xmax=getmax(num(40000:end,1));
    index_bifurcation=index_bifurcation+15;
    index_bifurcation1=index_bifurcation1+1;
end


hold on;
for i=1:1:length(1.72:0.015:1.874)
    QQ=alpha_inter(i)*ones(1,length(xmax(i).xmax));
    plot(QQ,xmax(i).xmax,'r.','MarkerSize',10);
    hold on;
end

clear; clc; %close all;
load x1_num_xcal_1_1.031_green.mat;
index_bifurcation=1;
index_bifurcation1=1;
for alpha=1:0.015:1.031
    alpha_inter(index_bifurcation1)=num_xcal(index_bifurcation).alpha1;
    num=num_xcal(index_bifurcation).x_cal;
    xmax(index_bifurcation1).xmax=getmax(num(40000:end,1));
    index_bifurcation=index_bifurcation+15;
    index_bifurcation1=index_bifurcation1+1;
end


hold on;
for i=1:1:length(1:0.015:1.031)
    QQ=alpha_inter(i)*ones(1,length(xmax(i).xmax));
    plot(QQ,xmax(i).xmax,'r.','MarkerSize',10);
    hold on;
end


