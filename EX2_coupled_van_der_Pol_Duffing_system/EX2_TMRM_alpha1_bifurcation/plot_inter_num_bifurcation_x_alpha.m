% clear; clc; %close all;
% load x_1_num_all.mat
% index_bifurcation=1;
% for alpha=0.6:0.002:0.684
%     alpha_inter(index_bifurcation)=num_x_cal(index_bifurcation).alpha;
%     num=num_x_cal(index_bifurcation).x_cal;
%     xmax(index_bifurcation).xmax=getmax(num(40000:end,1));
%     index_bifurcation=index_bifurcation+1;
% end
% 
% 
% figure;
% for i=1:1:length(0.6:0.002:0.684)
%     QQ=alpha_inter(i)*ones(1,length(xmax(i).xmax));
%     plot(QQ,xmax(i).xmax,'r.','MarkerSize',2);
%     hold on;
% end
% 
% index_bifurcation=44;
% index_bifurcation_alpha=1;
% for alpha=0.686:0.004:0.9
%     alpha_inter(index_bifurcation_alpha)=num_x_cal(index_bifurcation).alpha;
%     num=num_x_cal(index_bifurcation).x_cal;
%     xmax(index_bifurcation_alpha).xmax=getmax(num(40000:end,1));
%     index_bifurcation=index_bifurcation+2;
%     index_bifurcation_alpha=index_bifurcation_alpha+1;
% end
% 
% 
% hold on;
% for i=1:1:length(0.686:0.004:0.9)
%     QQ=alpha_inter(i)*ones(1,length(xmax(i).xmax));
%     plot(QQ,xmax(i).xmax,'b.','MarkerSize',6);
%     hold on;
% end

% 绘制x2
clear; clc; %close all;
% uiopen('E:\Dropbox\科研\博后第三年\王海苏\分数阶系统\wanghaisu_new\程序\EX3_fractional_van_der_Pol_TMRM_new\alpha1_bifurcation\bifurcation_x2_alpha1_1_1.95.fig',1);
load x2_num_xcal_1_1.906.mat
index_bifurcation=1;
index_bifurcation1=1;
for alpha=1:0.015:1.532
    alpha_inter(index_bifurcation1)=num_xcal(index_bifurcation).alpha1;
    num=num_xcal(index_bifurcation).x_cal;
    xmax(index_bifurcation1).xmax=getmax(num(48000:end,1));
    index_bifurcation=index_bifurcation+15;
    index_bifurcation1=index_bifurcation1+1;
end


hold on;
for i=1:1:length(1:0.015:1.532)
    QQ=alpha_inter(i)*ones(1,length(xmax(i).xmax));
    plot(QQ,xmax(i).xmax,'r.','MarkerSize',10);
    hold on;
end


index_bifurcation=680;
index_bifurcation1=1;
for alpha=1.679:0.015:1.726
    alpha_inter(index_bifurcation1)=num_xcal(index_bifurcation).alpha1;
    num=num_xcal(index_bifurcation).x_cal;
    xmax(index_bifurcation1).xmax=getmax(num(48000:end,1));
    index_bifurcation=index_bifurcation+15;
    index_bifurcation1=index_bifurcation1+1;
end


hold on;
for i=1:1:length(1.679:0.015:1.726)
    QQ=alpha_inter(i)*ones(1,length(xmax(i).xmax));
    plot(QQ,xmax(i).xmax,'r.','MarkerSize',10);
    hold on;
end


index_bifurcation=727;
index_bifurcation1=1;
for alpha=1.726:0.015:1.906
    alpha_inter(index_bifurcation1)=num_xcal(index_bifurcation).alpha1;
    num=num_xcal(index_bifurcation).x_cal;
    xmax(index_bifurcation1).xmax=getmax(num(48000:end,1));
    index_bifurcation=index_bifurcation+15;
    index_bifurcation1=index_bifurcation1+1;
end


hold on;
for i=1:1:length(1.726:0.015:1.906)
    QQ=alpha_inter(i)*ones(1,length(xmax(i).xmax));
    plot(QQ,xmax(i).xmax,'r.','MarkerSize',10);
    hold on;
end


index_bifurcation=534;
index_bifurcation1=1;
for alpha=1.533:0.005:1.678
    alpha_inter(index_bifurcation1)=num_xcal(index_bifurcation).alpha1;
    num=num_xcal(index_bifurcation).x_cal;
    xmax(index_bifurcation1).xmax=getmax(num(30000:end,1));
    index_bifurcation=index_bifurcation+5;
    index_bifurcation1=index_bifurcation1+1;
end


hold on;
for i=1:1:length(1.533:0.005:1.678)
    QQ=alpha_inter(i)*ones(1,length(xmax(i).xmax));
    plot(QQ,xmax(i).xmax,'r.','MarkerSize',2);
    hold on;
end






% clear; clc; %close all;
% load x2_num_xcal_1.401_1.65.mat
% index_bifurcation=1;
% for alpha=1.401:0.001:1.65
%     alpha_inter(index_bifurcation)=num_xcal(index_bifurcation).alpha1;
%     num=num_xcal(index_bifurcation).x_cal;
%     xmax(index_bifurcation).xmax=getmax(num(48000:end,1));
%     index_bifurcation=index_bifurcation+1;
% end
% 
% 
% hold on;
% for i=1:1:length(1.401:0.001:1.65)
%     QQ=alpha_inter(i)*ones(1,length(xmax(i).xmax));
%     plot(QQ,xmax(i).xmax,'b.','MarkerSize',2);
%     hold on;
% end
% 
% 
% 
% 
% 
% 
% 
% index_bifurcation=67;
% index_bifurcation_alpha=1;
% for alpha=0.732:0.004:0.9
%     alpha_inter(index_bifurcation_alpha)=num_x_cal(index_bifurcation).alpha;
%     num=num_x_cal(index_bifurcation).x_cal;
%     xmax(index_bifurcation_alpha).xmax=getmax(num(40000:end,1));
%     index_bifurcation=index_bifurcation+2;
%     index_bifurcation_alpha=index_bifurcation_alpha+1;
% end
% 
% 
% hold on;
% for i=1:1:length(0.732:0.004:0.9)
%     QQ=alpha_inter(i)*ones(1,length(xmax(i).xmax));
%     plot(QQ,xmax(i).xmax,'b.','MarkerSize',6);
%     hold on;
% end