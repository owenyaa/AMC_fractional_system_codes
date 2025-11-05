clear; clc; %close all;
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
load x_2_num_all1.mat
index_bifurcation=1;
for alpha=0.6:0.002:0.73
    alpha_inter(index_bifurcation)=num_x_cal(index_bifurcation).alpha;
    num=num_x_cal(index_bifurcation).x_cal;
    xmax(index_bifurcation).xmax=getmax(num(40000:end,1));
    index_bifurcation=index_bifurcation+1;
end


figure;
for i=1:1:length(0.6:0.002:0.73)
    QQ=alpha_inter(i)*ones(1,length(xmax(i).xmax));
    plot(QQ,xmax(i).xmax,'r.','MarkerSize',2);
    hold on;
end

index_bifurcation=67;
index_bifurcation_alpha=1;
for alpha=0.732:0.004:0.9
    alpha_inter(index_bifurcation_alpha)=num_x_cal(index_bifurcation).alpha;
    num=num_x_cal(index_bifurcation).x_cal;
    xmax(index_bifurcation_alpha).xmax=getmax(num(40000:end,1));
    index_bifurcation=index_bifurcation+2;
    index_bifurcation_alpha=index_bifurcation_alpha+1;
end


hold on;
for i=1:1:length(0.732:0.004:0.9)
    QQ=alpha_inter(i)*ones(1,length(xmax(i).xmax));
    plot(QQ,xmax(i).xmax,'b.','MarkerSize',6);
    hold on;
end