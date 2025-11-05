clear; clc; close all;

% 定义全局系统参数
global alpha1 alpha2 mu k alpha beta epsilon A mass_matrix
alpha1=1; alpha2=1; mu=0.8; k=3; alpha=0.8; beta=0; epsilon=1.2;
mass_matrix=[1,0;0,1];damping_matrix=[-mu 0;-k alpha+k];stiffness_matrix=[1+k -k;0 -1];
A=[zeros(2),eye(2);-mass_matrix\stiffness_matrix,-mass_matrix\damping_matrix];


% 初始条件
y0 = [1.0; 0.0; 0.5; 0.3];
% y0 = [displacement(1,1); displacement(2,1); velocity(1,1); velocity(2,1)];
dt=0.01;
% 时间点
t = 0:dt:500;

% 设置数值积分的选项
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
index_bifurcation=1;
for alpha=0.8:-0.001:0.6
    damping_matrix=[-mu 0;-k alpha+k];
    A=[zeros(2),eye(2);-mass_matrix\stiffness_matrix,-mass_matrix\damping_matrix];
    % 求解微分方程
    [t, num] = ode45('equations', t, y0, options);
    % 提取解
    x1 = num(:, 1);
    x2 = num(:, 2);
    numer(index_bifurcation).alpha=alpha;
    numer(index_bifurcation).num=num;
    index_bifurcation=index_bifurcation+1;
end
figure;
plot(num(10000:end,1), num(10000:end,3), 'r-', 'LineWidth', 1);
% hold on;
% plot(t, sol(:,2), 'b-', 'LineWidth', 1);
% legend_handle = legend('$$x_1$$', '$$x_2$$');
set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);


figure;
plot(t(10000:end,1), num(10000:end,2), 'r-', 'LineWidth', 1);
% hold on;
% plot(t, sol(:,2), 'b-', 'LineWidth', 1);
% legend_handle = legend('$$x_1$$', '$$x_2$$');
% set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

% 频谱分析
temp1 = fix(length(x1) * 1 / 4);
temp2 = length(x1);
fs=1/dt;%采样频率
% 采样频率与时间间隔之间的关系： fs=1/dt
% 采样定理告诉我们，采样频率要大于信号频率的两倍。 
N=2^20;  %采样点数2^17
% N个采样点，经过FFT之后，就可以得到N个点的FFT结果。为了方便进行FFT运算，通常N取2的整数次方。
% 要精确到xHz，则需要采样长度为1/x秒的信号，并做FFT。
% 要提高频率分辨率，就需要增加采样点数
n=0:N-1;
t=n/fs;  % dt=1/fs 表示时间间隔   fs=1/dt
y=fft(num(temp1:temp2,1),N);  % 进行fft变换
% 假设采样频率为Fs，信号频率F，采样点数为N。那么FFT之后结果就是一个为N点的复数。
% 每一个点就对应着一个频率点。这个点的模值，就是该频率值下的幅度特性。
% y % 输出y看看fft之后的结果。
m=abs(y(1:N/2))*2/N; % 求信号的真实幅值
m=log(m);
f=2*pi*n*fs/N;  % 注意 乘2*pi 和不乘2*pi的区别%m=log10(m);
figure;
% plot(f(1:N/2)./(2*pi),m(1:N/2),'r-','LineWidth',1);
plot(f(1:N/2),m(1:N/2),'r-','LineWidth',1);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
h1=legend('$$Iteration steps$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
xlim([0 5]);

% for i=1:1:length(x1)
%     acc1(1:2,i)=-mass_matrix\(damping_matrix*sol(i,3:4)'+stiffness_matrix*sol(i, 1:2)'+...
%         [eta*beta*cos(omega_cos*t(i))*sol(i, 1)'+omega_R1^2*sol(i, 3)';0.75*epsilon_omega*Omegas*sol(i, 2)'.^2.*sol(i, 4)']);
% end
% acc1=acc1';