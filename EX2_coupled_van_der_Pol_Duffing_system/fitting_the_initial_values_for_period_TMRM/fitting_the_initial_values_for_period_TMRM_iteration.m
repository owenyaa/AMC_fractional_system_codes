%% Author : GUANG_LIU  * owenyaa@gmail.com *
% Created Time : 2020-09-29 13:20
% Last Revised : GUANG_LIU ,2020-09-29
% Remark : 本程序的结构、功能和变量做一下说明%
%          本程序的主要功能是从稳态的数值响应(一般RK法获得)中，拟合出谐波系数。%
%          需要注意的是，拟合的周期(频率)会极大地影响拟合的精度。%
%          程序的主要部分是子函数time_to_fourier。%
%          对于仅含奇次谐波的响应，只要输入的频率精度够高，拟合出的偶次谐波系数会直接为0.%
clear;clc;close all;
% 数值计算结果
global num_degrees_freedom num_harmonics time_data temp1 temp2 num_harmonics_fitting
num_degrees_freedom=2;num_harmonics=25;
% load alpha_0.3_beta_0.3_omega_cos_2.1_Omegas_0.97_600_0.01.mat;
load matlab.mat;
dt=0.01;
tt=0:dt:500;
% num=x_cal(1:end,1:num_degrees_freedom);
temp1=fix(length(tt)*1/4);%temp1=20414;
temp2=length(tt);
% num=x_cal;

% fs=1/dt;%采样频率
% % 采样频率与时间间隔之间的关系： fs=1/dt
% % 采样定理告诉我们，采样频率要大于信号频率的两倍。 
% N=2^15;  %采样点数2^17
% % N个采样点，经过FFT之后，就可以得到N个点的FFT结果。为了方便进行FFT运算，通常N取2的整数次方。
% % 要精确到xHz，则需要采样长度为1/x秒的信号，并做FFT。
% % 要提高频率分辨率，就需要增加采样点数
% n=0:N-1;
% t=n/fs;  % dt=1/fs 表示时间间隔   fs=1/dt
% y=fft(num(temp1:temp2,1),N);  % 进行fft变换
% % 假设采样频率为Fs，信号频率F，采样点数为N。那么FFT之后结果就是一个为N点的复数。
% % 每一个点就对应着一个频率点。这个点的模值，就是该频率值下的幅度特性。
% % y % 输出y看看fft之后的结果。
% m=abs(y(1:N/2))*2/N; % 求信号的真实幅值
% m=log(m);
% fre=2*pi*n*fs/N;  % 注意 乘2*pi 和不乘2*pi的区别%m=log10(m);
% figure;
% % plot(f(1:N/2)./(2*pi),m(1:N/2),'r-','LineWidth',1);
% plot(fre(1:N/2),m(1:N/2),'r-','LineWidth',1);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% h1=legend('$$Iteration steps$$');
% set(h1,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% xlim([0 5]);

% num_harmonics_fitting=2*num_harmonics;
num_harmonics_fitting=num_harmonics;
wn=1.38058;
% wn=0.738288;

w_input=wn;
time_data=tt(temp1:temp2);
vector=[];
for nx=1:num_degrees_freedom
    temp=converse(time_to_fourier(w_input,time_data,num(temp1:temp2,nx)));
    vector=[vector temp];%....x_input, for Jacobi and residual
end
%vector 第一项为常数项
%开始将vector转换为parameter_a
parameter_a=zeros(num_harmonics+2,2*num_degrees_freedom);
parameter_a(1,1)=wn;
parameter_a(2,1)=vector(1,1);parameter_a(2,3)=vector(1,2);
for i=1:num_degrees_freedom
    jj=3;
    for ii=1:num_harmonics_fitting
        %         if mod(ii,2)==1
        parameter_a(jj,2*i-1)=vector(2*ii,i);
        parameter_a(jj,2*i)=vector(2*ii+1,i);
        jj=jj+1;
        %         end
    end
end

% initial_frequency = parameter_a(1, 1);
% harmonic_parameters = parameter_a(3:end, :);
% 
% % Initialize arrays for displacement, velocity, and acceleration
% displacement = zeros(num_degrees_freedom, length(time_data));
% velocity = zeros(num_degrees_freedom, length(time_data));
% acceleration = zeros(num_degrees_freedom, length(time_data));
% 
% % Compute displacement, velocity, and acceleration for each degree of freedom
% for dof_index = 1:num_degrees_freedom
%     for harm_index = 1:num_harmonics
%         harmonic_term = harm_index  * initial_frequency * time_data;
%         displacement(dof_index, :) = displacement(dof_index, :) + harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(harmonic_term) + harmonic_parameters(harm_index, 2 * dof_index) * sin(harmonic_term);
%         velocity(dof_index, :) = velocity(dof_index, :) - initial_frequency * harm_index  * harmonic_parameters(harm_index, 2 * dof_index - 1) * sin(harmonic_term) + initial_frequency * harm_index  * harmonic_parameters(harm_index, 2 * dof_index) * cos(harmonic_term);
%         acceleration(dof_index, :) = acceleration(dof_index, :) - (initial_frequency * harm_index )^2 * harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(harmonic_term) - (initial_frequency * harm_index )^2 * harmonic_parameters(harm_index, 2 * dof_index) * sin(harmonic_term);
%     end
%     displacement(dof_index, :) = displacement(dof_index, :) + parameter_a(2, 2*dof_index-1);
% end
% figure;
% plot(time_data,displacement(1,:),'r-','LineWidth',1.5);
% hold on;
% plot(time_data,displacement(2,:),'b-','LineWidth',1.5);
% 
% figure;
% plot(displacement(1,:),velocity(1,:),'r-','LineWidth',1.5);



%% NS迭代出频率和谐波系数
% 当指定频率之后，通过位移可以直接拟合出谐波系数，同样的道理，通过速度和加速度都可以拟合出系数
% 直接拟合的精度，和指定的频率的准确性密切相关
% 为此，可以通过NS迭代，直接精准计算出频率和谐波系数
% 本程序采用的是位移和速度进行NS迭代  其实也可以采用加速度，加速度需要对RK法的结果进一步计算获得
for i=1:1:num_degrees_freedom
    temp_fre=[0,0];
    temp_fre(1,1)=parameter_a(1,1);
    temp_parameter_a=parameter_a(2:end,2*i-1:2*i);
    dof(i).parameter_a=[temp_fre;temp_parameter_a];
end

for i=1:1:num_degrees_freedom
    temp_sensitivity_parameter_da=0.1;
    while norm(temp_sensitivity_parameter_da)>=1e-10
        %         aaa=dof(i).parameter_a;aaaa=;

        % temp_fre=[0,0];
        % temp_fre(1,1)=parameter_a(1,1);
        % temp_parameter_a=parameter_a(2:end,2*i-1:2*i);
        % temp_parameter_a=[temp_fre;temp_parameter_a];
        % temp_num=f(num(:,i),temp_parameter_a);
        % temp_sensitivity_parameter_da=-df(temp_parameter_a)\temp_num;

        temp_num=f(num(:,i),dof(i).parameter_a);
        temp_sensitivity_parameter_da=-df(dof(i).parameter_a)\temp_num;        
        temp_real_frequency = zeros(2, 2 * 1);
        temp_real_frequency(1, 1) = temp_sensitivity_parameter_da(1, 1);
        temp_real_frequency(2, 1) = temp_sensitivity_parameter_da(end, 1);
        coeff_updates = reshape(temp_sensitivity_parameter_da(2:end-1,:), 2, 1*num_harmonics);
        coeff_updates = coeff_updates';
        coeff_updates=[temp_real_frequency;coeff_updates];

        dof(i).parameter_a=dof(i).parameter_a+coeff_updates;
        aaa=dof(i).parameter_a;
        norm(temp_sensitivity_parameter_da)
    end
end

% time_data=0:0.01:150;
x=zeros(num_degrees_freedom,length(time_data));dx=zeros(num_degrees_freedom,length(time_data));ddx=zeros(num_degrees_freedom,length(time_data));
for ii=1:1:num_degrees_freedom
    temp_parameter_a=dof(ii).parameter_a;
    w0=temp_parameter_a(1,1);
    Harm_parameter_a=temp_parameter_a(3:end,:);

    %     w0=parameter_a(1,1);
    %     Harm_parameter_a=parameter_a(2:end,2*ii-1:2*ii);

    %% 计算方程残差
    for i=1:num_harmonics   % i=1,3,5
        x(ii,:)=x(ii,:)+Harm_parameter_a(i,1)*cos(i*w0*time_data)+Harm_parameter_a(i,2)*sin(i*w0*time_data);
        %         dx(ii,:)=dx(ii,:)-w0*(2*i-1)*Harm_parameter_a(i,1)*sin((2*i-1)*w0*time_data)+w0*(2*i-1)*Harm_parameter_a(i,2)*cos((2*i-1)*w0*time_data);
        %         ddx(ii,:)=ddx(ii,:)-(w0*(2*i-1))^2*Harm_parameter_a(i,1)*cos((2*i-1)*w0*time_data)-(w0*(2*i-1))^2*Harm_parameter_a(i,2)*sin((2*i-1)*w0*time_data);
    end
    x(ii,:)=x(ii,:)+temp_parameter_a(2, 1);
end

figure;
plot(time_data,x(1,:),'r-','LineWidth',1.5);
hold on;
plot(time_data,x(2,:),'b-','LineWidth',1.5);
% h1=legend('$$x_1$$','$$x_2$$','$$x_3$$','$$x_4$$','$$x_5$$');
% set(h1,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

hold on;
plot(time_data,num(temp1:temp2,1),'k-');
hold on
plot(time_data,num(temp1:temp2,2),'k-');
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);


parameter_a1=dof(1).parameter_a;
parameter_a2=dof(2).parameter_a;
harmonic_coefficients=zeros(num_harmonics+2,4);
harmonic_coefficients(1:num_harmonics+2,1:2)=parameter_a1;
harmonic_coefficients(2:num_harmonics+2,3:4)=parameter_a2(2:end,:);
save('fitting_para_alpha_0.8_alpha1_1_period_1.mat','harmonic_coefficients','-v7.3');
