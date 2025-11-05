clear;clc;close all;

% 周期
mu=0.8;k=3;epsilon=1.2;alpha=0.8;
% alpha1=1.5;
alpha2=0.9;beta=0.05;
% alpha1=1.5; alpha2=0.8; mu=0.8; k=3; alpha=0.8; beta=0.05; epsilon=1.2;
beta2=1/4;r=1/2;

tf=500;h=0.01;Tdata=(0:h:tf)';
iteration=length(Tdata);
x=zeros(2,iteration);dx=x;ddx=x;

% x(1,1)=0;dx(1,1)=0.5;ddx(1,1)=0;
% x(2,1)=0;dx(2,1)=0.3;ddx(2,1)=0;


load the_two_parameter_alpha_0.8_alpha1_1_1.033.mat;
num_degrees_freedom = 2;
num_harmonics = 25;
index_bifurcation=22;
index_bifurcation1=1;
% for alpha1=1.429:0.001:1.617
for alpha1=1.021:0.001:1.031
    alpha1
    %alpha=every_harmonic_coefficients(index_bifurcation).alpha;
    harmonic_coefficients=every_harmonic_coefficients(index_bifurcation).harmonic_coefficients;
    index_residual=every_harmonic_coefficients(index_bifurcation).index_residual;

    initial_frequency = harmonic_coefficients(1, 1);
    total_time = 2 * pi / initial_frequency;
    step_size = 2 * pi / (initial_frequency * 500);
    time_data = (0:step_size:total_time);
    harmonic_parameters = harmonic_coefficients(3:end, :);

    % Initialize arrays for displacement, velocity, and acceleration
    displacement = zeros(num_degrees_freedom, length(time_data));
    velocity = zeros(num_degrees_freedom, length(time_data));
    acceleration = zeros(num_degrees_freedom, length(time_data));

    % Compute displacement, velocity, and acceleration for each degree of freedom
    for dof_index = 1:num_degrees_freedom
        for harm_index = 1:num_harmonics
            harmonic_term = harm_index  * initial_frequency * time_data;
            displacement(dof_index, :) = displacement(dof_index, :) + harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(harmonic_term) + harmonic_parameters(harm_index, 2 * dof_index) * sin(harmonic_term);
            velocity(dof_index, :) = velocity(dof_index, :) - initial_frequency * harm_index  * harmonic_parameters(harm_index, 2 * dof_index - 1) * sin(harmonic_term) + initial_frequency * harm_index  * harmonic_parameters(harm_index, 2 * dof_index) * cos(harmonic_term);
            acceleration(dof_index, :) = acceleration(dof_index, :) - (initial_frequency * harm_index )^2 * harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(harmonic_term) - (initial_frequency * harm_index )^2 * harmonic_parameters(harm_index, 2 * dof_index) * sin(harmonic_term);
        end
        displacement(dof_index, :) = displacement(dof_index, :) + harmonic_coefficients(2, 2*dof_index-1);
    end

    x(1,1)=displacement(1,1);dx(1,1)=velocity(1,1);ddx(1,1)=acceleration(1,1);
    x(2,1)=displacement(2,1);dx(2,1)=velocity(2,1);ddx(2,1)=acceleration(2,1);
    % 周期-2

    for t=1:iteration
        for nj=1:t
            %         c_coef(nj)=(t-nj+2)^(2-alpha)+(t-nj)^(2-alpha)-2*(t-nj+1)^(2-alpha);
            c_coef1(nj)=(t-nj+2)^(2-alpha1)+(t-nj)^(2-alpha1)-2*(t-nj+1)^(2-alpha1);
            c_coef2(nj)=(t-nj+2)^(2-alpha2)+(t-nj)^(2-alpha2)-2*(t-nj+1)^(2-alpha2);
        end
        temp_CONS1=-r/(beta2*h)*x(1,t)+(1-r/beta2)*dx(1,t)+(1-0.5*r/beta2)*h*ddx(1,t);
        temp_CONS2=-r/(beta2*h)*x(2,t)+(1-r/beta2)*dx(2,t)+(1-0.5*r/beta2)*h*ddx(2,t);
        m11=1/(beta2*h^2)-mu*r*h^(1-alpha1)/(gamma(3-alpha1)*beta2*h)+(1+k);
        m12=-k;
        m21=-k*r/(beta2*h);
        m22=1/(beta2*h^2)+alpha*r*h^(1-alpha2)/(gamma(3-alpha2)*beta2*h)-1+k*r/(beta2*h);
        M=[m11,m12;m21,m22];
        CONS1=-1/(beta2*h^2)*x(1,t)-1/(beta2*h)*dx(1,t)-(0.5/beta2-1)*ddx(1,t)-...
            mu*h^(1-alpha1)/gamma(3-alpha1)*(dx(1,1:t)*c_coef1(1:t)'-r/(beta2*h)*x(1,t)+(1-r/beta2)*dx(1,t)+...
            (1-r/(2*beta2))*h*ddx(1,t));
        CONS2=-(1/(beta2*h^2)*x(2,t)+1/(beta2*h)*dx(2,t)+(0.5/beta2-1)*ddx(2,t))+...
            alpha*h^(1-alpha2)/gamma(3-alpha2)*(dx(2,1:t)*c_coef2(1:t)'-r/(beta2*h)*x(2,t)+(1-r/beta2)*dx(2,t)+...
            (1-r/(2*beta2))*h*ddx(2,t))+...
            k*temp_CONS2-k*temp_CONS1;
        k11=mu*h^(1-alpha1)/gamma(3-alpha1)*(dx(1,1:t)*c_coef1(1:t)'-r/(beta2*h)*x(1,t)+(1-r/beta2)*dx(1,t)+...
            (1-r/(2*beta2))*h*ddx(1,t));
        k12=mu*r*h^(1-alpha1)/(gamma(3-alpha1)*beta2*h);
        k21=beta;
        k22=epsilon;

        xn=x(1:2,t);dxn=[0.1;0.1];
        while abs(dxn)>=1e-8
            dxn=-df(xn(1),xn(2),M,k11,k12,k21,k22)\f(xn(1),xn(2),M,k11,k12,k21,k22,CONS1,CONS2);
            xn=xn+dxn;
        end
        x(1,t+1)=xn(1);
        dx(1,t+1)=r/(beta2*h)*(x(1,t+1)-x(1,t))+(1-r/beta2)*dx(1,t)+(1-0.5*r/beta2)*h*ddx(1,t);
        ddx(1,t+1)=1/(beta2*h^2)*(x(1,t+1)-x(1,t))-1/(beta2*h)*dx(1,t)-(0.5/beta2-1)*ddx(1,t);
        x(2,t+1)=xn(2);
        dx(2,t+1)=r/(beta2*h)*(x(2,t+1)-x(2,t))+(1-r/beta2)*dx(2,t)+(1-0.5*r/beta2)*h*ddx(2,t);
        ddx(2,t+1)=1/(beta2*h^2)*(x(2,t+1)-x(2,t))-1/(beta2*h)*dx(2,t)-(0.5/beta2-1)*ddx(2,t);
    end
    x_cal=[x;dx;ddx];
    x_cal=x_cal';

    num_xcal(index_bifurcation1).x_cal=x_cal;
    num_xcal(index_bifurcation1).alpha1=alpha1;
    index_bifurcation=index_bifurcation+1;
    index_bifurcation1=index_bifurcation1+1;
end
save('x1_num_xcal_1.021_1.031_red.mat', 'num_xcal','-v7.3');









x_cal=[x;dx;ddx];
x_cal=x_cal';
% x=x_cal(:,1:2);
%xmin(jj).min=getmin(x);
mm=1;
figure;
subplot(2,1,1);
plot(Tdata(1:mm:end),x_cal(1:mm:end-1,1),'k-','MarkerSize',10);
h1=legend('$$TMRM-x_1$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

subplot(2,1,2);
plot(Tdata(1:mm:end),x_cal(1:mm:end-1,2),'k-','MarkerSize',10);
h1=legend('$$Newmark-\beta-x_1$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
mm=1;
figure;
plot(x_cal(7000:mm:end-1,1),x_cal(7000:mm:end-1,3),'k-','MarkerSize',6);
h1=legend('$$Newmark-\beta$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% figure;
% plot(displacement(1,:),velocity(1,:),'k-','LineWidth',1.5);
% % hold on;
% % plot(displacement(2,:),velocity(2,:),'k-','LineWidth',1.5);
% % set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
%
% % tf=300;h=0.01;Tdata=(0:h:tf)';
% % initial_frequency=1.05461;
% Tdata1=(tf-2*pi/initial_frequency:h:tf)';
% Tdata2=length(Tdata1);
% % Tdata=0:0.04:500;
% x=x_cal(:,1:2);dx=x_cal(:,3:4);
% x=x';dx=dx';
% figure;
% plot(x(1,length(Tdata)-Tdata2:20:length(Tdata)-Tdata2+length(Tdata1)),dx(1,length(Tdata)-Tdata2:20:length(Tdata)-Tdata2+length(Tdata1)),'r.','MarkerSize',10);
% hold on;
% plot(x(2,length(Tdata)-Tdata2:20:length(Tdata)-Tdata2+length(Tdata1)),dx(2,length(Tdata)-Tdata2:20:length(Tdata)-Tdata2+length(Tdata1)),'b.','MarkerSize',10);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
%
%
% figure;
% plot(Tdata,x(1,1:end-1),'k-','LineWidth',1.5);
% hold on;
% plot(Tdata,x(2,1:end-1),'k-','LineWidth',1.5);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);





% x_cal=[x;dx;ddx];
% x_cal=x_cal';
% x=x_cal(:,1);dx=x_cal(:,2);ddx=x_cal(:,3);
% x_cal=[x,dx,ddx];
% savefile='simple_fre_datanewa.mat';
% save(savefile,'Tdata','x_cal');
% figure
% plot(Tdata,x_cal(1:end-1,1),'k-','LineWidth',1);
% hold on;
% % plot(Tdata,x_cal(1:end-1,2),'g-','LineWidth',1);
% % hold on;
% % plot(Tdata,x_cal(1:end-1,3),'r-','LineWidth',1);
% h=legend('浣嶇Щ');
% set(h,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','瀹嬩綋','FontSize',15,'LineWidth',1.5);
% %鐩稿浘;
% figure
% plot(x_cal(2000:end-1,1),x_cal(2000:end-1,2),'k-','LineWidth',1);
% figure
% plot(FF,dis,'.');
% title('鍙傛暟鍙樺寲涓嬬殑鍒嗗矓鍥?')
% xlabel('澶栨縺鍔盕鐨勫彉鍖?')
% ylabel('X鍊?')