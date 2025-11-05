clear;clc;close all;
tic;
% 周期
mu=0.8;k=3;epsilon=1.2;alpha=0.8;alpha1=1.5;alpha2=0.9;beta=0.05;
% alpha1=1.5; alpha2=0.8; mu=0.8; k=3; alpha=0.8; beta=0.05; epsilon=1.2;
beta2=1/4;r=1/2;

tf=150;h=0.01;Tdata=(0:h:tf)';
iteration=length(Tdata);
x=zeros(2,iteration);dx=x;ddx=x;

% 周期-2
x(1,1)=1.09038923768306;x(2,1)=0.835475468411433;
dx(1,1)=-2.81381542670743e-15;dx(2,1)=0.418646768041845;
ddx(1,1)=-1.85513054549973;ddx(2,1)=-1.59460672128468;


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
toc;
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

function y=f(x1,x2,M,k11,k12,k21,k22,CONS1,CONS2)
y(1,1)=M(1,1)*x1+M(1,2)*x2+k11*x1^2+k12*x1^3+CONS1;
y(2,1)=M(2,1)*x1+M(2,2)*x2+k21*x2^2+k22*x2^3+CONS2;
end

function y=df(x1,x2,M,k11,k12,k21,k22)
y(1,1)=M(1,1)+2*k11*x1+3*k12*x1^2;
y(1,2)=M(1,2);
y(2,1)=M(2,1);
y(2,2)=M(2,2)+2*k21*x2+3*k22*x2^2;
end