function x_cal=cal_x(parameter_a)
global h Tdata
e=parameter_a(1);alpha1=parameter_a(2);alpha2=parameter_a(3);k=parameter_a(4);w0=parameter_a(5);
r=1/2;
beta=0.25;
% tf=100;h=0.01;Tdata=(0:h:tf)';
iteration=length(Tdata);
x=zeros(1,iteration);dx=x;ddx=x;
% x=zeros(4,iteration);dx=x;ddx=x;
% 强迫振动

% x(1,1)=0.00;dx(1,1)=0.1;ddx(1,1)=0.0;
x(1,1)=0.087927636265798;dx(1,1)=0.023959804540231;ddx(1,1)=-0.056565054025158;

%计算正问题的响应

for t=1:iteration
    for nj=1:t
        c_coef1(nj)=(t-nj+2)^(2-alpha1)+(t-nj)^(2-alpha1)-2*(t-nj+1)^(2-alpha1);
        c_coef2(nj)=(t-nj+2)^(2-alpha2)+(t-nj)^(2-alpha2)-2*(t-nj+1)^(2-alpha2);
    end
        coeff_1=r*h^(1-alpha1)/(gamma(3-alpha1)*beta*h)+e*r*h^(1-alpha2)/(gamma(3-alpha2)*beta*h)+w0;
        coeff_2=-e*h^(1-alpha2)/gamma(3-alpha2)*(dx(1,1:t)*c_coef2(1:t)'-r/(beta*h)*x(1,t)+(1-r/beta)*dx(1,t)+...
        (1-r/(2*beta))*h*ddx(1,t));
        coeff_3=-e*r*h^(1-alpha2)/(gamma(3-alpha2)*beta*h)+k;
        f=@(coeff_1,coeff_2,coeff_3,cons,xn)coeff_1*xn+coeff_2*xn^2+coeff_3*xn^3+cons;
        df=@(coeff_1,coeff_2,coeff_3,xn)coeff_1+2*coeff_2*xn+3*coeff_3*xn^2;
      
        cons=h^(1-alpha1)/gamma(3-alpha1)*(dx(1,1:t)*c_coef1(1:t)'-r/(beta*h)*x(1,t)+(1-r/beta)*dx(1,t)+...
        (1-r/(2*beta))*h*ddx(1,t))+...
             e*h^(1-alpha2)/gamma(3-alpha2)*(dx(1,1:t)*c_coef2(1:t)'-r/(beta*h)*x(1,t)+(1-r/beta)*dx(1,t)+...
        (1-r/(2*beta))*h*ddx(1,t));
        xn=x(1,t);dxn=0.1;
    while abs(dxn)>=1e-8
        dxn=-f(coeff_1,coeff_2,coeff_3,cons,xn)/df(coeff_1,coeff_2,coeff_3,xn);
        xn=xn+dxn;
    end
        x(1,t+1)=xn;
        dx(1,t+1)=r/(beta*h)*(x(1,t+1)-x(1,t))+(1-r/beta)*dx(1,t)+(1-0.5*r/beta)*h*ddx(1,t);
        ddx(1,t+1)=1/(beta*h^2)*(x(1,t+1)-x(1,t))-1/(beta*h)*dx(1,t)-(0.5/beta-1)*ddx(1,t);
end
x_cal=[x;dx;ddx];