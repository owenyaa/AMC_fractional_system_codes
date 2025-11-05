function x_cal=cal_x_even(parameter_a)
global h Tdata
e=parameter_a(1);alpha=parameter_a(2);w0=parameter_a(3);k=parameter_a(4);para_constant=parameter_a(5);
r=1/2;beta=0.25;

iteration=length(Tdata);
x=zeros(1,iteration);dx=x;ddx=x;
% x=zeros(4,iteration);dx=x;ddx=x;
% 强迫振动
x(1,1)=1.128377299326538e+00;dx(1,1)=-2.967738146599960e+00;ddx(1,1)=-2.480015655621834e+00;
% x(1,1)=0.653167723305119;dx(1,1)=1.529692797995123;ddx(1,1)=0.295589656889249;
% x(1,1)=1.039311080139804;dx(1,1)=3.847953442625482;ddx(1,1)=-2.340006622697429;%k=1;
  
% x(1,1)=0.385142584552087;dx(1,1)=2.662372260455271;ddx(1,1)=0.855104709949905;%k=0;

% x(1,1)=1.062186294003650;dx(1,1)=2.656033565281555;ddx(1,1)=-1.300233690272597;
%计算正问题的响应


for t=1:iteration
    for nj=1:t
        c_coef(nj)=(t-nj+2)^(2-alpha)+(t-nj)^(2-alpha)-2*(t-nj+1)^(2-alpha);
    end
        coeff_1=1/(beta*h^2)+e*r*h^(1-alpha)/(gamma(3-alpha)*beta*h)+w0;
        coeff_2= -e*para_constant*h^(1-alpha)/gamma(3-alpha)*(dx(1,1:t)*c_coef(1:t)'-r/(beta*h)*x(1,t)+(1-r/beta)*dx(1,t)+...
        (1-r/(2*beta))*h*ddx(1,t));
        coeff_3=-e*para_constant*r*h^(1-alpha)/(gamma(3-alpha)*beta*h)+k;
        f=@(coeff_1,coeff_2,coeff_3,cons,xn)coeff_1*xn+coeff_2*xn^2+coeff_3*xn^3+cons;
        df=@(coeff_1,coeff_2,coeff_3,xn)coeff_1+2*coeff_2*xn+3*coeff_3*xn^2;
      
        cons=-1/(beta*h^2)*x(1,t)-1/(beta*h)*dx(1,t)-(0.5/beta-1)*ddx(1,t)+...
        e*h^(1-alpha)/gamma(3-alpha)*(dx(1,1:t)*c_coef(1:t)'-r/(beta*h)*x(1,t)+(1-r/beta)*dx(1,t)+...
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