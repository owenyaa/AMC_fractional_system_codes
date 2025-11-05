function y=time_to_fourier(w_input,t_input,x_input)%....tt是时间，xx是对于对应的时程，w是预设频率
global num_harmonics_fitting%....所考虑谐波个数

temp=size(t_input);
if temp(1)<=temp(2)
    t_input=t_input';
end

N=fix(3*num_harmonics_fitting);%........保留谐波数量
M=length(t_input);
tt=w_input*t_input;
%tt=tt/1.001;
%tt=M*rand(1,M)';

f_xy_t=x_input;
%tt=tt*1.001;
%xx=xx./(1+xx.^2);

Jacobi=[];
for i=1:N
    i_tt=i*tt;
    Jacobi=[Jacobi cos(i_tt) sin(i_tt)];
end
%....consider constant term 
Jacobi=[ones(M,1) Jacobi];

CS=Jacobi\f_xy_t;
y=anti_converse(CS);%.......对应于f_xy_t=(1-xt.^2-yt.^2).^(-n1);的展开式
y=y(1:(num_harmonics_fitting+1),:);

a=1;