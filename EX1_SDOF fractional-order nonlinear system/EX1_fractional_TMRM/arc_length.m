function Predictor=arc_length(every_a)
for j=1:3
    s(j)=sqrt(norm(every_a(j+1).harmonic_coefficients-every_a(j).harmonic_coefficients));
end
t0=0;t1=s(1);t2=t1+s(2);t3=t2+s(3);
%ds(j)=ds(j-1)*Nd/I(j-1),Nd为一个与非线性频率-振幅响应曲线相关的整数，一般取4或5,
%I(j-1)为完成上一次计算所进行的迭代次数  具体参见张丹伟博士论文P51页
ds=s(3)+(s(3)-s(2));
t4=t3+ds;t=[t0;t1;t2;t3;t4];
for i=1:4
    temp=0;
    for j=1:4
        if j~=i
            temp(j)=(t(5)-t(j))/(t(i)-t(j));
        end
    end
    temp(find(temp==0))=[];
    temp_every_a(i).harmonic_coefficients=temp(1)*temp(2)*temp(3)*every_a(i).harmonic_coefficients;
end

Predictor=zeros(size(every_a(1).harmonic_coefficients));
for i=1:4
    Predictor=Predictor+temp_every_a(i).harmonic_coefficients;
end
end