function y=anti_converse(A)%............将由A表示的F级数的系数矩阵转换为列向量
tempn=length(A);
tempm=(tempn-1)/2;
y=zeros(tempm,2);
y(1,1)=A(1);
for i=2:(tempm+1)
    y(i,1)=A(2*i-2);
    y(i,2)=A(2*i-1);
end