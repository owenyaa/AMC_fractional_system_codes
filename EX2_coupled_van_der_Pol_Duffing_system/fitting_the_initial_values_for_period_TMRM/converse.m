function y=converse(A)%............将由A表示的F级数的系数矩阵转换为列向量
tempn=length(A(:,1));
y=zeros(tempn,1);
y(1)=A(1,1);
for i=2:tempn
    y(2*i-2)=A(i,1);
    y(2*i-1)=A(i,2);
end