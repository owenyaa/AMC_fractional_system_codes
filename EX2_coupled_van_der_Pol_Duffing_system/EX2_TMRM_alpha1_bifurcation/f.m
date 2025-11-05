function y=f(x1,x2,M,k11,k12,k21,k22,CONS1,CONS2)
y(1,1)=M(1,1)*x1+M(1,2)*x2+k11*x1^2+k12*x1^3+CONS1;
y(2,1)=M(2,1)*x1+M(2,2)*x2+k21*x2^2+k22*x2^3+CONS2;
end