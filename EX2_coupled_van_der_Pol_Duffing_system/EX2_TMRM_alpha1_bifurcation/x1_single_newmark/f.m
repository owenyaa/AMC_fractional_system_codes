function y = f(x1, x2, M, k11, k12, k21, k22, CONS1, CONS2)
    % 确保每个计算都是按元素操作
    y = zeros(2, 1);  % 确保 y 是一个列向量
    
    % 计算 y(1) 和 y(2)
    y(1,1) = M(1,1)*x1 + M(1,2)*x2 + k11*x1^2 + k12*x1^3 + CONS1;
    y(2,1) = M(2,1)*x1 + M(2,2)*x2 + k21*x2^2 + k22*x2^3 + CONS2;
end

