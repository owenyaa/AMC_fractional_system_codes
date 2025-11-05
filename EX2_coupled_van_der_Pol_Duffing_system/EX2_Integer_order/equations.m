% 定义运动方程
function yprime = equations(t, y)
    global mu beta epsilon A mass_matrix
    % 外部激励项
    F=[0;0;-mass_matrix\[mu*y(1)^2*y(3);beta*y(2)^2+epsilon*y(2)^3]];
    yprime=A*y+F;
end