clear; clc; close all;
tic;
% 参数初始化
mu = 0.8; k = 3; epsilon = 1.2; alpha = 0.8; alpha1 = 1.5; alpha2 = 0.9; beta = 0.05;
beta2 = 1/4; r = 1/2;

tf = 150; h = 0.01; Tdata = (0:h:tf)';
iteration = length(Tdata);
x = zeros(2, iteration); dx = x; ddx = x;

% 初始条件
x(:, 1) = [1.09038923768306; 0.835475468411433];
dx(:, 1) = [-2.81381542670743e-15; 0.418646768041845];
ddx(:, 1) = [-1.85513054549973; -1.59460672128468];

% 预计算 c_coef1 和 c_coef2
c_coef1 = zeros(iteration, iteration);
c_coef2 = zeros(iteration, iteration);
for t = 1:iteration
    for nj = 1:t
        c_coef1(nj, t) = (t - nj + 2)^(2 - alpha1) + (t - nj)^(2 - alpha1) - 2 * (t - nj + 1)^(2 - alpha1);
        c_coef2(nj, t) = (t - nj + 2)^(2 - alpha2) + (t - nj)^(2 - alpha2) - 2 * (t - nj + 1)^(2 - alpha2);
    end
end

for t = 1:iteration
    % 计算 M 矩阵和 CONS1, CONS2
    temp_CONS1 = -r / (beta2 * h) * x(1, t) + (1 - r / beta2) * dx(1, t) + (1 - 0.5 * r / beta2) * h * ddx(1, t);
    temp_CONS2 = -r / (beta2 * h) * x(2, t) + (1 - r / beta2) * dx(2, t) + (1 - 0.5 * r / beta2) * h * ddx(2, t);
    
    m11 = 1 / (beta2 * h^2) - mu * r * h^(1 - alpha1) / (gamma(3 - alpha1) * beta2 * h) + (1 + k);
    m12 = -k;
    m21 = -k * r / (beta2 * h);
    m22 = 1 / (beta2 * h^2) + alpha * r * h^(1 - alpha2) / (gamma(3 - alpha2) * beta2 * h) - 1 + k * r / (beta2 * h);
    M = [m11, m12; m21, m22];
    
    CONS1 = -1 / (beta2 * h^2) * x(1, t) - 1 / (beta2 * h) * dx(1, t) - (0.5 / beta2 - 1) * ddx(1, t) - ...
        mu * h^(1 - alpha1) / gamma(3 - alpha1) * (dx(1, 1:t) * c_coef1(1:t, t))' + temp_CONS1;
    CONS2 = -1 / (beta2 * h^2) * x(2, t) - 1 / (beta2 * h) * dx(2, t) - (0.5 / beta2 - 1) * ddx(2, t) + ...
        alpha * h^(1 - alpha2) / gamma(3 - alpha2) * (dx(2, 1:t) * c_coef2(1:t, t))' + temp_CONS2 - k * temp_CONS1;
    
    % 非线性方程求解部分优化
    xn = x(:, t);
    for iteration_count = 1:10
        dxn = -df(xn(1), xn(2), M, mu, epsilon) \ f(xn(1), xn(2), M, mu, epsilon, CONS1, CONS2);
        xn = xn + dxn;
        if norm(dxn) < 1e-8
            break;
        end
    end
    
    x(:, t + 1) = xn;
    dx(:, t + 1) = r / (beta2 * h) * (x(:, t + 1) - x(:, t)) + (1 - r / beta2) * dx(:, t) + (1 - 0.5 * r / beta2) * h * ddx(:, t);
    ddx(:, t + 1) = 1 / (beta2 * h^2) * (x(:, t + 1) - x(:, t)) - 1 / (beta2 * h) * dx(:, t) - (0.5 / beta2 - 1) * ddx(:, t);
end
toc;
% 绘图
x_cal = [x; dx; ddx]';
figure;
subplot(2,1,1);
plot(Tdata, x_cal(1:end-1, 1), 'k-');
legend('$$TMRM-x_1$$', 'Interpreter', 'latex', 'FontSize', 15);
subplot(2,1,2);
plot(Tdata, x_cal(1:end-1, 2), 'k-');
legend('$$Newmark-\beta-x_1$$', 'Interpreter', 'latex', 'FontSize', 15);

figure;
plot(x_cal(7000:end-1, 1), x_cal(7000:end-1, 3), 'k-');
legend('$$Newmark-\beta$$', 'Interpreter', 'latex', 'FontSize', 15);

% 非线性方程组的函数
function y = f(x1, x2, M, mu, epsilon, CONS1, CONS2)
    y = [M(1,1)*x1 + M(1,2)*x2 + mu*x1^2 + epsilon*x1^3 + CONS1;
         M(2,1)*x1 + M(2,2)*x2 + mu*x2^2 + epsilon*x2^3 + CONS2];
end

function y = df(x1, x2, M, mu, epsilon)
    y = [M(1,1) + 2*mu*x1 + 3*epsilon*x1^2, M(1,2);
         M(2,1), M(2,2) + 2*mu*x2 + 3*epsilon*x2^2];
end
