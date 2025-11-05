clear; % 清除工作空间的所有变量
close all; % 关闭所有图形窗口

% 加载存储系统周期解的文件
load period_solution_0.8_1.2.mat;

% 定义全局变量，这些变量将在多个函数中使用
global A stiffness_coeff mass_matrix epsilon nonlinearity_coeff num_dof num_harmonics

% 设置系统参数
num_dof = 1; % 单自由度系统
num_harmonics = 10; % 使用的谐波数量
epsilon = 0.1; % 非线性阻尼的系数
stiffness_coeff = -1; % 线性刚度系数 k = -1
nonlinearity_coeff = 1; % 非线性刚度系数 k_n = 1

multipliers_matrix = []; % 用于存储Floquet乘子的矩阵

% 定义质量矩阵和系统矩阵
mass_matrix = 1; % 单位质量
A = [0, 1; -stiffness_coeff, 0]; % 初始线性系统矩阵

% 遍历不同频率的参数集，计算Floquet乘子
for j = 1:81
    % 获取当前频率和对应的参数集
    frequency = every_a(j).w;
    parameter_a = every_a(j).parameter_a;

    % 计算系统的周期
    period = 2 * pi / frequency;

    % 初始化状态矩阵和状态转移矩阵
    dx0 = eye(2); % 2x2单位矩阵
    dx_t = zeros(2, 2);

    % 定义时间步长
    time_steps = 0:period / 10000:period;

    % 对每个初始条件进行数值积分以计算Floquet乘子
    for i = 1:2
        [t, num] = ode45(@ode_floquet, time_steps, dx0(:, i));
        temp_dx = num(end, :);
        dx_t(:, i) = temp_dx';
    end

    % 计算Floquet乘子
    floquet_multipliers = eig(dx_t);
    disp(floquet_multipliers); % 输出Floquet乘子

    % 绘制Floquet乘子在复平面上的分布
    plot(real(floquet_multipliers), imag(floquet_multipliers), 'k.');
    hold on;

    % 存储Floquet乘子
    multipliers_matrix = [multipliers_matrix floquet_multipliers];
end

% 绘制单位圆用于参考
theta = 0:0.01:2 * pi;
radius = ones(size(theta));
polar(theta, radius);
axis equal;

% 用于计算Floquet乘子的子函数
function yprime = ode_floquet(t, y)
    global A epsilon nonlinearity_coeff parameter_a mass_matrix num_harmonics frequency

    % 初始化位移和速度
    displacement = 0;
    velocity = 0;

    % 使用傅里叶级数展开计算位移和速度
    for i = 1:num_harmonics
        displacement = displacement + parameter_a(i, 1) * cos(i * frequency * t) + parameter_a(i, 2) * sin(i * frequency * t);
        velocity = velocity - i * frequency * parameter_a(i, 1) * sin(i * frequency * t) + i * frequency * parameter_a(i, 2) * cos(i * frequency * t);
    end

    % 计算非线性项的偏导数
    nonlinear_damping = epsilon * (1 - 1.1 * displacement^2); % 非线性阻尼项的系数
    nonlinear_stiffness = 3 * nonlinearity_coeff * displacement^2; % 非线性刚度项的线性化形式

    % 构建状态空间矩阵的非线性部分
    nonlinear_matrix = [0, 0; -nonlinear_stiffness, -nonlinear_damping];

    % 计算状态导数
    yprime = (A + nonlinear_matrix) * y;
end
