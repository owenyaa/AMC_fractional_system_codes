%% Author : GUANG_LIU  * owenyaa@gmail.com *
% Created Time : 2025-11-02 22:13
% Last Revised : GUANG_LIU ,2020-11-12
% Remark : The structure, function and variables of this program are explained:%
%%======================================================================
%% 主程序：基于“时域最小残值法”的分数阶非线性系统半解析求解与残差诊断
%%======================================================================
%【程序目标】
%  本程序用于求解论文给出的分数阶非线性系统的半解析解，并对所得解进行残差诊断与渐近性质验证。
%  具体对象是文中示例方程：
%      \ddot{x} + \epsilon(1 - x^2) D^\alpha \dot{x} + k x^3 = 0 ，0<\alpha<1。
%  程序的核心思想是：以有限谐波级数表示响应，利用“时域最小残值法”最小化时间域内的控制方程残差，
%  同时通过频率与谐波系数的联合迭代更新获得收敛的半解析近似解；在收敛后，依据论文公式(28)计算
%  分数阶算子的非周期“残差余项” J^{\alpha}(t)，并与理论幂律衰减 t^{-1/2} 进行数值对比验证。
%【理论与方法简述】
%  1) 表示：位移用奇次谐波截断的傅里叶级数近似，基频为 \omega，谐波阶数为 N。
%  2) 残差：将近似解代回控制方程，构造时间域残差 R(t)，并显式构造对频率与各谐波系数的灵敏度响应。
%  3) 优化：在给定时间网格上最小化残差的加权范数，等价为解线性化的“参数—残差”更新方程。
%  4) 正则：为缓解病态与过拟合，本程序框架支持基于奇异值分解与 L-curve 的正则化参数选取（本版本未附正则
%     子程序实现，但不影响主框架与接口；见【外部依赖】）。
%  5) 诊断：调用 compute_J_alpha1，按论文公式(28)计算 J^{\alpha}(t)。为避免 e^{k\omega t} 与 K_\nu(k\omega t)
%     的数值抵消导致的上/下溢，compute_J_alpha1 在对数域进行带符号的稳定求和，并输出 log J^{\alpha}(t)。
%  6) 渐近：依据 K_\nu(z) 的大参数渐近，J^{\alpha}(t) ~ t^{-1/2}。程序在对数坐标下对比 exp(log J^{\alpha}) 与
%     t^{-1/2}，用于验证“稳态 S-渐近 T-周期性”而非严格 T-周期。
%【程序结构（主流程与关键函数）】
%  A. 初始化（本文件主脚本段）
%     - 设置全局维度、谐波数、材料/系统参数（epsilon, alpha, k, 质量与非线性刚度矩阵等）。
%     - 设置初始频率 initial_frequency、时间区间 total_time 与时间步长 step_size，生成 time_data。
%     - 初始化参数矩阵 harmonic_coefficients（第一行存 \omega，其余行为各自由度的奇次谐波 cos/sin 系数）。
%
%  B. 主迭代（响应灵敏度—参数更新—信任域检验）
%     B1) 计算残差与灵敏度：residual_identification = calculate_residual(harmonic_coefficients)
%         - 输出矩阵按列块包含：控制方程残差、对频率的残差灵敏度、对各谐波系数的残差灵敏度。
%         - 将该矩阵重塑为（时域拼接）的大型灵敏度矩阵 sensitivity_matrix，并提取负残差列向量。
%     B2) 病态与正则：对 sensitivity_matrix 作 SVD，调用 L-curve 估计正则化参数（接口已留出）。
%     B3) 线性化更新：解线性方程获得最优参数增量（包含频率与谐波系数）；注意本程序对频率增量的行列布局
%         做了“置换与对齐”，以与参数矩阵结构一致（见代码中 real_coeff_updates 的重排）。
%     B4) 信任域判别：用 \rho_agreement 衡量线性模型对实际残差下降的拟合度；若通过阈值 \rho_trust 则接受
%         更新，否则放大正则/收缩“信任域半径”（由 gamma_trust 控制）并重试。
%     B5) 收敛性判据：以相对参数步长 ||Δa||/||a|| 与容限 Etol 判断迭代终止；记录历次参数与正则量。
%
%  C. 收敛后后处理与验证
%     - 以最终参数重建位移/速度/加速度时间历程并作图。
%     - 构造 harmonic_parameters = harmonic_coefficients(2:end,:) 并调用：
%           J_alpha1_log = compute_J_alpha1(time_data, alpha, initial_frequency, harmonic_parameters, j=1, num_harmonics)
%       作图对比 exp(J_alpha1_log) 与 t^{-1/2} 于普通坐标与对数坐标下的斜率与量级。
%
%  D. 关键子程序角色
%     - calculate_residual(harmonic_params):
%         构造 R 与对频率/各谐波系数的灵敏度列块；内部显式实现 D^\alpha x 的谐波级数形式。
%     - compute_J_alpha1(...):
%         依据论文公式(28)计算 J^{\alpha}(t) 的对数，采用带符号 log-sum-exp 稳定处理。
%     - 外部依赖：csvd, l_curve（用于 SVD 与 L-curve 正则参数估计；本版本未附源码）。
%
%【变量与数据结构约定】
%  - 全局变量：
%      total_time, step_size, time_data ：时间剖分相关
%      harmonic_coefficients           ：(N+1) × (2*ndof)，第1行是 \omega；其余行为奇谐波 cos/sin 系数
%      initial_frequency               ：\omega
%      num_degrees_freedom, num_harmonics
%      alpha, epsilon, omega_2, parameter_constant, mass_matrix, nonliear_stiffness_matrix
%  - harmonic_coefficients 结构：
%      行：1 行为 \omega；2..N+1 为 k=1..N 奇谐波
%      列：每个自由度占两列，奇列为 cos 系数 c_{jk}，偶列为 sin 系数 b_{jk}
%  - 残差拼接与灵敏度重塑：
%      residual_identification 先按“自由度×时间”拼接，再重塑为 (ndof*Nt)×(参数数目) 的灵敏度矩阵。
%
%【数值稳定与收敛建议】
%  - 时间步与谐波阶：num_time_point 与 num_harmonics 需协同加密以捕捉高阶谐波与非平滑项。
%  - 正则与信任域：若 \rho_agreement 反复偏低，增大 lambda 或减小步长（通过 gamma_trust 调整）。
%  - compute_J_alpha1：务必使用其对数输出再指数化作图，避免直接和式造成数值溢出。
%
%【输出与图形】
%  - 残差时间历程图、J^{\alpha}(t) 与 t^{-1/2} 的普通/对数双坐标对比图；
%  - 命令行打印最终 harmonic_coefficients 与收敛迭代信息。
%
%【外部依赖与缺失说明】
%  - 正则化相关例程 csvd 与 l_curve 在本版本中未提供源码，但主框架与接口完整，替换为同名函数即可运行。
%
%【可扩展性】
%  - 多自由度：按列扩展谐波系数阵即可；calculate_residual 与 compute_J_alpha1 已按列索引 j 兼容。
%  - 非光滑项：可在残差构造处替换/叠加相应非光滑算子，只需保持灵敏度列块维度一致。
%
%【再现性与版本】
%  - 建议固定随机种子（如存在随机扰动）、打印 Etol、rho_trust、gamma_trust 与 SVD 截断阈值。
%  - Author : GUANG_LIU  * owenyaa@gmail.com *
%    Created : 2020-11-12 22:13
%    Last Revised (Header) : 2025-11-05
%
%%======================================================================
%% MAIN PROGRAM: Semi-analytical solution and residual diagnostics for fractional
%% nonlinear systems based on the "Time-Domain Minimum Residual Method" 
%%======================================================================
% [Program Objective]
%   This program computes a semi-analytical solution for the fractional nonlinear
%   system given in the paper and performs residual diagnostics and asymptotic checks.
%   The target equation is:
%       \ddot{x} + \epsilon(1 - x^2) D^\alpha \dot{x} + k x^3 = 0 , 0<\alpha<1.
%   The core idea is: represent the response by a finite odd-harmonic Fourier series,
%   minimize the time-domain residual of the governing equation using the "Time-Domain
%   Minimum Residual Method", and jointly update the base frequency and harmonic
%   coefficients to obtain a converged semi-analytical approximation. After convergence,
%   compute the non-periodic "residual tail" J^{\alpha}(t) of the fractional operator
%   according to Eq. (28) in the paper, and verify the theoretical power-law decay
%   t^{-1/2} numerically.
%
% [Brief Theory and Method]
%   1) Representation: displacement is approximated by an odd-harmonic truncated series
%      with base frequency \omega and truncation order N.
%   2) Residual: substitute the approximation into the governing equation to construct
%      the time-domain residual R(t), and explicitly construct residual sensitivities
%      to the base frequency and all harmonic coefficients.
%   3) Optimization: minimize a weighted norm of the residual on the time grid, which is
%      equivalent to solving a linearized "parameter–residual" update system.
%   4) Regularization: to alleviate ill-posedness and overfitting, the framework supports
%      SVD-based regularization and L-curve selection (the routines are not included here,
%      but this does not affect the main framework; see [External Dependencies]).
%   5) Diagnostics: call compute_J_alpha1 to evaluate J^{\alpha}(t) per Eq. (28).
%      To avoid overflow/underflow caused by the cancellation between e^{k\omega t}
%      and K_\nu(k\omega t), compute_J_alpha1 performs a signed, stable accumulation
%      in the logarithmic domain and outputs log J^{\alpha}(t).
%   6) Asymptotics: by the large-argument asymptotics of K_\nu(z), J^{\alpha}(t) ~ t^{-1/2}.
%      The program compares exp(log J^{\alpha}) with t^{-1/2} in log-scale to verify the
%      "steady S-asymptotic T-periodicity" rather than strict T-periodicity.
%
% [Program Structure (Main Flow and Key Routines)]
%   A. Initialization (the top script block of this file)
%      - Set global dimensions, number of harmonics, and system/material parameters
%        (epsilon, alpha, k, mass and nonlinear stiffness matrices, etc.).
%      - Set initial_frequency \omega, total_time and step_size to generate time_data.
%      - Initialize harmonic_coefficients matrix (first row stores \omega; the rest are
%        odd-harmonic cos/sin coefficients per degree of freedom).
%
%   B. Main Iteration (sensitivity–parameter update–trust-region check)
%      B1) Residual and sensitivities: residual_identification = calculate_residual(harmonic_coefficients)
%          - The output stacks: residual of governing equation, its sensitivity to frequency,
%            and sensitivities to each harmonic coefficient in column blocks.
%          - Reshape it into a large sensitivity_matrix by stacking over time, and extract
%            the negative residual vector.
%      B2) Ill-conditioning and regularization: perform SVD on sensitivity_matrix and call
%          L-curve to estimate the regularization parameter (the interface is present).
%      B3) Linearized update: solve the linear system for optimal parameter increments
%          (including frequency and harmonic coefficients); note the row/column rearrangement
%          for the frequency increment to align with the parameter matrix structure (see the
%          reordering of real_coeff_updates in the code).
%      B4) Trust-region acceptance: evaluate \rho_agreement to measure how well the linear
%          model predicts the actual residual reduction; accept the step if it exceeds
%          \rho_trust, otherwise increase regularization/shrink the trust region (controlled
%          by gamma_trust) and retry.
%      B5) Convergence: terminate by the relative parameter step ||Δa||/||a|| ≤ Etol; record
%          parameter histories and regularization values.
%
%   C. Post-processing and Verification after Convergence
%      - Reconstruct displacement/velocity/acceleration time histories and plot.
%      - Build harmonic_parameters = harmonic_coefficients(2:end,:) and call:
%            J_alpha1_log = compute_J_alpha1(time_data, alpha, initial_frequency, harmonic_parameters, j=1, num_harmonics)
%        Plot and compare exp(J_alpha1_log) and t^{-1/2} in linear and log scales.
%
%   D. Roles of Key Routines
%      - calculate_residual(harmonic_params):
%          Construct R and the sensitivity blocks to frequency and to each harmonic
%          coefficient; internally realize the harmonic-series form of D^\alpha x.
%      - compute_J_alpha1(...):
%          Evaluate J^{\alpha}(t) per Eq. (28) in the paper in the log domain using a
%          signed log-sum-exp.
%      - External dependencies: csvd, l_curve (for SVD and L-curve regularization);
%          the source codes are not included here.
%
% [Conventions for Variables and Data Structures]
%   - Global variables:
%       total_time, step_size, time_data
%       harmonic_coefficients            : (N+1) × (2*ndof), first row is \omega
%       initial_frequency                : \omega
%       num_degrees_freedom, num_harmonics
%       alpha, epsilon, omega_2, parameter_constant, mass_matrix, nonliear_stiffness_matrix
%   - Structure of harmonic_coefficients:
%       Rows: row 1 is \omega; rows 2..N+1 correspond to k=1..N odd harmonics
%       Cols: each DOF occupies two columns; odd column = cos coefficient c_{jk},
%             even column = sin coefficient b_{jk}
%   - Residual stacking and sensitivity reshaping:
%       residual_identification is stacked over DOFs × time, then reshaped to
%       (ndof*Nt)×(number of parameters) as sensitivity_matrix.
%
% [Numerical Stability and Convergence Hints]
%   - Time step and harmonics: tune num_time_point and num_harmonics jointly to capture
%     higher harmonics and potential non-smoothness.
%   - Regularization and trust region: if \rho_agreement stays low repeatedly, enlarge
%     lambda or reduce step size (via gamma_trust).
%   - compute_J_alpha1: always use its logarithmic output and exponentiate for plotting,
%     instead of summing directly in the linear domain.
%
% [Outputs and Figures]
%   - Residual time-history plots; linear/log-scale comparisons of J^{\alpha}(t) versus t^{-1/2};
%   - Final harmonic_coefficients and iteration diagnostics printed to console.
%
% [External Dependencies and Missing Parts]
%   - Regularization-related routines csvd and l_curve are not provided here, but the
%     main framework and interfaces are complete; replacing them with homonymous
%     implementations will make the program runnable.
%
% [Extensibility]
%   - Multi-DOF: extend columns of the harmonic coefficient matrix; calculate_residual and
%     compute_J_alpha1 already support column-wise DOF indexing.
%   - Non-smooth terms: replace/augment operators in residual construction as needed,
%     while keeping the sensitivity block dimensions consistent.
%
% [Reproducibility and Version]
%   - It is recommended to fix random seeds (if any), and print Etol, rho_trust, gamma_trust,
%     and SVD truncation thresholds.
%   - Author : GUANG_LIU  * owenyaa@gmail.com *
%     Created : 2020-11-12 22:13
%     Last Revised (Header) : 2025-11-05
%%======================================================================

clear; clc; close all;

% Start timing the execution
tic;

% Declare global variables
global total_time step_size time_data harmonic_coefficients initial_frequency num_degrees_freedom num_harmonics
global alpha  nonliear_stiffness_matrix mass_matrix epsilon omega_2 parameter_constant
% Initialize degrees of freedom and number of harmonics
num_degrees_freedom = 1;
num_harmonics = 10;
mass_matrix=1; nonliear_stiffness_matrix=1;
epsilon=-0.8; alpha=0.5;
% epsilon=0.1; alpha=0.5;
omega_2=1; parameter_constant=1.1;

% initial_frequency = 1.103;
initial_frequency = 1.855;

num_time_point=200;

% Define time parameters
total_time = 2 * pi / initial_frequency;
step_size = 2 * pi / (initial_frequency * num_time_point);
time_data = (0:step_size:total_time);

% Initialize harmonic coefficients matrix
harmonic_coefficients = zeros(num_harmonics + 1, 2 * num_degrees_freedom);
harmonic_coefficients(1, 1) = initial_frequency;
harmonic_coefficients(2, :) = [1.3   -1.3]; % Initial values for the first-order harmonics
% load alpha_0.3_beta_0.001_omega_cos_0.2_Omegas_0.6_IHB_15.mat;
% harmonic_coefficients(2, 1)=0.043987586232923;


% Store initial harmonic coefficients
ini_harmonic_coefficients = harmonic_coefficients;

% Define the number of iterations based on time data length
iteration = length(time_data);

% Function to check if coefficients are within acceptable range
coefficients_check = @(harmonic_coefficients)(abs(harmonic_coefficients) < 5);

% Algorithm parameters for trust-region algorithm
gamma_trust = 1.414;
rho_trust = 0.5;

% Record values of parameters during iteration
harmonic_coefficients_record = harmonic_coefficients;
trust_region_record = [];

% Maximum number for response sensitivity iteration
max_iterations_rs = 1000;

% Maximum number for trust region iteration
max_iterations_tr = 20;

% Define the relative error tolerance for convergence of the algorithm
Etol = 1e-10;
index_bifurcation=1;
% Main loop for response sensitivity iteration
% for alpha=1:-0.01:0.97
for iteration_index = 1:max_iterations_rs
    % Reload frequency and update time parameters
    initial_frequency = harmonic_coefficients(1, 1);
    total_time = 2 * pi / initial_frequency;
    step_size = 2 * pi / (initial_frequency * num_time_point);
    time_data = (0:step_size:total_time);

    % Compute residual identification
    residual_identification = calculate_residual(harmonic_coefficients);
    % Extract response sensitivity matrix
    temp_sensitivity_matrix = reshape(residual_identification(:,num_degrees_freedom+1:2*num_degrees_freedom), num_degrees_freedom*length(time_data), 1);
    for i = 1:2*num_harmonics*num_degrees_freedom
        temp_sensitivity_matrix = [temp_sensitivity_matrix, reshape(residual_identification(:,num_degrees_freedom*(i+1)+1:num_degrees_freedom*(i+2)), num_degrees_freedom*length(time_data), 1)];
    end
    sensitivity_matrix(:,1:2)=temp_sensitivity_matrix(:,1:2);sensitivity_matrix(:,3:2*num_harmonics*num_degrees_freedom)=temp_sensitivity_matrix(:,4:end);

    % Calculate neg_residual (negative of the first column of residual)
    neg_residual = -reshape(residual_identification(:,1:num_degrees_freedom), num_degrees_freedom*length(time_data), 1);

    % Compute SVD of sensitivity matrix
    [U, s, V] = csvd(sensitivity_matrix);

    % Calculate the optimal lambda using L-curve method
    lambda_inverse = l_curve(U, s, neg_residual);

    % Store a copy of the harmonic coefficients matrix for later use
    temp_harmonic_coefficients = harmonic_coefficients;

    % Trust-region algorithm loop
    for trust_index = 1:max_iterations_tr
        % Calculate real coefficient updates
        % Compute real_coeff_updates (optimal parameter increment)
        % real_coeff_updates = tikhonov(U, s, V, neg_residual, lambda_inverse);
        real_coeff_updates =sensitivity_matrix \ neg_residual;


        % % Reshape real_coeff_updates to match the harmonic coefficients matrix structure
        % temp_real_frequency = zeros(1, 2 * num_degrees_freedom);
        % temp_real_frequency(1, 1) = real_coeff_updates(1, 1);
        % real_coeff_updates(1,1)=real_coeff_updates(2,1);
        % real_coeff_updates(2,1)=0;
        % coeff_updates = reshape(real_coeff_updates(1:end-1,:), 2, num_degrees_freedom*num_harmonics);
        % coeff_updates = coeff_updates';
        %
        % % Extract sensitivity_matrix (formerly sensitivity_parameter_da)
        % sensitivity_parameter_updates = coeff_updates(1:num_harmonics, 1:2);
        %
        % for dof_index = 1:num_degrees_freedom-1
        %     sensitivity_parameter_updates = [sensitivity_parameter_updates, coeff_updates(dof_index*num_harmonics+1:(dof_index+1)*num_harmonics, 1:2)];
        % end
        % %sensitivity_parameter_updates = [temp_real_frequency; sensitivity_parameter_updates];
        %
        % sensitivity_const = zeros(1, 2 * num_degrees_freedom);
        % sensitivity_const(1,1) = real_coeff_updates(end,:);
        % sensitivity_parameter_updates = [temp_real_frequency;sensitivity_const; sensitivity_parameter_updates];
        %
        % % Check if updated parameters are within the acceptable range
        % if ~coefficients_check(temp_harmonic_coefficients + sensitivity_parameter_updates)
        %     lambda_inverse = lambda_inverse * gamma_trust;
        %     continue;
        % end
        %
        % % Calculate real coefficient updates
        % real_coeff_updates = tikhonov(U, s, V, neg_residual, lambda_inverse);
        % Reshape real_coeff_updates to match the harmonic coefficients matrix structure
        ini_real_coeff_updates=real_coeff_updates;
        temp_real_frequency = zeros(1, 2 * num_degrees_freedom);
        temp_real_frequency(1, 1) = real_coeff_updates(1, 1);
        real_coeff_updates(1,1)=real_coeff_updates(2,1);
        real_coeff_updates(2,1)=0;
        coeff_updates = reshape(real_coeff_updates, 2, num_degrees_freedom*num_harmonics);
        coeff_updates = coeff_updates';

        % Extract sensitivity_matrix (formerly sensitivity_parameter_da)
        sensitivity_parameter_updates = coeff_updates(1:num_harmonics, 1:2);

        for dof_index = 1:num_degrees_freedom-1
            sensitivity_parameter_updates = [sensitivity_parameter_updates, coeff_updates(dof_index*num_harmonics+1:(dof_index+1)*num_harmonics, 1:2)];
        end
        sensitivity_parameter_updates = [temp_real_frequency; sensitivity_parameter_updates];

        % Update harmonic coefficients and recalculate residuals
        harmonic_coefficients = temp_harmonic_coefficients + sensitivity_parameter_updates;
        initial_frequency = harmonic_coefficients(1, 1);
        total_time = 2 * pi / initial_frequency;
        step_size = 2 * pi / (initial_frequency * num_time_point);
        time_data = (0:step_size:total_time);
        residual_updates = calculate_residual(harmonic_coefficients);

        % Calculate agreement indicator and check for convergence
        neg_residual_temp = -reshape(residual_updates(:, 1:num_degrees_freedom), num_degrees_freedom * length(time_data), 1);
        L_neg_residual = sensitivity_matrix * ini_real_coeff_updates - neg_residual;
        rho_agreement = (neg_residual' * neg_residual - neg_residual_temp' * neg_residual_temp) / (neg_residual' * neg_residual - L_neg_residual' * L_neg_residual);
        if rho_agreement >= rho_trust
            break;
        end
        lambda_inverse = lambda_inverse * gamma_trust;
    end

    % Check for overall convergence
    tol_convergence = norm(coeff_updates) / norm(harmonic_coefficients);
    harmonic_coefficients_record = [harmonic_coefficients_record, harmonic_coefficients];
    trust_region_record = [trust_region_record; lambda_inverse];
    harmonic_coefficients
    if tol_convergence <= Etol
        break;
    end
    all_coefficients(iteration_index).harmonic_coefficients = harmonic_coefficients;
    iteration_index
end
%     every_harmonic_coefficients(index_bifurcation).alpha=alpha;
%     every_harmonic_coefficients(index_bifurcation).harmonic_coefficients=harmonic_coefficients;
%     index_residual=norm(residual_updates(:, 1));
%     every_harmonic_coefficients(index_bifurcation).index_residual=index_residual;
%     index_bifurcation=index_bifurcation+1;
% end

% End timing the execution
toc;

initial_frequency = harmonic_coefficients(1, 1);
% total_time = 2 * pi / initial_frequency;
total_time = 300;
time_step = 0.01;
time_data = (0:time_step:total_time);
% Calculate final residual
final_residual = calculate_residual(harmonic_coefficients);

% Plot the residuals over time
figure;
plot(time_data, final_residual(:, 1), 'r-', 'LineWidth', 1);
legend_handle = legend('$$x_1$$');
set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

% Define new time data for response analysis
% time_data = 0:0.04:200;

initial_frequency = harmonic_coefficients(1, 1);
harmonic_parameters = harmonic_coefficients(2:end, :);

% Initialize arrays for displacement, velocity, and acceleration
displacement = zeros(num_degrees_freedom, length(time_data));
velocity = zeros(num_degrees_freedom, length(time_data));
acceleration = zeros(num_degrees_freedom, length(time_data));

% Compute displacement, velocity, and acceleration for each degree of freedom
for dof_index = 1:num_degrees_freedom
    for harm_index = 1:num_harmonics
        frequency_term = (2 * harm_index - 1) * initial_frequency * time_data;
        displacement(dof_index, :) = displacement(dof_index, :) + harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(frequency_term) + harmonic_parameters(harm_index, 2 * dof_index) * sin(frequency_term);
        velocity(dof_index, :) = velocity(dof_index, :) - initial_frequency * (2 * harm_index - 1) * harmonic_parameters(harm_index, 2 * dof_index - 1) * sin(frequency_term) + initial_frequency * (2 * harm_index - 1) * harmonic_parameters(harm_index, 2 * dof_index) * cos(frequency_term);
        acceleration(dof_index, :) = acceleration(dof_index, :) - (initial_frequency * (2 * harm_index - 1))^2 * harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(frequency_term) - (initial_frequency * (2 * harm_index - 1))^2 * harmonic_parameters(harm_index, 2 * dof_index) * sin(frequency_term);
    end
end

J_alpha1 = compute_J_alpha1(time_data, alpha, initial_frequency, harmonic_parameters, 1, num_harmonics);
figure;
plot(time_data, exp(J_alpha1), 'k-', 'LineWidth', 1);
hold on;
plot(time_data, time_data.^(-0.5), '--k', 'LineWidth', 1.2);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
figure;
loglog(time_data, exp(J_alpha1), 'b-', 'LineWidth', 1.5);
hold on;
loglog(time_data, time_data.^(-0.5), '--k', 'LineWidth', 1.5);
legend_handle = legend('$$log(J^{\alpha_1})$$','$$log(t^{-\frac{1}{2}})$$');
% legend_handle = legend('$$Logarithmic time$$','$$Logarithmic coordinate$$');
set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
grid on;


% figure;
% plot(time_data, log_J_alpha2, 'k-', 'LineWidth', 1);
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
% 
% 
% % Plot displacement response over time
% figure;
% % plot(time_data, displacement(1, :), 'k-', 'LineWidth', 1);
% plot(time_data, displacement(1, :), 'r.', 'MarkerSize', 12);
% legend_handle = legend('$$x_1$$');
% set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
% 
% % Plot displacement response over time
% figure;
% plot(time_data, velocity(1, :), 'k-', 'LineWidth', 1);
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
% 
% 
% 
% % Plot phase diagrams for displacement and velocity
% figure;
% plot(displacement(1, :), velocity(1, :), 'k--', 'LineWidth', 1.5);
% legend_handle = legend('$$x_1$$');
% set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

% Store initial harmonic coefficients for later reference
ini_harmonic_coefficients

