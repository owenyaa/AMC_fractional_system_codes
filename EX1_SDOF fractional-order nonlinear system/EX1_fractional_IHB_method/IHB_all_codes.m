%% Incremental Harmonic Balance (IHB) for Nonlinear Fractional-Order Systems
% This program is a supporting program for the first example in paper 
% "Hai-su Wang, Zhong-rong Lu, Ji-ke Liu, Guang Liu(*),
% Semi-analytical solution and nonlinear characterization analysis of fractional-order 
% nonlinear systems based on the time-domain minimum residual method",
% which solves the following equation using the IHB method
% M x''(t) + K x(t) + ε * (1 - c * x(t).^2) * D^{α}{x}(t) + k3 * x(t).^3 = 0
%% ========================================================================
%
% I. Problem Statement
% --------------------
% We consider a nonlinear system (scalar or MDOF) with fractional damping
% of van der Pol type:
%
%   M x''(t) + K x(t) + ε * (1 - c * x(t).^2) * D^{α}{x}(t) + k3 * x(t).^3 = 0
%
% where:
%   - x(t)        ∈ R^{N_dof} is the displacement vector;
%   - M, K        are mass and stiffness (scalars or N_dof×N_dof matrices);
%   - ε < 0       is the van der Pol–type coefficient (named vdp_coefficient);
%   - c > 0       is the parameter coefficient (parameter_coefficient);
%   - k3          is the cubic stiffness (cubic_stiffness_matrix; scalar/matrix; scalar here);
%   - D^{α}       is the fractional derivative of order α∈(0,2] (α∈(0,1] in this example).
%
% Goal: compute the periodic response x(t) and the unknown angular frequency ω
% such that the equation holds in the time domain.
%
% Note: although the example is SDOF for demonstration, all data structures,
% packing/unpacking conventions, IHB linearization, and spectral operations
% are designed for MDOF and can be directly extended to N_dof>1 without changing
% the core interfaces.
%
%
% II. Method Overview
% -------------------
% 1) Spectral approximation: represent each DOF with a truncated Fourier series
%       x_d(t) ≈ Σ_{k=0}^{H-1} [ A_{k,d} cos(k ω t) + B_{k,d} sin(k ω t) ].
%    Stack all DOFs into a “series coefficient matrix”:
%       series ∈ R^{H × (2*N_dof)}, each DOF occupies two columns:
%       [cos-series, sin-series]. Row 1 is the DC term (k=0); rows 2..H are k=1..H-1.
%
% 2) IHB linearized iteration: perform Newton-like corrections on {coefficients, ω}.
%    At the n-th iteration, assemble residual R(a, ω) and Jacobians:
%       Ja = ∂R/∂a,   Jw = ∂R/∂ω,
%    solve increments [Δa; Δω], then update (a, ω) ← (a, ω) + (Δa, Δω).
%
% 3) Fractional derivative in coefficient space:
%    For each harmonic k≥1, D^{α} acting on cos/sin follows the smooth spectral rule:
%       c1 = k^α cos(απ/2),  c2 = k^α sin(απ/2)
%       D^{α}{cos} →  c1 * cos  +  c2 * sin
%       D^{α}{sin} →  c1 * sin  −  c2 * cos
%    Hence D^{α}x is obtained by a per-DOF linear transform on the series matrix columns.
%
% 4) Nonlinear products via “series-space convolution”:
%    Time-domain products such as x^2, x^3, x^2 D^{α}x correspond to convolutions
%    in the spectral domain. We implement two-series convolution by
%    harmonic_product_convolution(·,·), and chain multiple convolutions by
%    multiply_series_coeffs(·).
%
%
% III. Data Layout & Naming
% -------------------------
% * Unified series matrix convention: H × (2*N_dof)
%     [cos_1  sin_1 | cos_2  sin_2 | ... | cos_{N_dof} sin_{N_dof}],
%     row 1 = DC (k=0); rows 2..H correspond to k=1..H-1.
%
% * Packing / Unpacking conventions (compatible with legacy code):
%     - pack_series_coeffs_per_dof(series)  →  packed ∈ R^{(2H-1) × N_dof}
%           Each DOF is one column: cosine block (incl. DC) stacked on top of
%           sine block (starting from k=1, no DC).
%     - stack_columns_vector(packed)        →  vec ∈ R^{(2H-1)*N_dof}
%     - unstack_columns_by_dof(vec)         →  packed (inverse)
%     - columns_to_series_matrix(packed)    →  series (inverse)
%
% * Globals use descriptive names: number_of_harmonics, mass_matrix,
%   stiffness_matrix, vdp_coefficient, parameter_coefficient, cubic_stiffness_matrix,
%   fractional_order, num_dofs, etc., and are reusable across files.
%
%
% IV. Execution Workflow
% ----------------------
% Step 0. Set globals and initial guesses (harmonic truncation H, α, M/K/ε/k3, etc.;
%         initial a and ω).
%
% Step 1. IHB main loop (see the while block):
%   1.1 Compute ihb_blocks = compute_ihb_residual_and_jacobians(a, ω)
%       - ihb_blocks(1).vector = R (already packed by DOF and stacked)
%       - ihb_blocks(2).vector = Ja (square Jacobian)
%       - ihb_blocks(3).vector = Jw (frequency Jacobian, vector)
%   1.2 Augment Ja by appending Jw as the last column; solve Δall = [Δa; Δω].
%   1.3 Extract Δω; convert Δa from “stacked vector” → “packed by DOF”
%       → “series layout”.
%   1.4 Update unknowns: a ← a + Δa, ω ← ω + Δω; compute ‖R‖ and check
%       convergence/divergence guard.
%
% Step 2. Time-Domain Reconstruction
%   Using the converged (a, ω), synthesize x(t), x'(t), x''(t), and D^{α}x(t)
%   from the series coefficients; evaluate the time history of the equation LHS
%   as a residual check; plot phase portraits, etc.
%
%
% V. Key Building Blocks
% ----------------------
% 1) compute_ihb_residual_and_jacobians(a, ω)
%    - Use derive_series_coeffs to get x'' in series form;
%      use fractional_derive_series_coeffs to get D^{α}x in series form.
%    - Residual series:
%         R_series = -[ M*(ω^2 x'') + K*x + ε*(ω^α D^{α}x) ] - R_nl(a, ω),
%      where R_nl is the first component from build_nonlinear_terms.
%    - Frequency Jacobian series:
%         d/dω[ M*(ω^2 x'') ] = 2ω M x'',
%         d/dω[ ε*(ω^α D^{α}x) ]  →  ε * D^{α}x  (implemented per the original convention),
%      plus the nonlinear ω-term (third component from build_nonlinear_terms).
%    - Convert series to packed and stacked formats and return to the main loop.
%
% 2) assemble_coefficients_jacobian(a, ω)
%    - For each unit increment basis δa(i,j) (set one entry to 1 in series layout),
%      construct the Ja column:
%         Ja_col =  M*(ω^2 δa'') + ε*(ω^α δ(D^{α}a)) + K*δa + Ja_nl(a, δa, ω),
%      where Ja_nl is the second component from build_nonlinear_terms.
%    - After packing and stacking, assemble the full square Jacobian (identical
%      to the legacy “jacobi”).
%
% 3) build_nonlinear_terms(ω, a, Δa)
%    - For the vdp–Duffing nonlinearity, build three components per DOF
%      (all in series layout):
%        (1) Residual nonlinearity:  k3*x^3  −  ε*c*ω^α * x^2 * D^{α}x
%        (2) Ja nonlinearity:       3*k3*x^2*Δx  −  ε*c*ω^α * x^2 * D^{α}Δx
%                                   − 2*ε*c*ω^α * x * D^{α}x * Δx
%        (3) Jw nonlinearity:       (−ε) * ( x * x * D^{α}x )
%      All time-domain products are implemented as spectral convolutions via
%      harmonic_product_convolution and multiply_series_coeffs; lastly each DOF’s
%      result is written back into the global H × (2*N_dof) container.
%
% 4) derive_series_coeffs(series, n)
%    - Implement integer-order time derivatives in coefficient space: odd/even
%      orders swap/keep cos/sin with appropriate signs; DC row (k=0) is zero
%      for n≥1; sine’s DC is kept at zero.
%
% 5) fractional_derive_series_coeffs(series, α)
%    - Implement fractional derivative D^{α} in coefficient space using
%      the spectral rule (c1, c2).
%
% 6) harmonic_product_convolution(A, B) / multiply_series_coeffs(list)
%    - Implement spectral convolution for time-domain products: four blocks
%      (cos·cos, cos·sin, sin·cos, sin·sin) accumulate sum/difference frequencies,
%      and respect sin(0)=0.
%
% 7) apply_linear_operator_to_series(L, series)
%    - Apply scalar or square matrix L to a series matrix by per-DOF 2-column
%      linear combinations; suitable for M, K, or scaled ε.
%
% 8) sum_series_coeffs(list) / add_series_coeff_mats(A, B)
%    - Row-aligned series summation with automatic row padding; enforce
%      DC sine (row 1, even columns) = 0.
%
% 9) pack_series_coeffs_per_dof / columns_to_series_matrix /
%    stack_columns_vector / unstack_columns_by_dof
%    - Define reversible transforms between “series matrix ↔ packed columns ↔
%      stacked vector”, which are essential for assembling IHB linear algebra
%      (kept compatible with legacy code).
%
%
% VI. Convergence & Numerics
% --------------------------
% * Residual measure: the main loop monitors ||R||_2; adjust tolerance to your scale.
% * Divergence guard: if iterations > 500 and ||R||_2 > 1, a simple divergence
%   detection is triggered to avoid wasting iterations.
% * Frequency derivative convention: in Jw, d/dω[ε*(ω^α D^{α}x)] is implemented as
%   “ε * D^{α}x”, consistent with the original pipeline (if you need a stricter
%   α-dependent form, modify compute_ihb_residual_and_jacobians accordingly).
% * Harmonic truncation: too small H causes aliasing; large H increases convolution
%   cost roughly O(H^2).
%
%
% VII. Extending to New Nonlinear Models
% --------------------------------------
% * Only modify build_nonlinear_terms:
%   - Residual component: express your nonlinear term as time-domain products
%     and convert them via convolution to series layout.
%   - Ja component: first-order linearization of the nonlinear term w.r.t. Δa.
%   - Jw component: first-order derivative of the nonlinear term w.r.t. ω
%     (if explicit ω-dependence exists).
% * Other modules (derivatives, packing, linear operators, summations) are reusable
%   without changes.
%
%
% VIII. Troubleshooting
% ---------------------
% * Dimension mismatch (most common): convolution yields H_out = H_A + H_B − 1,
%   which **must be trimmed** to H rows when writing back to global arrays.
%   This code consistently applies `x = x(1:H,:)` after every convolution.
% * Mixed old/new names: if you see “function/variable not found”, check for
%   legacy calls such as:
%     - unstack_by_dof        → should be unstack_columns_by_dof
%     - vec_stack_columns     → should be stack_columns_vector
%   You may add small alias shims at the end of the file for backward compatibility.
% * Validity of α: warnings are issued for α∉[0,2]. Please correct α first.
%
%
% IX. Symbols & Globals
% ---------------------
%   number_of_harmonics(H)     number of harmonic rows (including DC);
%                              row 1=DC, rows 2..H correspond to k=1..H-1
%   num_dofs                   number of DOFs
%   mass_matrix (M)            scalar or matrix
%   stiffness_matrix (K)       scalar or matrix
%   vdp_coefficient (ε)        van der Pol coefficient (negative here)
%   parameter_coefficient (c)  parameter inside (1 - c x^2)
%   cubic_stiffness_matrix(k3) cubic stiffness (scalar here)
%   fractional_order (α)       fractional derivative order
%
%
% X. Post-Processing
% ------------------
% * With converged (a, ω), synthesize x(t), x'(t), x''(t), D^{α}x(t) and plot;
% * Compute the equation LHS(t) as a time-domain residual check (closer to 0 is better);
% * Plot phase portraits (x vs x') to observe limit cycles/periodic responses.
%
% The above conventions ensure smooth reuse from SDOF to MDOF and facilitate
% replacing the nonlinear model. If you need a stricter α-dependent form in Jw
% or to introduce external excitation/coupled nonlinearities, only local edits
% in build_nonlinear_terms and compute_ihb_residual_and_jacobians are required.
%
%% ========================================================================
%
% 一、要解决的数学问题（Problem Statement）
% ------------------------------------------
% 我们考虑如下带有分数阶阻尼（van der Pol 型）的非线性系统（标量或 MDOF）：
%
%   M x''(t) + K x(t) + ε * (1 - c * x(t).^2) * D^{α}{x}(t) + k3 * x(t).^3 = 0
%
% 其中：
%   - x(t)        ∈ R^{N_dof} 为位移向量；
%   - M, K        分别为质量与刚度（可为标量或 N_dof×N_dof 矩阵）；
%   - ε < 0       为 van der Pol 型系数（此处命名 vdp_coefficient）；
%   - c > 0       为参数系数（parameter_coefficient）；
%   - k3          三次刚度（cubic_stiffness_matrix；可为标量/矩阵，本算例为标量）；
%   - D^{α}       为阶次 α∈(0,2] 的分数阶导数（本算例取 α∈(0,1] 常用情形）。
%
% 目标：求解系统的**周期性响应** x(t) 及其**未知角频率** ω，使得上式在时间域成立。
%
% 注：代码以 SDOF 演示为主，但所有数据结构、打包/解包约定、IHB 线性化与谱运算
%     都面向 MDOF 设计，直接可扩展到 N_dof>1 的模型（不需要改核心接口）。
%
%
% 二、核心思想（Method Overview）
% ------------------------------
% 1) 频域近似：用截断的 Fourier 级数表示每个自由度的位移：
%       x_d(t) ≈ Σ_{k=0}^{H-1} [ A_{k,d} cos(k ω t) + B_{k,d} sin(k ω t) ]
%    将所有 DOF 的系数装成一个“级数系数矩阵”：
%       series ∈ R^{H × (2*N_dof)}, 每个 DOF 占两列：[cos 系列列, sin 系列列]，
%       第 1 行是 DC（k=0）；第 2..H 行对应谐波 k=1..H-1。
%
% 2) IHB 线性化迭代：对“系数 + 频率”做牛顿式校正。
%    在第 n 次迭代构造残差 R(a, ω) 及其关于 a 与 ω 的雅可比：
%       Ja = ∂R/∂a,   Jw = ∂R/∂ω
%    然后求解增量 [Δa; Δω]，并更新 (a, ω) ← (a, ω) + (Δa, Δω)。
%
% 3) 分数阶导数在“系数空间”的实现：
%    对于每一阶谐波 k≥1，D^{α} 作用到 cos/sin 的系数变换遵循光滑的谱规则：
%       c1 = k^α cos(απ/2),  c2 = k^α sin(απ/2)
%       D^{α}{cos} → c1 * cos  +  c2 * sin
%       D^{α}{sin} → c1 * sin  −  c2 * cos
%    故可直接在“级数系数矩阵”上做列内线性变换，得到 D^{α}x 的系数矩阵。
%
% 4) 非线性项的“系数空间卷积”：
%    x^2、x^3、x^2 D^{α}x 等时间域乘积，对应频域的卷积。我们用
%    harmonic_product_convolution(·,·) 实现两列级数的三角恒等式卷积，
%    multiply_series_coeffs(·) 顺序串接多重卷积（如 x·x·x）。
%
%
% 三、数据布局与统一命名（Data Layout & Naming）
% ----------------------------------------------
% * series 矩阵统一约定：H × (2*N_dof)
%     [cos_1  sin_1 | cos_2  sin_2 | ... | cos_{N_dof} sin_{N_dof}]
%     行 1 = DC（k=0）；行 2..H = k=1..H-1。
%
% * “打包 / 解包”约定（与原始旧代码兼容）：
%     - pack_series_coeffs_per_dof(series)  →  packed ∈ R^{(2H-1) × N_dof}
%           每个 DOF 一列：先堆 cos（含 DC），再堆 sin（从 k=1 起，无 DC）
%     - stack_columns_vector(packed)        →  vec ∈ R^{(2H-1)*N_dof}
%     - unstack_columns_by_dof(vec)         →  packed（逆操作）
%     - columns_to_series_matrix(packed)    →  series（逆操作）
%
% * 全局（global）命名：number_of_harmonics, mass_matrix, stiffness_matrix,
%   vdp_coefficient, parameter_coefficient, cubic_stiffness_matrix,
%   fractional_order, num_dofs 等，均为**可跨文件复用**的统一名。
%
%
% 四、主流程（Execution Workflow）
% -------------------------------
% Step 0. 设定全局参数与初始猜测（harmonics 截断阶 H、α、M/K/ε/k3 等；a 与 ω 初值）。
%
% Step 1. IHB 迭代主循环（见 while 块）：
%   1.1 计算：ihb_blocks = compute_ihb_residual_and_jacobians(a, ω)
%       - ihb_blocks(1).vector = R（已按 DOF 打包并堆成列向量）
%       - ihb_blocks(2).vector = Ja（大雅可比，方阵）
%       - ihb_blocks(3).vector = Jw（关于 ω 的雅可比，向量）
%   1.2 组装增广线性系统：将 Jw 写入 Ja 的最后一列，解出 Δall = [Δa; Δω]
%   1.3 解增量：拆出 Δω，并将 Δa 从“列向量”→“packed by DOF”→“series 布局”
%   1.4 更新未知量：a ← a + Δa， ω ← ω + Δω；计算残差范数并检查收敛/发散保护。
%
% Step 2. 时域重构（Time-Domain Reconstruction）
%   使用收敛后的 (a, ω)，通过谐波合成得到 x(t)、x'(t)、x''(t) 与 D^{α}x(t)。
%   进一步可评估“方程左端 LHS”的时间历程 residual(t)，并绘制相图等。
%
%
% 五、IHB 关键构件（Key Building Blocks）
% -------------------------------------
% 1) compute_ihb_residual_and_jacobians(a, ω)
%    - 先用 derive_series_coeffs 求 x'' 的系数矩阵；
%      用 fractional_derive_series_coeffs 求 D^{α}x 的系数矩阵；
%    - 残差 R 的系数矩阵：
%         R_series = -[ M*(ω^2 x'') + K*x + ε*(ω^α D^{α}x) ] - R_nl(a, ω)
%      其中 R_nl 来自 build_nonlinear_terms 的第 1 个分量；
%    - 频率雅可比 Jw 的系数矩阵：
%         d/dω[ M*(ω^2 x'') ] = 2ω M x''，
%         d/dω[ ε*(ω^α D^{α}x) ] → ε * D^{α}x  （按原工程约定实现）
%      再加上非线性对 ω 的项（build_nonlinear_terms 第 3 个分量）；
%    - 系数矩阵经 pack_series_coeffs_per_dof 与 stack_columns_vector
%      变为向量/方阵，返回给主循环。
%
% 2) assemble_coefficients_jacobian(a, ω)
%    - 对每个“单元增量基” δa(i,j)（在 series 布局中把某个位置置 1，其余为 0），
%      构造 Ja 的列：
%         Ja_col =  M*(ω^2 δa'') + ε*(ω^α δ(D^{α}a)) + K*δa + Ja_nl(a, δa, ω)
%      其中 Ja_nl 由 build_nonlinear_terms 的第 2 个分量提供；
%    - 列向量打包与堆叠后，组装为大雅可比（与原“jacobi”完全等价）。
%
% 3) build_nonlinear_terms(ω, a, Δa)
%    - 针对 vdp–Duffing 型非线性，逐 DOF 构造三个分量（均为 series 布局）：
%        (1) Residual 非线性：  k3*x^3  −  ε*c*ω^α * x^2 * D^{α}x
%        (2) Ja 非线性：       3*k3*x^2*Δx  −  ε*c*ω^α * x^2 * D^{α}Δx
%                               − 2*ε*c*ω^α * x * D^{α}x * Δx
%        (3) Jw 非线性：       (−ε) * ( x * x * D^{α}x )
%      其中时间域乘法均用 harmonic_product_convolution → multiply_series_coeffs
%      在系数空间做卷积；最后逐 DOF 写回全局 H × (2*N_dof) 容器。
%
% 4) derive_series_coeffs(series, n)
%    - 在系数空间实现**整阶**时间导数 n：奇偶阶分别做 cos/sin 的交叉和符号；
%      DC 行（k=0）对 n≥1 的导数为 0；保持 sin 的 DC 为 0。
%
% 5) fractional_derive_series_coeffs(series, α)
%    - 在系数空间实现**分数阶**导数 D^{α}，用上面的谱规则（c1, c2）。
%
% 6) harmonic_product_convolution(A, B) / multiply_series_coeffs(list)
%    - 实现时间域乘积在谐波系数上的卷积：分四类（cos·cos、cos·sin、sin·cos、sin·sin）
%      分块累加，自动处理“和/差频”，并遵循 sin(0)=0 的约束。
%
% 7) apply_linear_operator_to_series(L, series)
%    - 将标量或方阵 L 作用于 series（逐 DOF 的 2 列块进行线性组合），
%      适用于 M、K 或缩放后的 ε 等。
%
% 8) sum_series_coeffs(list) / add_series_coeff_mats(A, B)
%    - series 逐行对齐求和：自动做行数补齐；强制输出的“DC 行的 sin 列”为 0。
%
% 9) pack_series_coeffs_per_dof / columns_to_series_matrix /
%    stack_columns_vector / unstack_columns_by_dof
%    - 这组函数定义了“系数矩阵 ↔ 打包列 ↔ 列向量”的互逆变换，
%      是 IHB 线性代数拼装的关键约定（与旧工程保持兼容）。
%
%
% 六、收敛与数值细节（Convergence & Numerics）
% ---------------------------------------------
% * 残差度量：主循环用 ||R||_2 监控；你可根据问题尺度调整 tolerance。
% * 发散保护：若迭代数 > 500 且 ||R||_2 > 1，触发简单“发散检测”以避免无意义迭代。
% * 频率导数约定：Jw 中 d/dω[ε*(ω^α D^{α}x)] 采用“ε * D^{α}x”的实现，
%   与原始代码链保持一致（如需严格 α 依赖项，可在 compute_ihb_residual_and_jacobians 定制）。
% * 谐波截断：H 过小会带来混叠误差；H 较大时卷积成本随之增加（O(H^2) 级）。
%
%
% 七、如何替换/拓展非线性（Extending to New Models）
% ---------------------------------------------------
% * 仅需在 build_nonlinear_terms 中替换非线性构造：
%   - Residual 分量：把目标模型的非线性项写成“时间域乘积”，再用卷积得到 series。
%   - Ja 分量：对“非线性项”对系数 Δa 的一阶线性化（可按物理推导写出）。
%   - Jw 分量：对“非线性项”对 ω 的一阶导数（若存在显式 ω 依赖）。
% * 其他模块（导数、打包、线性算子、汇总）均可复用，无需改动。
%
%
% 八、常见排错（Troubleshooting）
% -------------------------------
% * 维度不匹配（最常见）：卷积产生 H_out = H_A + H_B − 1 的更长 series，**必须裁剪**
%   到 H 行再写回全局容器；本工程已在所有卷积后统一做 `x = x(1:H,:)` 的裁剪。
% * 旧名/新名混用：若报“找不到函数/变量”，核对是否还有旧名调用：
%     - unstack_by_dof → 应为 unstack_columns_by_dof
%     - vec_stack_columns → 应为 stack_columns_vector
%   可以在工具尾部加“别名垫片”以向后兼容。
% * 分数阶 α 的合法性：函数会对 α∈[0,2] 做提示（warning），越界请先纠正。
%
%
% 九、符号与全局（Symbols & Globals）
% -----------------------------------
%   number_of_harmonics(H)     谐波行数（含 DC）；行 1=DC，行 2..H=k=1..H-1
%   num_dofs                   自由度数
%   mass_matrix (M)            标量或矩阵
%   stiffness_matrix (K)       标量或矩阵
%   vdp_coefficient (ε)        van der Pol 系数（本算例为负）
%   parameter_coefficient (c)  非线性参数（在 1 - c x^2 内部）
%   cubic_stiffness_matrix(k3) 三次刚度（标量演示）
%   fractional_order (α)       分数阶阶次
%
%
% 十、可视化与后处理（Post-Processing）
% ------------------------------------
% * 利用收敛的 (a, ω) 合成 x(t)、x'(t)、x''(t)、D^{α}x(t) 并作图；
% * 计算方程左端“LHS(t)”作为时域残差检验（越接近 0 越好）；
% * 相图（x 对 x'）用于观察极限环/周期响应的形态。
%
% 以上约定确保本程序可在 SDOF 到 MDOF 间平滑复用，并能方便替换非线性模型。
% 若你需将 Jw 的 α-依赖推导为更严格的形式或引入外激励/耦合非线性，只需在
% build_nonlinear_terms 与 compute_ihb_residual_and_jacobians 内部做局部修改。
%
%% ========================================================================

clear; clc; close all;
tic

%% ----------------------------- Global settings -----------------------------
global number_of_harmonics mass_matrix vdp_coefficient stiffness_matrix ...
       num_dofs cubic_stiffness_matrix fractional_order parameter_coefficient

number_of_harmonics    = 30;
parameter_coefficient  = 1.1;   % (1 - c*x^2)
vdp_coefficient        = -0.8;  % van der Pol–type coefficient
natural_frequency_sq   = 1;     % w0^2
cubic_stiffness_matrix = 1;     % cubic stiffness
fractional_order       = 0.5;   % α in (0,2], typically (0,1]
mass_matrix            = 1;     % scalar or matrix
stiffness_matrix       = natural_frequency_sq;
num_dofs               = size(mass_matrix, 1);

%% ---------------------- Initial guess (coeffs & frequency) -----------------
angular_frequency  = 1.855;
coefficient_matrix = zeros(number_of_harmonics, 2*num_dofs);
coefficient_matrix(2, :) = [1.3, -1.3];   % initial 1st-harmonic cos/sin

%% ---------------------- IHB iterative correction loop ----------------------
residual_norm   = 1;
tolerance       = 1e-10;
delta_frequency = 0;
iteration       = 0;

while (residual_norm >= tolerance && iteration < 1001)
    % Assemble residual R, Jacobian wrt coefficients (Ja), and wrt frequency (Jw)
    ihb_blocks = compute_ihb_residual_and_jacobians(coefficient_matrix, angular_frequency);
    residual_vector       = ihb_blocks(1).vector;   % R
    coefficients_jacobian = ihb_blocks(2).vector;   % Ja = dR/da (stacked layout)
    frequency_jacobian    = ihb_blocks(3).vector;   % Jw = dR/dω (stacked layout)

    % Augment Ja with Jw as the last column and solve for increments
    coefficients_jacobian(:, number_of_harmonics + 1) = frequency_jacobian;
    delta_all = coefficients_jacobian \ residual_vector;

    % Extract frequency increment; keep only coefficient increments in delta_all
    delta_frequency = delta_all(number_of_harmonics + 1);
    delta_all(number_of_harmonics + 1) = 0;

    % Reshape stacked increments -> [H x (2*num_dofs)] series matrix
    delta_all = unstack_columns_by_dof(delta_all);
    delta_all = columns_to_series_matrix(delta_all);

    % Update unknowns
    residual_norm      = norm(residual_vector);
    angular_frequency  = angular_frequency + delta_frequency;
    coefficient_matrix = coefficient_matrix + delta_all;

    iteration = iteration + 1;

    % Simple divergence guard
    if residual_norm > 1 && iteration > 500
        disp('Divergence detected: residual > 1 after 500 iterations.');
        break
    end
end

toc

%% -------------------------- Time-domain reconstruction ---------------------
time_step   = 0.01;
time_vector = 0:time_step:20;

displacement                       = zeros(num_dofs, length(time_vector));
velocity                           = zeros(num_dofs, length(time_vector));
acceleration                       = zeros(num_dofs, length(time_vector));
fractional_derivative_displacement = zeros(num_dofs, length(time_vector));

% Exclude DC row; harmonics 1..(H-1)
harmonic_coeffs = coefficient_matrix(2:end, :);

% Fractional derivative in coefficient space (same layout as coefficient_matrix)
coeffs_fractional_derivative = fractional_derive_series_coeffs(coefficient_matrix, fractional_order);

% Build combined series coefficients for D^{α}(x) and related terms
fractional_operator_terms(1).matrix = ...
    -vdp_coefficient * parameter_coefficient * angular_frequency^fractional_order * coefficient_matrix;
fractional_operator_terms(2).matrix = coefficient_matrix;
fractional_operator_terms(3).matrix = coeffs_fractional_derivative;

combined_fractional_series(1).matrix = multiply_series_coeffs(fractional_operator_terms);
combined_fractional_series(2).matrix = apply_linear_operator_to_series( ...
                                         vdp_coefficient, ...
                                         angular_frequency^fractional_order * coeffs_fractional_derivative);
combined_fractional_series = sum_series_coeffs(combined_fractional_series);

% Remove DC row to align with harmonic indices
combined_fractional_series = combined_fractional_series(2:end, :);

% Synthesize time histories from series coefficients
for dof = 1:num_dofs
    for h = 1:(number_of_harmonics - 1)
        cos_term = cos(h * angular_frequency * time_vector);
        sin_term = sin(h * angular_frequency * time_vector);

        % x(t)
        displacement(dof, :) = displacement(dof, :) ...
            + harmonic_coeffs(h, 2*dof - 1) * cos_term ...
            + harmonic_coeffs(h, 2*dof)     * sin_term;

        % D^{α}x(t)
        fractional_derivative_displacement(dof, :) = fractional_derivative_displacement(dof, :) ...
            + combined_fractional_series(h, 2*dof - 1) * cos_term ...
            + combined_fractional_series(h, 2*dof)     * sin_term;

        % x'(t)
        velocity(dof, :) = velocity(dof, :) ...
            - angular_frequency * h * harmonic_coeffs(h, 2*dof - 1) * sin_term ...
            + angular_frequency * h * harmonic_coeffs(h, 2*dof)     * cos_term;

        % x''(t)
        acceleration(dof, :) = acceleration(dof, :) ...
            - (angular_frequency * h)^2 * harmonic_coeffs(h, 2*dof - 1) * cos_term ...
            - (angular_frequency * h)^2 * harmonic_coeffs(h, 2*dof)     * sin_term;
    end

    % Add DC component to displacement
    displacement(dof, :) = displacement(dof, :) + coefficient_matrix(2*dof - 1, 1);
end

%% -------------------------- Governing-equation residual --------------------
% For each DOF, the equation (LHS) reads approximately:
%   mass_matrix*x'' + stiffness_matrix*x
% + vdp_coefficient*(1 - parameter_coefficient*x.^2)*D^{fractional_order}(x)
% + cubic_stiffness_matrix*x.^3  ≈ 0
equation_residual = mass_matrix * acceleration ...
                  + fractional_derivative_displacement ...
                  + stiffness_matrix * displacement ...
                  + cubic_stiffness_matrix * displacement.^3;

%% --------------------------------- Plots -----------------------------------
figure;
plot(time_vector, equation_residual, 'k-', 'LineWidth', 1.5);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
xlabel('Time'); ylabel('Residual (equation LHS)');

sample_stride = 20;

figure;
plot(time_vector(10:sample_stride:end), displacement(1, 10:sample_stride:end), ...
     'bo', 'MarkerSize', 5, 'LineWidth', 1.5);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
xlabel('Time'); ylabel('Displacement');

figure;
plot(displacement(1, 8:sample_stride:end), velocity(1, 8:sample_stride:end), 'ro', 'MarkerSize', 4);
legend_handle = legend('$$TMRM$$','$$Newmark-\beta$$','$$IHB\ method$$');
set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
xlabel('Displacement'); ylabel('Velocity');


function big_jacobian = assemble_coefficients_jacobian(coefficient_matrix, angular_frequency)
%% Assemble Jacobian dR/da for self-excited periodic solution (IHB)
% Contract matches original `jacobi`:
%   Input:
%     coefficient_matrix : [H x (2*num_dofs)] current Fourier coeffs
%     angular_frequency  : ω
%   Output:
%     big_jacobian       : [(2H-1)*num_dofs  x  (2H-1)*num_dofs]
%
% Method
%   For each unit increment basis δa(i,j):
%     - Build derivatives δa_2 and δa_alpha in series space
%     - Linear operator pieces: M*ω^2*δa_2,  ε*ω^α*δa_alpha,  K*δa
%     - Nonlinear Ja contribution via build_nonlinear_terms(…, δa)
%     - Sum → trim to first H rows → pack per DOF → stack columns → place into S
%   Finally, rearrange S columns to form the square Jacobian (same as original).

global number_of_harmonics mass_matrix vdp_coefficient stiffness_matrix ...
       num_dofs fractional_order

H = number_of_harmonics;

% S will store column-wise assembled blocks before final rearrangement
S = zeros( (2*H-1)*num_dofs, 2*H*num_dofs );

for dof_col = 1:(2*num_dofs)            % j-index over columns of a
    for row_h = 1:H                      % i-index over harmonic rows
        % Unit increment basis δa
        delta_coeffs = zeros(H, 2*num_dofs);
        delta_coeffs(row_h, dof_col) = 1;

        % Derivatives in series space
        delta_ddot   = derive_series_coeffs(delta_coeffs, 2);
        delta_frac   = fractional_derive_series_coeffs(delta_coeffs, fractional_order);

        % Linear operator contributions: M*ω^2*δa_2 + ε*ω^α*δa_α + K*δa
        Ja_terms(1).matrix = apply_linear_operator_to_series(mass_matrix,      angular_frequency^2       * delta_ddot);
        Ja_terms(2).matrix = apply_linear_operator_to_series(vdp_coefficient,  angular_frequency^fractional_order * delta_frac);
        Ja_terms(3).matrix = apply_linear_operator_to_series(stiffness_matrix,                                 delta_coeffs);

        % Nonlinear contribution (Ja part) via build_nonlinear_terms
        nl = build_nonlinear_terms(angular_frequency, coefficient_matrix, delta_coeffs);
        Ja_terms(4).matrix = nl(2).part;

        % Sum, trim to first H rows, pack per DOF, then stack into a vector column
        summed = sum_series_coeffs(Ja_terms);
        summed = summed(1:H, :);
        packed = pack_series_coeffs_per_dof(summed);
        colvec = stack_columns_vector(packed);

        S(:, row_h + (dof_col-1)*H) = colvec;
    end
end

% Rearrangement to final square Jacobian (same as original indexing scheme)
big_jacobian = zeros( (2*H-1)*num_dofs, (2*H-1)*num_dofs );
for dof = 1:num_dofs
    % cos columns block
    big_jacobian(:, (dof-1)*(2*H-1)+1 : (dof-1)*(2*H-1)+H) = ...
        S(:, (dof-1)*(2*H) + 1 : (dof-1)*(2*H) + H);
    % sin columns block (skip the DC sine column inside S)
    big_jacobian(:, (dof-1)*(2*H-1) + H + 1 : (dof-1)*(2*H-1) + 2*H - 1) = ...
        S(:, (dof-1)*(2*H) + H + 2 : (dof-1)*(2*H) + 2*H);
end
end

function nonlinear_terms = build_nonlinear_terms(angular_frequency, coefficient_matrix, increment_matrix)
%% Build nonlinear contributions for IHB (Residual / Ja / Jw) in MDOF form
% Outputs (contract kept compatible with original `nonterm`):
%   nonlinear_terms(1).part : nonlinear contribution to Residual   (R)
%   nonlinear_terms(2).part : nonlinear contribution to dR/d(a)    (Ja)
%   nonlinear_terms(3).part : nonlinear contribution to dR/d(ω)    (Jw)
%
% Model (van der Pol–Duffing type, per-DOF uncoupled nonlinearity):
%   x'' + vdp*(1 - param*x^2)*D^{alpha}(x) + k*x + k3*x^3 = 0
% Nonlinear pieces assembled here (consistent with previous pipeline):
%   R_nl  =  k3*x^3  -  vdp*param*ω^α * x^2 * D^{α}x
%   Ja_nl =  3*k3*x^2*Δx  -  vdp*param*ω^α * x^2 * D^{α}Δx  -  2*vdp*param*ω^α * x * D^{α}x * Δx
%   Jw_nl =  (-vdp) * ( x * x * D^{α}x )

global num_dofs cubic_stiffness_matrix fractional_order ...
       vdp_coefficient parameter_coefficient number_of_harmonics

H = number_of_harmonics;

% Split per DOF
for dof_index = 1:num_dofs
    disp_series_by_dof(dof_index).series = coefficient_matrix(:, 2*dof_index-1 : 2*dof_index); %#ok<AGROW>
    incr_series_by_dof(dof_index).series = increment_matrix(:,   2*dof_index-1 : 2*dof_index); %#ok<AGROW>
end

% Precompute fractional derivatives (series space)
for dof_index = 1:num_dofs
    frac_deriv_disp_by_dof(dof_index).series = ...
        fractional_derive_series_coeffs(disp_series_by_dof(dof_index).series, fractional_order); %#ok<AGROW>
    frac_deriv_incr_by_dof(dof_index).series = ...
        fractional_derive_series_coeffs(incr_series_by_dof(dof_index).series, fractional_order); %#ok<AGROW>
end

% Allocate accumulators (final size must be H x (2*num_dofs))
residual_nonlinear_series            = zeros(H, 2*num_dofs);
jacobian_wrt_coeffs_nonlinear_series = zeros(H, 2*num_dofs);
jacobian_wrt_frequency_nonlinear     = zeros(H, 2*num_dofs);

%% ------------------------ (1) Residual nonlinearity ------------------------
for dof_index = 1:num_dofs
    x_series   = disp_series_by_dof(dof_index).series;
    Dx_series  = frac_deriv_disp_by_dof(dof_index).series;

    % k3 * x^3
    blk1(1).matrix = cubic_stiffness_matrix * x_series;
    blk1(2).matrix = x_series;
    blk1(3).matrix = x_series;
    x_cubed = multiply_series_coeffs(blk1);
    x_cubed = x_cubed(1:H, :);  % trim

    % - vdp*param*ω^α * x^2 * D^{α}x
    blk2(1).matrix = -vdp_coefficient * parameter_coefficient * angular_frequency^fractional_order * x_series;
    blk2(2).matrix =  x_series;
    blk2(3).matrix =  Dx_series;
    x2_Dx = multiply_series_coeffs(blk2);
    x2_Dx = x2_Dx(1:H, :);      % trim

    sum_list(1).matrix = x_cubed;
    sum_list(2).matrix = x2_Dx;
    residual_for_dof   = sum_series_coeffs(sum_list);
    residual_for_dof   = residual_for_dof(1:H, :); % safety trim

    residual_nonlinear_series(:, 2*dof_index-1 : 2*dof_index) = residual_for_dof;
end

%% ------------------ (2) Ja nonlinearity (w.r.t. coefficients) --------------
for dof_index = 1:num_dofs
    x_series    = disp_series_by_dof(dof_index).series;
    Dx_series   = frac_deriv_disp_by_dof(dof_index).series;
    da_series   = incr_series_by_dof(dof_index).series;
    Dda_series  = frac_deriv_incr_by_dof(dof_index).series;

    % 3*k3*x^2*Δx
    jb1(1).matrix = 3 * cubic_stiffness_matrix * x_series;
    jb1(2).matrix =     x_series;
    jb1(3).matrix =     da_series;
    t1 = multiply_series_coeffs(jb1); t1 = t1(1:H, :);

    % - vdp*param*ω^α * x^2 * D^{α}Δx
    jb2(1).matrix = -vdp_coefficient * parameter_coefficient * angular_frequency^fractional_order * x_series;
    jb2(2).matrix =  x_series;
    jb2(3).matrix =  Dda_series;
    t2 = multiply_series_coeffs(jb2); t2 = t2(1:H, :);

    % - 2*vdp*param*ω^α * x * D^{α}x * Δx
    jb3(1).matrix = -2 * vdp_coefficient * parameter_coefficient * angular_frequency^fractional_order * x_series;
    jb3(2).matrix =      Dx_series;
    jb3(3).matrix =      da_series;
    t3 = multiply_series_coeffs(jb3); t3 = t3(1:H, :);

    s(1).matrix = t1; s(2).matrix = t2; s(3).matrix = t3;
    ja_for_dof  = sum_series_coeffs(s);
    ja_for_dof  = ja_for_dof(1:H, :);

    jacobian_wrt_coeffs_nonlinear_series(:, 2*dof_index-1 : 2*dof_index) = ja_for_dof;
end

%% --------------------- (3) Jw nonlinearity (w.r.t. ω) ----------------------
for dof_index = 1:num_dofs
    x_series   = disp_series_by_dof(dof_index).series;
    Dx_series  = frac_deriv_disp_by_dof(dof_index).series;

    fb(1).matrix = -vdp_coefficient * x_series;
    fb(2).matrix =  x_series;
    fb(3).matrix =  Dx_series;
    jw_for_dof   = multiply_series_coeffs(fb);
    jw_for_dof   = jw_for_dof(1:H, :);

    jacobian_wrt_frequency_nonlinear(:, 2*dof_index-1 : 2*dof_index) = jw_for_dof;
end

%% ----------------------------- Pack outputs --------------------------------
nonlinear_terms(1).part = residual_nonlinear_series;            % R (nonlinear)
nonlinear_terms(2).part = jacobian_wrt_coeffs_nonlinear_series; % Ja (nonlinear)
nonlinear_terms(3).part = jacobian_wrt_frequency_nonlinear;     % Jw (nonlinear)
end

function derivative_matrix = derive_series_coeffs(series_coeff_matrix, derivative_order)
%% Derive Fourier-series coefficient matrix for trigonometric terms
% Purpose:
%   Compute the n-th derivative of a Fourier-series coefficient matrix
%   with respect to time, for a representation:
%       x(t) = Σ [a_i*cos(iωt) + b_i*sin(iωt)]
%
% Inputs:
%   series_coeff_matrix : [H x (2*num_dofs)] matrix of Fourier coefficients
%                         Each DOF has two columns [cos_coeffs, sin_coeffs].
%                         Row 1 is DC term; rows 2..H correspond to harmonics.
%   derivative_order    : order of differentiation (integer ≥ 0)
%
% Outputs:
%   derivative_matrix   : coefficient matrix after differentiation
%                         Same dimension as input (H x 2*num_dofs)
%
% Method:
%   - Differentiation in time introduces factors of (iω)^n.
%   - For odd derivatives: cos→±sin, sin→∓cos.
%   - For even derivatives: cos→±cos, sin→±sin.
%   - DC term (row 1) has zero derivative for n ≥ 1.

% -------------------------------------------------------------------------

[num_harmonics, num_columns] = size(series_coeff_matrix);
derivative_matrix = zeros(num_harmonics, num_columns);

is_odd_derivative = mod(derivative_order, 2);  % 1=odd, 0=even

for harmonic_index = 1:num_harmonics
    for col_index = 1:num_columns
        % DC term has zero derivative when order >= 1
        if harmonic_index == 1 && derivative_order >= 1
            derivative_matrix(harmonic_index, col_index) = 0;
            continue;
        end

        harmonic_number = harmonic_index - 1; % actual harmonic order

        if is_odd_derivative
            % Odd derivative → cos↔sin with sign alternation
            % (-1) term switches between + and - depending on parity of column
            sign_factor = (-1)^((derivative_order - 3 + 2*(2 - mod(col_index,2))) / 2);
            paired_col = col_index + (-1)^(mod(col_index,2) + 1);
            derivative_matrix(harmonic_index, col_index) = ...
                (harmonic_number)^derivative_order * sign_factor * ...
                series_coeff_matrix(harmonic_index, paired_col);
        else
            % Even derivative → same function type, possible sign inversion
            sign_factor = (-1)^(derivative_order / 2);
            derivative_matrix(harmonic_index, col_index) = ...
                (harmonic_number)^derivative_order * sign_factor * ...
                series_coeff_matrix(harmonic_index, col_index);
        end
    end
end

% Ensure DC component of sine coefficients is zero (cos derivative vanishes)
for dof = 1:(num_columns/2)
    derivative_matrix(1, 2*dof) = 0;
end

end

function fractional_derivative_matrix = fractional_derive_series_coeffs(series_coeff_matrix, fractional_order)
%% Fractional derivative of Fourier-series coefficient matrix (series layout)
% Purpose
%   Compute D^{alpha}{x(t)} in the *coefficient space* for a trigonometric
%   Fourier expansion:
%       x(t) = sum_{k>=0} [ a_k * cos(k*.) + b_k * sin(k*.) ]
%   Here we operate directly on the coefficient "series matrix":
%       size(series_coeff_matrix) = [H  x  (2*num_dofs)]
%       per DOF: [cos_coeffs, sin_coeffs], row 1 is DC (k=0), rows 2..H are k=1..H-1.
%
% Inputs
%   series_coeff_matrix : [H x (2*num_dofs)] coefficients (DC row + harmonics)
%   fractional_order    : alpha in [0, 2]
%
% Output
%   fractional_derivative_matrix : same size as input, coefficients of D^{alpha}x
%
% Spectral rule used (kept consistent with your original codebase):
%   For harmonic k>=1:
%     c1 = k^alpha * cos(alpha*pi/2)
%     c2 = k^alpha * sin(alpha*pi/2)
%   Then, per DOF:
%     D^{alpha}{cos} ↔  c1 * cos  ±  c2 * sin
%     D^{alpha}{sin} ↔  c1 * sin  ∓  c2 * cos
%   DC (k=0) sine column remains 0.

% ----------------------------- argument checks -----------------------------
if fractional_order < 0 || fractional_order > 2
    warning('fractional_derive_series_coeffs:alphaOutOfRange', ...
        'fractional_order (alpha) should be within [0, 2].');
end

[num_harmonics, num_columns] = size(series_coeff_matrix);
fractional_derivative_matrix = zeros(num_harmonics, num_columns);

% ---------------------------- main (0 < alpha <= 1) ------------------------
% NOTE: Your original code used "if alpha>0 || alpha<=1" which is always true
% (for alpha in [0,2]) and made the else-branch unreachable. We keep the *same*
% behavior here (i.e., always execute this branch) to avoid any functional change.
for harmonic_idx = 1:num_harmonics
    k = harmonic_idx - 1;                 % harmonic order (k=0 is DC)
    c1 = k^fractional_order * cos(fractional_order*pi/2);
    c2 = k^fractional_order * sin(fractional_order*pi/2);

    for col_idx = 1:num_columns
        if mod(col_idx, 2) == 0
            % even column → sine coefficients
            % D^{alpha}{sin} = c1 * sin  -  c2 * cos
            fractional_derivative_matrix(harmonic_idx, col_idx) = ...
                c1 * series_coeff_matrix(harmonic_idx, col_idx) ...
              - c2 * series_coeff_matrix(harmonic_idx, col_idx - 1);
        else
            % odd column  → cosine coefficients
            % D^{alpha}{cos} = c1 * cos  +  c2 * sin
            fractional_derivative_matrix(harmonic_idx, col_idx) = ...
                c1 * series_coeff_matrix(harmonic_idx, col_idx) ...
              + c2 * series_coeff_matrix(harmonic_idx, col_idx + 1);
        end
    end
end

% ------------------------------ DC sine row --------------------------------
% Ensure DC (k=0) sine columns are exactly zero
for dof = 1:(num_columns/2)
    fractional_derivative_matrix(1, 2*dof) = 0;
end

% ----------------------------- historical note -----------------------------
% The original function contained an "else" branch with different algebra and
% references to undefined symbols (alpha1, Y1_N_harm_alpha1, etc.). Because the
% original condition was `if alpha>0 || alpha<=1`, that else-block was never
% reachable for alpha in [0,2]. To preserve behavior and avoid introducing bugs,
% we have *not* re-enabled that branch here. If you later need a distinct rule
% for alpha in (1,2], we can add it explicitly and update all callers together.

end

function series_matrix = columns_to_series_matrix(stacked_by_dof)
%% Rebuild series coefficient matrix from DOF-stacked column layout
% Purpose
%   Convert the stacked-by-DOF layout (each DOF is one column made by
%   stacking its cosine series followed by sine series-without-DC)
%   back to the standard series matrix layout:
%       [H x (2*num_dofs)] with per-DOF columns [cos, sin].
%
% Input
%   stacked_by_dof : [ (ceil(H/?) * 2 - 1) x num_dofs ]  layout produced by your solver:
%                    For each DOF column:
%                      - Top L rows are cosine coefficients (including DC).
%                      - Next (L-1) rows are sine coefficients (starting at k=1).
% Output
%   series_matrix  : [L x (2*num_dofs)] where
%                    - column 2*d-1 : cosine coefficients for DOF d
%                    - column 2*d   : sine   coefficients for DOF d (row 1 = 0)
%
% Notes
%   - This is the exact inverse of your packing convention where sine DC is omitted
%     and implicitly considered zero.

[num_rows_stacked, num_dofs] = size(stacked_by_dof);

% The input format uses L cosine rows and (L-1) sine rows stacked,
% so total rows per DOF column = L + (L-1) = 2L - 1  =>  L = ceil(num_rows_stacked/2)
num_harmonics = ceil(num_rows_stacked / 2);

series_matrix = zeros(num_harmonics, 2*num_dofs);

for d = 1:num_dofs
    % cosine coefficients (including DC at row 1)
    series_matrix(:, 2*d-1) = stacked_by_dof(1:num_harmonics, d);

    % sine coefficients start from harmonic 1 (row 2..end)
    if num_harmonics > 1
        series_matrix(2:end, 2*d) = stacked_by_dof(num_harmonics+1:end, d);
    end
    % row 1, sine column stays 0 by construction
end
end

function ihb_results = compute_ihb_residual_and_jacobians(coefficient_matrix, angular_frequency)
%% Compute residual R, Ja=dR/da, and Jw=dR/dω for the IHB linearization

global number_of_harmonics mass_matrix vdp_coefficient stiffness_matrix ...
       num_dofs fractional_order

H = number_of_harmonics;

% Zero increment used by the nonlinear builder for R term (interface stable)
zero_increment_coeffs = zeros(H, 2*num_dofs);

% Spectral derivatives (series-coefficient layout)
first_deriv_coeffs      = derive_series_coeffs(coefficient_matrix, 1); %#ok<NASGU>
second_deriv_coeffs     = derive_series_coeffs(coefficient_matrix, 2);
fractional_deriv_coeffs = fractional_derive_series_coeffs(coefficient_matrix, fractional_order);

% Nonlinear contributions (R / Ja / Jw)
nonlinear_terms = build_nonlinear_terms(angular_frequency, coefficient_matrix, zero_increment_coeffs);

%% ------------------------------- Residual R --------------------------------
% residual_series = -[ M*(ω^2 x'') + K*x + vdp*(ω^α D^α x) ] - nonlinear_terms(1)
rblk(1).matrix = apply_linear_operator_to_series(mass_matrix,     -angular_frequency^2          * second_deriv_coeffs);
rblk(2).matrix = apply_linear_operator_to_series(vdp_coefficient, -angular_frequency^fractional_order * fractional_deriv_coeffs);
rblk(3).matrix = apply_linear_operator_to_series(stiffness_matrix, -coefficient_matrix);
rblk(4).matrix = -nonlinear_terms(1).part;

residual_series = sum_series_coeffs(rblk);
residual_series = residual_series(1:H, :);  % trim to H rows

residual_packed = pack_series_coeffs_per_dof(residual_series);
residual_vector = stack_columns_vector(residual_packed);

%% ---------------------------- Ja: dR/d(coeffs) -----------------------------
jacobian_wrt_coeffs = assemble_coefficients_jacobian(coefficient_matrix, angular_frequency);

%% ----------------------------- Jw: dR/dω -----------------------------------
% d/dω[ M*(ω^2 x'') ] = 2ω M x''
% d/dω[ vdp*(ω^α D^α x) ]  →  vdp * D^α x  (kept consistent with original pipeline)
wblk(1).matrix = apply_linear_operator_to_series(mass_matrix,     2*angular_frequency * second_deriv_coeffs);
wblk(2).matrix = apply_linear_operator_to_series(vdp_coefficient, fractional_deriv_coeffs);
wblk(3).matrix = nonlinear_terms(3).part;

freq_jac_series = sum_series_coeffs(wblk);
freq_jac_series = freq_jac_series(1:H, :);  % trim

freq_packed            = pack_series_coeffs_per_dof(freq_jac_series);
jacobian_wrt_frequency = stack_columns_vector(freq_packed);

%% ------------------------------ Pack outputs -------------------------------
ihb_results(1).vector = residual_vector;          % R
ihb_results(2).vector = jacobian_wrt_coeffs;      % Ja
ihb_results(3).vector = jacobian_wrt_frequency;   % Jw
end

function series_out = harmonic_product_convolution(series_A, series_B)
%% Spectral convolution for time-domain product of two series coefficient matrices
% Inputs
%   series_A, series_B : [H x (2*num_dofs)] series coeffs (per DOF: [cos, sin])
% Output
%   series_out         : [H_A+H_B-1 x (2*num_dofs)] convolved series coeffs
%
% Implements trigonometric product identities in coefficient space:
%   cos⋅cos, cos⋅sin, sin⋅cos, sin⋅sin → four blocks combined (n1..n4 in original)
% Notes
%   - Uses sign_nonzero(k) which returns sign(k) with sign(0)=0.
%   - DC row sine is kept implicit by construction elsewhere.

len_A   = size(series_A, 1);
len_B   = size(series_B, 1);
num_dofs = size(series_A, 2) / 2;
len_out = len_A + len_B - 1;

% Preallocate four partial blocks for (cos*cos), (cos*sin), (sin*cos), (sin*sin)
cc_block = zeros(len_out, 2*num_dofs);
cs_block = zeros(len_out, 2*num_dofs);
sc_block = zeros(len_out, 2*num_dofs);
ss_block = zeros(len_out, 2*num_dofs);

%% 1) cos*cos → cos(sum) + cos(diff)
for i = 1:len_A
    for j = 1:len_B
        for d = 1:num_dofs
            cc_block(i+j-1, 2*d-1) = cc_block(i+j-1, 2*d-1) + 0.5 * series_A(i, 2*d-1) * series_B(j, 2*d-1);
            cc_block(abs(i-j)+1, 2*d-1) = cc_block(abs(i-j)+1, 2*d-1) + 0.5 * series_A(i, 2*d-1) * series_B(j, 2*d-1);
        end
    end
end

%% 2) cos*sin → sin(sum) - sign(diff)*sin(diff)
% j starts at 2 because sin(0)=0 (first row sine column is zero)
for i = 1:len_A
    for j = 2:len_B
        for d = 1:num_dofs
            cs_block(i+j-1, 2*d) = cs_block(i+j-1, 2*d) + 0.5 * series_A(i, 2*d-1) * series_B(j, 2*d);
            cs_block(abs(i-j)+1, 2*d) = cs_block(abs(i-j)+1, 2*d) ...
                - 0.5 * sign_nonzero(i-j) * series_A(i, 2*d-1) * series_B(j, 2*d);
        end
    end
end

%% 3) sin*cos → sin(sum) + sign(diff)*sin(diff)
for i = 2:len_A
    for j = 1:len_B
        for d = 1:num_dofs
            sc_block(i+j-1, 2*d) = sc_block(i+j-1, 2*d) + 0.5 * series_A(i, 2*d) * series_B(j, 2*d-1);
            sc_block(abs(i-j)+1, 2*d) = sc_block(abs(i-j)+1, 2*d) ...
                + 0.5 * sign_nonzero(i-j) * series_A(i, 2*d) * series_B(j, 2*d-1);
        end
    end
end

%% 4) sin*sin → -cos(sum) + cos(diff)
for i = 2:len_A
    for j = 2:len_B
        for d = 1:num_dofs
            ss_block(i+j-1, 2*d-1) = ss_block(i+j-1, 2*d-1) - 0.5 * series_A(i, 2*d) * series_B(j, 2*d);
            ss_block(abs(i-j)+1, 2*d-1) = ss_block(abs(i-j)+1, 2*d-1) + 0.5 * series_A(i, 2*d) * series_B(j, 2*d);
        end
    end
end

% Combine all partial blocks
series_out = cc_block + cs_block + sc_block + ss_block;

end

function product_series = multiply_series_coeffs(series_operands)
%% Multiply multiple Fourier series (in coefficient-matrix form) via convolution
% Input
%   series_operands : struct array with field .matrix
%                     each .matrix is a series coefficient matrix:
%                       size = [H x (2*num_dofs)], columns per DOF = [cos, sin]
% Output
%   product_series  : resulting series coefficient matrix after sequential multiplications
%
% Behavior
%   ((...((A ⊗ B) ⊗ C) ⊗ D) ... ), where ⊗ denotes spectral convolution
%   implemented by harmonic_product_convolution.

num_operands = numel(series_operands);
product_series = series_operands(1).matrix;
for k = 2:num_operands
    product_series = harmonic_product_convolution(product_series, series_operands(k).matrix);
end
end

function packed_cols = pack_series_coeffs_per_dof(series_matrix)
%% Pack per-DOF Fourier series coefficients into stacked columns
% Purpose
%   Convert [H x (2*num_dofs)] series matrix (per DOF: [cos, sin])
%   into [2H-1 x num_dofs] stacked columns:
%     column d = [cos_coeffs(:,d); sin_coeffs(2:end,d)]
%
% Input
%   series_matrix : [H x (2*num_dofs)] with [cos, sin] for each DOF
%
% Output
%   packed_cols   : [2*H-1 x num_dofs], per-DOF stacked column

[num_harmonics, two_ndofs] = size(series_matrix);
num_dofs = two_ndofs / 2;

packed_cols = zeros(2*num_harmonics - 1, num_dofs);
for dof = 1:num_dofs
    cos_col = series_matrix(:, 2*dof - 1);
    sin_col = series_matrix(:, 2*dof    );  % sin DC = row1 (usually 0 by convention)
    packed_cols(1:num_harmonics, dof)               = cos_col;
    packed_cols(num_harmonics+1:end, dof)           = sin_col(2:end);
end
end

function s = sign_nonzero(x)
%% Sign function with sign_nonzero(0) = 0
% Returns:
%   -1 if x < 0,  0 if x == 0,  +1 if x > 0

if x == 0
    s = 0;
else
    s = x / abs(x);
end
end

function colvec = stack_columns_vector(matrix_in)
%% Stack a 2D matrix column-wise into a single column vector
% Input : matrix_in [m x n]
% Output: colvec    [m*n x 1]  (column-major stacking)
[m, n] = size(matrix_in);
colvec = zeros(m*n, 1);
for j = 1:n
    colvec( (j-1)*m + 1 : j*m, 1 ) = matrix_in(:, j);
end
end

function summed_series = sum_series_coeffs(series_list)
%% Sum multiple Fourier series coefficient matrices (cos/sin paired by DOF)
% Input
%   series_list : struct array with field .matrix
%                 each .matrix is [H x (2*num_dofs)] with per-DOF columns [cos, sin]
% Output
%   summed_series : the sum over all .matrix, row-aligned (zero-padded as needed)
%
% Notes
%   - Handles different row counts by expanding to the max row length.
%   - Preserves the convention: DC sine (row 1, even columns) set to 0.

num_terms = numel(series_list);
max_rows  = 0;
num_cols  = size(series_list(1).matrix, 2);

% find max number of rows among inputs
for k = 1:num_terms
    max_rows = max(max_rows, size(series_list(k).matrix, 1));
end

summed_series = zeros(max_rows, num_cols);

% accumulate using the two-matrix adder (handles different heights)
for k = 1:num_terms
    summed_series = add_series_coeff_mats(summed_series, series_list(k).matrix);
end
end

function matrix_out = unstack_columns_by_dof(vector_in)
%% Unstack a vector back to [n_rows x num_dofs] by DOF columns
% Purpose
%   Inverse of stack_columns_vector (for DOF-wise column matrices).
% Inputs
%   vector_in  : [(n_rows*num_dofs) x 1]
% Globals
%   num_dofs   : number of DOFs
% Output
%   matrix_out : [n_rows x num_dofs]

global num_dofs
L  = numel(vector_in);
nr = L / num_dofs;

matrix_out = zeros(nr, num_dofs);
for dof = 1:num_dofs
    matrix_out(:, dof) = vector_in( 1 + (dof-1)*nr : dof*nr );
end
end

function series_out = add_scalar_to_series_dc(scalar_value, series_in)
%% Add a scalar to the DC cosine entry of a series coefficient matrix
% Inputs
%   scalar_value : scalar to add
%   series_in    : [H x (2*num_dofs)] series coeffs (DC row = row 1)
% Output
%   series_out   : same as series_in, with (1,1) increased by scalar_value

series_in(1,1) = series_in(1,1) + scalar_value;
series_out = series_in;
end

function out_series = add_series_coeff_mats(series_A, series_B)
%% Add two Fourier series coefficient matrices with row alignment
% Inputs
%   series_A, series_B : [H_A x (2*num_dofs)], [H_B x (2*num_dofs)]
% Output
%   out_series         : [max(H_A,H_B) x (2*num_dofs)] element-wise sum
%
% Behavior
%   - Pads the shorter matrix in rows.
%   - Keeps per-DOF [cos, sin] layout.
%   - Enforces DC sine to be zero in the output (row 1, even columns).

rows_A  = size(series_A, 1);
rows_B  = size(series_B, 1);
num_dofs = size(series_A, 2) / 2;

rows_out   = max(rows_A, rows_B);
out_series = zeros(rows_out, 2*num_dofs);

% ensure A is the taller (swap if needed)
if rows_A < rows_B
    tmp      = series_B;
    series_B = series_A;
    series_A = tmp;
    rows_B   = rows_A;
    rows_A   = rows_out;
end

% rows where both have data
if rows_B > 0
    out_series(1:rows_B, :) = series_A(1:rows_B, :) + series_B(1:rows_B, :);
end

% remaining rows from the taller A
if rows_A > rows_B
    out_series(rows_B+1:rows_A, :) = series_A(rows_B+1:rows_A, :);
end

% enforce DC sine (row 1, even columns) = 0
for d = 1:num_dofs
    out_series(1, 2*d) = 0;
end
end

function series_out = apply_linear_operator_to_series(operator_matrix, series_in)
%% Apply a linear operator (scalar or square matrix) to a series coefficient matrix
% Inputs
%   operator_matrix : scalar or [num_dofs x num_dofs] matrix (e.g., mass, stiffness)
%   series_in       : [H x (2*num_dofs)] series coeffs, columns per DOF = [cos, sin]
% Output
%   series_out      : same size as series_in; operator applied DOF-wise
%
% Notes
%   For each DOF block (two columns: cos/sin), we form a linear combination
%   with operator rows. This preserves the [cos, sin] pairing per DOF.

[num_harmonics, total_cols] = size(series_in);
num_dofs = total_cols / 2;

% Split per DOF into 2-column blocks
for d = 1:num_dofs
    dof_block(d).part = series_in(:, 2*d-1:2*d); %#ok<AGROW>
end

% First column combination
for i = 1:num_dofs
    series_blocks(i).part = operator_matrix(i,1) * dof_block(1).part; %#ok<AGROW>
end

% Accumulate remaining columns
for i = 1:num_dofs
    for j = 1:(num_dofs - 1)
        series_blocks(i).part = add_series_coeff_mats( ...
            series_blocks(i).part, operator_matrix(i, j+1) * dof_block(j+1).part);
    end
end

% Stitch back to full series matrix
series_out = zeros(num_harmonics, 2*num_dofs);
for d = 1:num_dofs
    series_out(:, 2*d-1:2*d) = series_blocks(d).part;
end
end
