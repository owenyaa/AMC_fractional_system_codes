%%=========================================================================
%  compute_J_alpha1  ——  Residual tail J^{alpha1}(t) for fractional operator
%%==========================================================================
% Purpose (做什么)
%   Given a semi-analytical quasi-periodic/periodic approximation from the
%   “时域最小残值法”(Time-Domain Minimum Residual Method, TMRM), this routine
%   computes the residual tail of the fractional-order operator, J^{alpha1}(t),
%   defined in the paper’s Eq. (28) for 0<alpha1<1, using numerically stable
%   log-domain accumulation. We return log(J^{alpha1}(t)) to prevent overflow.
%
%   核心目标：在获得 TMRM 半解析解的谐波系数 {b_{jk}, c_{jk}} 与基频 ω 后，
%   依据论文公式 (28) 评估分数阶算子 D^{alpha1} x_j(t) 的“非周期余项”J^{alpha1}(t)，
%   并通过对数域稳定求和避免 e^{kωt} 与 K_ν(kωt) 的数值灾难，从而用于验证
%   J^{alpha1}(t) ~ t^{-1/2} 的幂律衰减与 S-渐近 T-周期性。
%
% Mathematical background (数学背景——对应论文公式(28))
%   For 0<alpha1<1, the residual tail reads (paper Eq. (28)):
%
%     J^{alpha1}(t) = [sin(pi*alpha1)/(2*pi)] * sum_{k=1}^∞ (k*ω)^{alpha1} e^{kωt} * ...
%        { c_{jk} * Γ((alpha1+1)/2) * K_{(1 - alpha1)/2}(kωt) ...
%        + b_{jk} * Γ(alpha1/2)     * K_{1 - alpha1/2}(kωt) } .
%
%   Asymptotically, K_ν(z) ~ sqrt(pi/(2z)) * e^{-z}, so e^{kωt}*K_ν(kωt) ~ O(t^{-1/2}),
%   hence J^{alpha1}(t) decays like t^{-1/2}. This justifies computing/plotting
%   exp(log_J_alpha1) against t^{-1/2} in the main script for validation.
%
% Numerical strategy (数值策略——为何取对数与带符号求和)
%   • Direct summation in (28) is ill-conditioned because e^{kωt} can overflow
%     while K_ν(kωt) underflows; their product is well-scaled but not in finite
%     precision. We therefore:
%       - Assemble per-harmonic terms in the LOG domain;
%       - Track the algebraic sign from the linear combination
%         [c_{jk}*Γ1*K_{(1-alpha1)/2} + b_{jk}*Γ2*K_{1-alpha1/2}];
%       - Use a signed log-sum-exp reducer to combine many terms robustly;
%       - Add the constant log-factor: log(sin(pi*alpha1)) - log(2*pi).
%   • We compute LOG(J^{alpha1}) directly to avoid overflow/underflow; the caller
%     can exponentiate if needed, e.g., for log-log plots vs. t^{-1/2}.
%
% Inputs (输入)
%   time_data          : 1×Nt vector of time samples t >= 0 used for evaluation.
%   alpha1             : fractional order in (0,1).
%   initial_frequency  : scalar ω (the base frequency from the TMRM solution).
%   harmonic_parameters: N×(2*nDOF) matrix of Fourier coefficients from TMRM;
%                        for DOF j, column (2*j-1) stores c_{jk}, column (2*j)
%                        stores b_{jk}, each row k=1..N the k-th harmonic.
%   j                  : target DOF index (1-based) for which J^{alpha1}(t) is computed.
%   num_harmonics      : number of retained harmonics N used in the truncated sum.
%
% Outputs (输出)
%   log_J_alpha1       : 1×Nt vector, log( J^{alpha1}(t) ) evaluated at time_data.
%                        Use exp(log_J_alpha1) when you need J^{alpha1}(t).
%
% Conventions & units (约定与单位)
%   • ω = initial_frequency; harmonic ν_k = k*ω.
%   • c_{jk} multiplies cos(kωt), b_{jk} multiplies sin(kωt) in your main code.
%   • Γ1 = Γ((alpha1+1)/2), Γ2 = Γ(alpha1/2).
%   • K_ν(·) is the modified Bessel function of the 2nd kind (Matlab: besselk).
%
% Algorithm (求解流程)
%   1) Precompute constant terms in LOG domain:
%        log_factor = log(sin(pi*alpha1)) - log(2*pi);
%        Γ1 = gamma((alpha1+1)/2); Γ2 = gamma(alpha1/2).
%   2) For k = 1..N:
%        ν = k*ω; z = ν * time_data;
%        read c_{jk} = harmonic_parameters(k, 2*j - 1);
%             b_{jk} = harmonic_parameters(k, 2*j);
%        evaluate K1 = K_{(1-alpha1)/2}(z), K2 = K_{1 - alpha1/2}(z);
%        form B_k(t) = c_{jk}*Γ1*K1 + b_{jk}*Γ2*K2   (track sign(B_k));
%        accumulate per-harmonic log-term:
%             log_term_k(t) = alpha1*log(ν) + z + log(|B_k(t)|).
%   3) Reduce across k via a signed log-sum-exp:
%        Let s_k(t) = sign(B_k(t)); perform
%           LOG |Σ_k s_k * exp(log_term_k - max_log(t))| + max_log(t),
%        with max_log(t) = max_k log_term_k(t). Handle exact cancellation by
%        tagging points where net sum is ≈ 0.
%   4) Add the global factor:
%        log_J_alpha1(t) = logsum(t) + log_factor.
%   5) Return log_J_alpha1. (Caller may plot exp(log_J_alpha1) vs t^{-1/2}.)
%
% Stability & edge cases (稳定性与边界情形)
%   • Very small |B_k(t)| → use floor epsilon before log to avoid -Inf,
%     but keep sign tracking for cancellation patterns.
%   • Large t or ω: z = ν t can be huge; besselk handles large-argument decay,
%     yet log-domain accumulation is still mandatory for stability.
%   • If positive/negative contributions cancel to machine zero at some t,
%     we mark those time slots (internally) to avoid spurious -Inf; the caller
%     may interpret them as “effectively zero residual tail”.
%
% Complexity (复杂度)
%   O(Nt * N) Bessel calls + O(Nt * N) elementary ops; memory O(Nt).
%
% Validation checklist (数值验证要点)
%   • As t grows, exp(log_J_alpha1) should align with Ct^{-1/2} in log-log
%     scale (slope ≈ -1/2), consistent with the paper’s asymptotics.
%   • Sensitivity to N: enlarge num_harmonics until the curve stabilizes.
%
% Integration with TMRM (与主程序的接口)
%   • The function consumes the TMRM-updated ω, b_{jk}, c_{jk}; call it after
%     current iteration converges (or at the end) to diagnose the residual tail.
%   • For multi-DOF, loop over j if needed; keep columns [2*j-1, 2*j].
%
% References
%   • See paper Eq. (28) and surrounding derivation for J^{alpha1}(t), and
%     asymptotic analysis showing t^{-1/2} decay of the residual tail.
%
% Author: GUANG_LIU  * owenyaa@gmail.com *
% First created: 2020-11-12 ; Last revised: 2025-11-05
%==========================================================================

function log_J_alpha1 = compute_J_alpha1(time_data, alpha1, initial_frequency, harmonic_parameters, j, num_harmonics)
% 稳定计算 log(J^{alpha1}(t))，保留符号追踪并避免 NaN 问题

log_factor = log(sin(pi * alpha1)) - log(2 * pi);
gamma1 = gamma((alpha1 + 1) / 2);
gamma2 = gamma(alpha1 / 2);
epsilon = 1e-400;

Nt = length(time_data);
log_terms = -Inf(num_harmonics, Nt);  % 每个谐波对应一个 log 值
sign_terms = zeros(num_harmonics, Nt);  % 记录正负符号 (+1 / -1)

for k = 1:num_harmonics
    nu = k * initial_frequency;
    c_jk = harmonic_parameters(k, 2*j - 1);
    b_jk = harmonic_parameters(k, 2*j);
    z = nu * time_data;

    % 计算 Bessel 函数项
    K1 = besselk((1 - alpha1)/2, z);
    K2 = besselk(1 - alpha1/2, z);

    % 构造 B_k(t)
    B_k = c_jk * gamma1 .* K1 + b_jk * gamma2 .* K2;

    % 对于极小的 B_k，避免 log(0)
    abs_B_k = abs(B_k);
    abs_B_k(abs_B_k < epsilon) = epsilon;

    % 对数项
    log_B_k = log(abs_B_k);
    common_log = alpha1 * log(nu) + z + log_B_k;

    % 符号记录：正数记 +1，负数记 -1
    sign_terms(k, :) = sign(B_k);
    log_terms(k, :) = common_log;
end

% 总和：按符号分组进行 logsumexp
[log_J_alpha1, nan_mask] = signed_logsumexp(log_terms, sign_terms);

% 乘上 log_factor，最终结果
log_J_alpha1 = log_J_alpha1 + log_factor;

% 避免不可定义情况
log_J_alpha1(nan_mask) = NaN;
end

% ------------- 带符号的 logsumexp ----------------
function [logsum, nan_mask] = signed_logsumexp(logX, signs)
% logX: M x N log 值矩阵
% signs: M x N (+1 or -1)
% 返回 logsum：log(abs(总和))，以及 nan_mask：0值点

max_log = max(logX, [], 1);  % 对每个时间点求最大项，避免下溢
X_shifted = exp(logX - max_log);  % 变成非爆炸项

% 分别加总正负
pos_sum = sum(X_shifted .* (signs > 0), 1);
neg_sum = sum(X_shifted .* (signs < 0), 1);
net_sum = pos_sum - neg_sum;

% 符号判断（最终结果正还是负）
valid = abs(net_sum) > 1e-300;
logsum = -Inf(size(net_sum));
logsum(valid) = log(abs(net_sum(valid))) + max_log(valid);

% 返回哪些时间点由于完全抵消而导致 NaN（如 net_sum = 0）
nan_mask = ~valid;
end