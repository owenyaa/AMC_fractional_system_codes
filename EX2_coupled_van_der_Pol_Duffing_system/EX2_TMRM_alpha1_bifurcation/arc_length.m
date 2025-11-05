function Predictor=arc_length(every_a)
s = zeros(1, 3); % Arc lengths
for j = 1:3
    s(j) = sqrt(norm(every_a(j+1).harmonic_coefficients - every_a(j).harmonic_coefficients));
    if s(j) < 1e-15
        s(j) = 1e-15; % Avoid zero arc length
    end
end

% Step 2: Define time parameters t based on arc lengths
t = zeros(5, 1); % Time parameters
t(1) = 0;
t(2) = s(1);
t(3) = t(2) + s(2);
t(4) = t(3) + s(3);

%ds(j)=ds(j-1)*Nd/I(j-1),Nd为一个与非线性频率-振幅响应曲线相关的整数，一般取4或5,
%I(j-1)为完成上一次计算所进行的迭代次数  具体参见张丹伟博士论文P51页
ds=s(3)+(s(3)-s(2));
t(5)=t(4)+ds;

Predictor = zeros(size(every_a(1).harmonic_coefficients)); % Initialize Predictor
for i = 1:4
    L = 1; % Initialize Lagrange basis
    for j = 1:4
        if i ~= j
            L = L * (t(5) - t(j)) / (t(i) - t(j));
        end
    end
    Predictor = Predictor + L * every_a(i).harmonic_coefficients;
end
end