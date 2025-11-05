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
