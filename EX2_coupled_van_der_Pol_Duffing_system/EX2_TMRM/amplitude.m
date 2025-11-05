clear; clc; close all;
load para_alpha_0.8_0.6.mat;
for i=1:1:length(0.8:-0.002:0.6)
    harmonic_coefficients=every_harmonic_coefficients(i).harmonic_coefficients;
    a1=harmonic_coefficients(2:end,1);b1=harmonic_coefficients(2:end,2);
    a2=harmonic_coefficients(2:end,3);b2=harmonic_coefficients(2:end,4);
    A_total_x1(i) = sqrt(sum(a1.^2 + b1.^2));
    A_total_x2(i) = sqrt(sum(a2.^2 + b2.^2));
    alpha(i)=every_harmonic_coefficients(i).alpha;
end
figure;
plot(alpha, A_total_x1, 'k-', 'LineWidth', 1);
% hold on;
% plot(Omegas, A_total_x2, 'b-', 'LineWidth', 1);
% legend_handle = legend('$$x_1$$', '$$x_2$$');
% set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);




