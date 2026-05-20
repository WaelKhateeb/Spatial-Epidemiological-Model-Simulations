%% Define symbolic variables and parameters
clear; clc;
syms I mu1_sym A_sym

% Parameters
beta     = 0.1;
d        = 0.1;
mu0      = 0.1;
mu1_val  = 0.2;
b        = 0.3;
A_val    = 0.4;
l_val    = 5;
r2_val   = 0.01;

% Define h(I) and its derivative h'(I)
h       = mu0 + (mu1_sym - mu0)*b/(I + b);
h_prime = diff(h, I);

% Solve for I from: beta * A/(d + beta*I) - d - h(I) = 0
S_expr = A_sym / (d + beta * I);
eq     = beta * S_expr - d - h == 0;
I_sol  = solve(eq, I);

% Substitute numeric values
I_sol_num = double(subs(I_sol, [mu1_sym, A_sym], [mu1_val, A_val]));

% Filter for real, positive roots
tol = 1e-8;
real_pos_idx = (abs(imag(I_sol_num)) < tol) & (real(I_sol_num) > tol);
if ~any(real_pos_idx)
    error('No positive real endemic equilibria found.');
end

% Among those, pick the one with the largest I-value
valid_I_vals = real(I_sol_num(real_pos_idx));
[~, idx_max] = max(valid_I_vals);
valid_roots  = I_sol(real_pos_idx);
I2_sol       = valid_roots(idx_max);

% Compute Sigma values
Sigma1 = -beta * I2_sol - d;
Sigma2 = -subs(h, I, I2_sol) - d;
Sigma3 = beta * I2_sol;
Sigma4 = -I2_sol * subs(h_prime, I, I2_sol);

Sigma1_num = double(subs(Sigma1, [mu1_sym, A_sym], [mu1_val, A_val]));
Sigma2_num = double(subs(Sigma2, [mu1_sym, A_sym], [mu1_val, A_val]));
Sigma3_num = double(subs(Sigma3, [mu1_sym, A_sym], [mu1_val, A_val]));
Sigma4_num = double(subs(Sigma4, [mu1_sym, A_sym], [mu1_val, A_val]));

%% Define r1(k) function
r1_func = @(k) (l_val^4 * (Sigma1_num * Sigma4_num - Sigma2_num * Sigma3_num - Sigma1_num * r2_val * (k.^2 / l_val^2))) ...
             ./ (k.^2 .* (Sigma4_num * l_val^2 - r2_val * k.^2));

%% Compute asymptote value
k_VA = sqrt(Sigma4_num * l_val^2 / r2_val);  % Vertical Asymptote
fprintf('Vertical asymptote at k = %.4f\n', k_VA);

%% Define domain up to just before asymptote
k_real = linspace(0, k_VA - 1e-10, 1e6);      % ℝ domain
k_nat  = 1:floor(k_VA - 0.01);                % ℕ domain

% Evaluate r1
r1_real_vals = r1_func(k_real);
r1_nat_vals  = r1_func(k_nat);

% Clean up undefined values
r1_real_vals(~isfinite(r1_real_vals)) = NaN;
r1_nat_vals(~isfinite(r1_nat_vals))   = NaN;

%% Plotting
figure('Position', [100, 100, 900, 500]);
hold on;

plot(k_real, r1_real_vals, 'b-', 'LineWidth', 2, 'DisplayName', 'r_1(k) for k ∈ ℝ');
plot(k_nat,  r1_nat_vals,  'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'r_1(k) for k ∈ ℕ');

xline(k_VA, ':k', 'LineWidth', 1.5, 'HandleVisibility', 'off');  % asymptote

xlabel('k');
ylabel('r_1(k)');
title(sprintf('Plot of r_1 vs k (only up to VA at k = %.4f), r_2 = %.2f', k_VA, r2_val));
legend('Location', 'best');
legend('$r_1(k)$ for $k \in$ ℝ', '$r_1(k)$ for $k \in$ ℕ', 'Interpreter', 'latex');

grid on;
hold off;
