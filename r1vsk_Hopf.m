%% Pick the positive endemic root with the largest I-value
clear; clc;
syms I mu1_sym A_sym

% Parameter values
beta     = 0.009407413516;     % Transmission rate
d        = 0.01;               % Natural death rate
mu0      = 0.1;                % Minimum disease‐induced mortality
mu1_val  = 10;                 % Maximum disease‐induced mortality
b        = 0.03;               % Half‐saturation constant
A_val    = 1;                  % Recruitment rate
r2_val   = 0.01;               % Diffusion coefficient for I
l_val    = 5;                  % Domain length [0,L]

% Define h(I) and its derivative h'(I)
h        = mu0 + (mu1_sym - mu0)*b/(I + b);
h_prime  = diff(h, I);

% Equilibrium condition: beta * A/(d + beta*I) - d - h(I) = 0
S_expr   = A_sym / (d + beta * I);
eq       = beta * S_expr - d - h == 0;

% Solve symbolically for I
I_sol    = solve(eq, I);

% Substitute numeric values and pick largest positive real root
I_sol_num = double(subs(I_sol, [mu1_sym, A_sym], [mu1_val, A_val]));
valid     = (imag(I_sol_num)==0) & (I_sol_num > eps);

if ~any(valid)
    error('No positive endemic equilibria found.');
end

% Among the valid roots, find the one with the maximal I-value
valid_I_vals = I_sol_num(valid);
[~, max_idx] = max(valid_I_vals);
valid_roots  = I_sol(valid);
I2_sol       = valid_roots(max_idx);

% Compute the four Sigma entries
Sigma1 = -beta * I2_sol - d;
Sigma2 = -subs(h, I, I2_sol) - d;
Sigma3 = beta * I2_sol;
Sigma4 = -I2_sol * subs(h_prime, I, I2_sol);

% Numeric evaluation of Sigmas
Sigma1_num = double(subs(Sigma1, [mu1_sym, A_sym], [mu1_val, A_val]));
Sigma2_num = double(subs(Sigma2, [mu1_sym, A_sym], [mu1_val, A_val]));
Sigma3_num = double(subs(Sigma3, [mu1_sym, A_sym], [mu1_val, A_val]));
Sigma4_num = double(subs(Sigma4, [mu1_sym, A_sym], [mu1_val, A_val]));

%% Define r1(k) function
r1_func = @(k) (l_val^4 * (Sigma1_num * Sigma4_num - Sigma2_num * Sigma3_num - Sigma1_num * r2_val * (k.^2 / l_val^2))) ...
             ./ (k.^2 .* (Sigma4_num * l_val^2 - r2_val * k.^2));

%% Compute vertical asymptote
k_VA = sqrt(Sigma4_num * l_val^2 / r2_val);
fprintf('Vertical asymptote at k = %.4f\n', k_VA);

%% Define domains
k_real = linspace(0, k_VA - 1e-10, 1e6);  % ℝ domain up to asymptote
k_nat  = 1:floor(k_VA - 0.01);            % ℕ domain within safe range

% Evaluate r1
r1_real_vals = r1_func(k_real);
r1_nat_vals  = r1_func(k_nat);

% Clean up non‐finite values
r1_real_vals(~isfinite(r1_real_vals)) = NaN;
r1_nat_vals(~isfinite(r1_nat_vals))   = NaN;

%% Plotting with axes limited to [0,10] x [0,10]
figure('Position', [100, 100, 900, 500]); hold on;
plot(k_real, r1_real_vals, 'b-', 'LineWidth', 2, 'DisplayName', 'r_1(k), k ∈ ℝ');
plot(k_nat,  r1_nat_vals,  'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'r_1(k), k ∈ ℕ');

% Asymptote line
xline(k_VA, ':k', 'LineWidth', 1.5, 'HandleVisibility','off');

% Limit axes to 10x10
axis([0 10 0 10]);

xlabel('k');
ylabel('r_1(k)');
title(sprintf('r_1 vs k (axes limited to [0,10]×[0,10]), r_2 = %.2f', r2_val));
legend('Location','best','Interpreter','none');
grid on;
hold off;
