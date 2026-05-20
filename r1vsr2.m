%% Define parameters (all numeric)
clear; clc;
beta    = 0.1;      % Transmission rate
d       = 0.1;      % Natural death rate
mu0     = 0.1;      % Minimum disease‐induced mortality
mu1_val = 0.2;      % Maximum disease‐induced mortality
b       = 0.3;      % Half‐saturation constant
A_val   = 0.4;      % Recruitment rate
l_val   = 5;        % Domain length
    

%% Solve the endemic condition numerically:
%    beta * A/(d + beta*I) - d - [mu0 + (mu1-mu0)*b/(I+b)] = 0
syms I
eq_num = beta*A_val/(d + beta*I) - d - (mu0 + (mu1_val - mu0)*b/(I + b)) == 0;

% find all roots
I_roots_sym = solve(eq_num, I);
I_roots     = double(I_roots_sym);

% pick only real, positive
valid = (imag(I_roots)==0) & (I_roots > eps);
if ~any(valid)
    error('No positive endemic equilibria found.');
end

I_vals = I_roots(valid);
I2     = max(I_vals);    % the endemic I* with highest value
fprintf('Selected endemic I* = %.6f\n', I2);

%% Compute the Sigma entries numerically
h       = @(I) mu0 + (mu1_val - mu0)*b./(I + b);
h_prime = @(I) - (mu1_val - mu0)*b ./ (I + b).^2;  % derivative of h

Sigma1 = -beta * I2         - d;
Sigma2 = -h(I2)             - d;
Sigma3 =  beta * I2;
Sigma4 = -I2 * h_prime(I2);

fprintf('Sigma1 = %.4f\n', Sigma1);
fprintf('Sigma2 = %.4f\n', Sigma2);
fprintf('Sigma3 = %.4f\n', Sigma3);
fprintf('Sigma4 = %.4f\n', Sigma4);

%% Compute R1 and the line r1 = r2 / R1
numerator_R1   = Sigma1*Sigma4 ...
               - 2*Sigma2*Sigma3 ...
               - 2*sqrt(Sigma2*Sigma3*(Sigma2*Sigma3 - Sigma1*Sigma4));
denominator_R1 = Sigma1^2;
R1             = numerator_R1 / denominator_R1;
fprintf('R1 = %.4f\n', R1);

r1_line = @(r2) r2 / R1;

%% Define the mode functions r_k(r2)
r_k = @(r2, k) ...
    (l_val^4 * (Sigma1*Sigma4 - Sigma2*Sigma3 - Sigma1*r2*(k^2/l_val^2))) ...
   ./ (k^2 * (Sigma4*l_val^2 - r2*k^2));

%% Plot everything
r2_vals        = linspace(0,0.5,2000);
largeThreshold = 1e3;
colors         = lines(6);

figure('Position',[100,100,900,550]); hold on;
for k = 1:5
    r_vals = r_k(r2_vals, k);
    % clean up
    bad = ~isfinite(r_vals) | r_vals<=0 | abs(r_vals)>largeThreshold;
    r_vals(bad) = NaN;
    plot(r2_vals, r_vals, 'Color',colors(k,:), 'LineWidth',1.5, ...
         'DisplayName',sprintf('r_{%d}(r_2)',k));
    % asymptote
    r2_asym = (l_val^2 * Sigma4)/(k^2);
    if r2_asym>0 && r2_asym<=max(r2_vals)
        xline(r2_asym,':','Color',colors(k,:), 'LineWidth',1.2, ...
              'DisplayName',sprintf('Asymptote k=%d',k));
    end
end

% plot r1 = r2 / R1
r1_vals = r1_line(r2_vals);
bad = ~isfinite(r1_vals) | r1_vals<=0 | abs(r1_vals)>largeThreshold;
r1_vals(bad) = NaN;
plot(r2_vals, r1_vals, '--','Color',colors(6,:), 'LineWidth',2, ...
     'DisplayName',sprintf('r_1 = r_2/R_1 (%.4f)',R1));

xlabel('r_2'); ylabel('r(k)');
title('r_k(r_2) for k=1:5 and line r_1 = r_2/R_1');
yline(0,'--','Color',[0.5 0.5 0.5]);
legend('Location','northwest');
grid on;
hold off;
