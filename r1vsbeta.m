%% Complete MATLAB script: r₁^{(k)} vs. β with overlays and direct β–solving
clear; clc;

%% 1) Model parameters
d        = 0.01;
mu0      = 0.1;
mu1_val  = 10;
b        = 0.03;
A_val    = 1;
l_val    = 5;    %

%% 2) Domain / scan settings
r2_fixed  = 0.01;
beta_vals = linspace(0.001,0.013,10000);
Kmax      = 50;

%% 3) Preallocate
r1_vals   = nan(numel(beta_vals),Kmax);
r1_cond   = nan(numel(beta_vals),Kmax);

%% 4) First, directly solve for β such that Σ₁+Σ₄=0 AND endemic‐equilibrium holds
syms beta_sym I_sym real
% Endemic equilibrium:  β A/(d+βI) − d − [μ₀ + (μ₁−μ₀)b/(I+b)] = 0
eq1 = beta_sym*A_val/(d + beta_sym*I_sym) ...
      - d ...
      - (mu0 + (mu1_val - mu0)*b/(I_sym + b)) == 0;
% Σ₁+Σ₄ = −βI − d + I*(μ₁−μ₀)b/(I+b)^2 = 0
eq2 = -beta_sym*I_sym - d ...
      + I_sym*(mu1_val - mu0)*b/(I_sym + b)^2 == 0;

% Use vpasolve with an initial guess near mid β and moderate I
init_beta = mean(beta_vals);
init_I    = 1;  
sol       = vpasolve([eq1, eq2], [beta_sym, I_sym], [init_beta, init_I]);

% Extract real, positive β‐roots
beta_zero = double(sol.beta_sym);
beta_zero = beta_zero(imag(beta_zero)==0 & beta_zero>0);

%% 5) Scan and compute r₁^{(k)} + cond‐curves exactly as before
syms I
for i = 1:numel(beta_vals)
    beta = beta_vals(i);
    % --- endemic I
    eq = beta*A_val/(d + beta*I) ...
         - d ...
         - (mu0 + (mu1_val - mu0)*b/(I + b)) == 0;
    Ir = double(solve(eq, I));
    valid = (imag(Ir)==0) & (Ir>eps);
    if ~any(valid), continue; end
    I2 = max(Ir(valid));

    % --- Sigmas
    h       = @(I) mu0 + (mu1_val - mu0)*b./(I + b);
    h_p     = @(I) - (mu1_val - mu0)*b ./ (I + b).^2;
    S1      = -beta*I2 - d;
    S2      = -h(I2)   - d;
    S3      =  beta*I2;
    S4      = -I2*h_p(I2);

    % --- r₁ and cond
    for k = 1:Kmax
        numr1       = l_val^4 * (S1*S4 - S2*S3 - S1*r2_fixed*(k^2/l_val^2));
        denr1       = k^2     * (S4*l_val^2 - r2_fixed*k^2);
        r1_vals(i,k)= numr1./denr1;
        r1_cond(i,k)= (l_val^2/(k^2))*(S1 + S4) - r2_fixed;
    end
end

%% 6) Plot everything
figure('Position',[100,100,900,550]); hold on;
cols = lines(Kmax);

% (a) r₁^{(k)} curves
for k = 1:Kmax
    y = r1_vals(:,k);
    y(y<=0 | ~isfinite(y)) = NaN;
    plot(beta_vals, y, 'Color',cols(k,:), 'LineWidth',1.5, ...
         'DisplayName',sprintf('r_1^{(%d)}(β)',k));
end

% (b) condition curves
for k = 1:Kmax
    y = r1_cond(:,k);
    y(y<=0 | ~isfinite(y)) = NaN;
    plot(beta_vals, y, '--', 'Color',cols(k,:), 'LineWidth',1.5, ...
         'DisplayName',sprintf('cond, k=%d',k));
end

% (c) vertical line(s) at Σ₁+Σ₄=0
for bz = beta_zero(:)'
    xline(bz, '--k', 'LineWidth',2, 'DisplayName','Σ₁+Σ₄=0');
end

xlabel('\beta','FontSize',12);
ylabel('r_1','FontSize',12);
title(sprintf('r_1^{(k)} vs. \\beta  (r_2=%.2f)',r2_fixed),'FontSize',14);
legend('Location','northeast');
grid on; hold off;
