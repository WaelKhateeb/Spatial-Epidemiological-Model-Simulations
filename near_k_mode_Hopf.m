function solve_SI_pde()
    %------------------------------------------------------------------
    % User-defined plotting interval (edit these before running)
    %------------------------------------------------------------------
    t_min_plot = 00000;         % minimum t to include in surface plot
    t_max_plot = 20000;% maximum t to include in surface plot

    %------------------------------------------------------------------
    % 1. Model Parameters
    %------------------------------------------------------------------
    beta = 0.0088;           % Transmission rate
    d    = 0.01;             % Natural death rate
    mu0  = 0.1;              % Min disease-induced mortality
    mu1  = 10;               % Max disease-induced mortality
    b    = 0.03;             % Half-saturation constant
    A    = 1;                % Recruitment rate
    r1   = 0.05;% Diffusion coefficient for S
    r2   = 0.01;             % Diffusion coefficient for I
    L    = 5*pi;             % Domain length [0, L]
    Nx   = 200;              % Number of spatial points

    %------------------------------------------------------------------
    % 2. Endemic equilibrium (largest I*)
    %------------------------------------------------------------------
    a2 = beta * (d + mu0);
    a1 = d*(d + mu0) + beta*b*(d + mu1) - beta*A;
    a0 = d*b*(d + mu1) - beta*A*b;
    I_roots = roots([a2, a1, a0]);
    I_pos   = I_roots(imag(I_roots)==0 & real(I_roots)>eps);
    if isempty(I_pos)
        error('No positive endemic equilibrium found.');
    end
    I_star = max(I_pos);
    S_star = A / (d + beta*I_star);
    fprintf('Selected equilibrium: S^* = %.6f, I^* = %.6f\n', S_star, I_star);

    %------------------------------------------------------------------
    % 2a. Reaction‐term Jacobian at (S*, I*)
    %------------------------------------------------------------------
    % Define h(I) and its derivative h'(I) at I = I_star
    hI_star = mu0 + (mu1 - mu0) * (b / (I_star + b));
    dhI_dI  = - (mu1 - mu0) * b / (I_star + b)^2;

    % Compute Jacobian entries of the reaction ODE at equilibrium
    J11 = -d - beta * I_star;
    J12 = -beta * S_star;
    J21 =  beta * I_star;
    J22 = beta * S_star - d - (dhI_dI * I_star + hI_star);

    J_react = [J11, J12;
               J21, J22];

    %------------------------------------------------------------------
    % 2b. Print eigenvalues of the full RD linearization (all modes)
    %------------------------------------------------------------------
    fprintf('\nEigenvalues of the reaction Jacobian (k = 0):\n');
    ev0 = eig(J_react);
    fprintf('  λ₁ = %.6f %+.6fi,  λ₂ = %.6f %+.6fi\n', ...
            real(ev0(1)), imag(ev0(1)), real(ev0(2)), imag(ev0(2)));

    fprintf('\nEigenvalues of the RD model including diffusion:\n');
    fprintf(' Mode   (nπ/L)^2        λ₁ (with diff)                λ₂ (with diff)\n');
    fprintf('-----------------------------------------------------------------------\n');
    for n = 0:(Nx-1)
        k2 = (n * pi / L)^2; 
        % Construct mode‐dependent Jacobian: J_react - diag(r1*k2, r2*k2)
        M = J_react - diag([r1 * k2, r2 * k2]);
        ev = eig(M);
        fprintf('  %3d    %10.4e    % .6f%+.6fi    % .6f%+.6fi\n', ...
                n, k2, ...
                real(ev(1)), imag(ev(1)), ...
                real(ev(2)), imag(ev(2)));
    end

    %------------------------------------------------------------------
    % 3. Spatial grid & initial perturbation
    %------------------------------------------------------------------
    x  = linspace(0, L, Nx);
    dx = x(2) - x(1);
    perturb_scale = 0.00001;
    cos_perturb   = cos(0.2 * x);
    S0 = S_star + perturb_scale * cos_perturb;
    I0 = I_star + perturb_scale * cos_perturb;
    U0 = [S0, I0]';

    %------------------------------------------------------------------
    % 4. Time integration setup
    %------------------------------------------------------------------
    tspan  = [0 20000];
    nt     = 5000;
    t_eval = linspace(tspan(1), tspan(2), nt);

    opts = odeset('RelTol',1e-8, 'AbsTol',1e-8);
    sol  = ode15s(@(t,U) pde_rhs(t, U, Nx, dx, beta, d, mu0, mu1, b, A, r1, r2), ...
                  tspan, U0, opts);
    U    = deval(sol, t_eval);
    S_sol = U(1:Nx, :);
    I_sol = U(Nx+1:end, :);

    %------------------------------------------------------------------
    % 5. Determine indices for the user-defined time window
    %------------------------------------------------------------------
    idx = find(t_eval >= t_min_plot & t_eval <= t_max_plot);
    if isempty(idx)
        error('No solution points in the specified plot interval [%g, %g].', t_min_plot, t_max_plot);
    end
    t_plot = t_eval(idx);

    %------------------------------------------------------------------
    % 6. Prompt for variable choice and plot
    %------------------------------------------------------------------
    varChoice = input('Enter variable to plot (S or I): ', 's');
    switch upper(varChoice)
        case 'S'
            figure
            surf(x, t_plot, S_sol(:, idx)', 'EdgeColor', 'none')
            xlabel('Position (x)'), ylabel('Time (t)'), zlabel('S(x,t)')
            title(sprintf('Susceptible Dynamics (t \\in [%.1f, %.1f])', t_min_plot, t_max_plot))
            colorbar, shading interp, colormap parula, view(40,30)

            figure
            plot(x, S_sol(:,1), 'b', x, S_star*ones(size(x)), '--r', 'LineWidth', 1.5)
            xlabel('Position (x)'), ylabel('S(x,0)')
            legend('Initial S (perturbed)', 'Equilibrium S^*', 'Location', 'best')
            title('Initial Condition Verification for S')
            grid on

        case 'I'
            figure
            surf(x, t_plot, I_sol(:, idx)', 'EdgeColor', 'none')
            xlabel('Position (x)'), ylabel('Time (t)'), zlabel('I(x,t)')
            title(sprintf('Infected Dynamics (t \\in [%.1f, %.1f])', t_min_plot, t_max_plot))
            colorbar, shading interp, colormap parula, view(40,30)

            figure
            plot(x, I_sol(:,1), 'b', x, I_star*ones(size(x)), '--r', 'LineWidth', 1.5)
            xlabel('Position (x)'), ylabel('I(x,0)')
            legend('Initial I (perturbed)', 'Equilibrium I^*', 'Location', 'best')
            title('Initial Condition Verification for I')
            grid on

        otherwise
            error('Invalid choice. Please enter S or I.');
    end
end

% -------------------------------------------------------------------------
function dUdt = pde_rhs(~, U, Nx, dx, beta, d, mu0, mu1, b, A, r1, r2)
    S = U(1:Nx);
    I = U(Nx+1:end);

    S_xx = laplacian_neumann(S, dx);
    I_xx = laplacian_neumann(I, dx);
    hI   = mu0 + (mu1 - mu0) .* (b ./ (I + b));

    dSdt = A - d .* S - beta .* S .* I + r1 * S_xx;
    dIdt = beta .* S .* I    - d .* I    - hI .* I    + r2 * I_xx;

    dUdt = [dSdt; dIdt];
end

% -------------------------------------------------------------------------
function u_xx = laplacian_neumann(u, dx)
    N    = numel(u);
    u_xx = zeros(size(u));
    u_xx(2:N-1) = (u(1:N-2) - 2*u(2:N-1) + u(3:N)) / dx^2;
    u_xx(1)     = (u(2)     - 2*u(1)     + u(1))   / dx^2;  % Neumann at x=0
    u_xx(N)     = (u(N-1)   - 2*u(N)     + u(N))   / dx^2;  % Neumann at x=L
end
