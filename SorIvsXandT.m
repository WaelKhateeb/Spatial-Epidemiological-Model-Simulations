function solve_SI_pde()
    % Parameters
    beta = 0.1;     % Transmission rate
    d    = 0.1;     % Natural death rate
    mu0  = 0.1;     % Min disease‐induced mortality
    mu1  = 0.2;     % Max disease‐induced mortality
    b    = 0.3;     % Half‐saturation constant
    A    = 0.4;     % Recruitment rate
    r1   = 21.7222; % Susceptible diffusion (fixed to match pdepe run)
    r2   = 0.09897; % Infected diffusion     % Diffusion coefficient for I
    L = 5*pi;       % Domain length [0, L]
    Nx = 600;       % Number of spatial grid points

    % Find equilibrium (I_star, S_star)
    fun = @(I) beta.*(A/(d + beta.*I)).*I - d.*I - I.*(mu0 + (mu1 - mu0).*(b./(I + b)));
    I_star = fzero(fun, 0.75); % Numerically find I_star
    S_star = A / (d + beta*I_star);
    fprintf('Equilibrium: S* = %.4f, I* = %.4f\n', S_star, I_star);

    % Spatial grid
    x = linspace(0, L, Nx);
    dx = x(2) - x(1);

    % Initial conditions with small perturbation (zero-mean cosine)
    perturb_scale = 0.01;
    cos_perturb = cos(0.4*x);
    S0 = S_star + perturb_scale*cos_perturb;
    I0 = I_star + perturb_scale*cos_perturb;
    U0 = [S0, I0]'; % Combined as a column vector (first Nx for S, next Nx for I)

    % Time parameters
    tspan = [0 10000000000]; 
    nt = 10000;
    t_eval = linspace(tspan(1), tspan(2), nt);

    % Solve the PDE system
    opts = odeset('RelTol',1e-5, 'AbsTol',1e-5);
    sol = ode15s(@(t,U) pde_rhs(t, U, Nx, dx, beta, d, mu0, mu1, b, A, r1, r2), tspan, U0, opts);
    U = deval(sol, t_eval);
    S_sol = U(1:Nx, :);
    I_sol = U(Nx+1:end, :);

    % Let the user choose which variable to plot
    varChoice = input('Enter variable to plot (S or I): ', 's');
    if strcmpi(varChoice, 'S')
        % 3D surface plot of S(x,t)
        figure
        surf(x, t_eval, S_sol', 'EdgeColor', 'none')
        xlabel('Position (x)'), ylabel('Time (t)'), zlabel('S(x,t)')
        title('Susceptible Population Dynamics')
        colorbar, shading interp, colormap parula, view(40,30)

        % Equilibrium verification at t=0 for S
        figure
        plot(x, S_sol(:,1), 'b', x, S_star*ones(size(x)), '--r', 'LineWidth', 1.5)
        xlabel('Position (x)'), ylabel('S(x,0)')
        legend('Initial S', 'Equilibrium S^*', 'Location', 'best')
        title('Initial Condition Verification for S')
        grid on
    elseif strcmpi(varChoice, 'I')
        % 3D surface plot of I(x,t)
        figure
        surf(x, t_eval, I_sol', 'EdgeColor', 'none')
        xlabel('Position (x)'), ylabel('Time (t)'), zlabel('I(x,t)')
        title('Infected Population Dynamics')
        colorbar, shading interp, colormap parula, view(40,30)

        % Equilibrium verification at t=0 for I
        figure
        plot(x, I_sol(:,1), 'b', x, I_star*ones(size(x)), '--r', 'LineWidth', 1.5)
        xlabel('Position (x)'), ylabel('I(x,0)')
        legend('Initial I', 'Equilibrium I^*', 'Location', 'best')
        title('Initial Condition Verification for I')
        grid on
    else
        error('Invalid choice. Please enter S or I.');
    end
end

function dUdt = pde_rhs(~, U, Nx, dx, beta, d, mu0, mu1, b, A, r1, r2)
    % Split U into S and I components
    S = U(1:Nx);
    I = U(Nx+1:end);

    % Compute Laplacians with Neumann BCs
    S_xx = laplacian_neumann(S, dx);
    I_xx = laplacian_neumann(I, dx);

    % Mortality function for I
    hI = mu0 + (mu1 - mu0).*(b./(I + b));

    % PDE equations
    dSdt = A - d.*S - beta.*S.*I + r1*S_xx;
    dIdt = beta.*S.*I - d.*I - hI.*I + r2*I_xx;
    dUdt = [dSdt; dIdt];
end

function u_xx = laplacian_neumann(u, dx)
    % Second derivative with Neumann BCs (du/dx=0 at boundaries)
    N = length(u);
    u_xx = zeros(size(u));
    
    % Interior points
    u_xx(2:N-1) = (u(1:N-2) - 2*u(2:N-1) + u(3:N)) / dx^2;
    
    % Boundary points (using ghost nodes)
    u_xx(1) = (u(2) - 2*u(1) + u(1)) / dx^2;      % Left boundary
    u_xx(N) = (u(N-1) - 2*u(N) + u(N)) / dx^2;      % Right boundary
end
