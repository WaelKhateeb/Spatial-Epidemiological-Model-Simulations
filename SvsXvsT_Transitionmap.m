function solve_SI_pde()
    %------------------------------------------------------------------
    % User-defined plotting interval (edit these before running)
    %------------------------------------------------------------------
    t_min_plot = 1000;         
    t_max_plot = 5000;         

    %------------------------------------------------------------------
    % 1. Model Parameters
    %------------------------------------------------------------------
    beta = 0.0094;   d = 0.01;  mu0 = 0.1;  mu1 = 10;
    b = 0.03;        A = 1;
    r1 = 0.07;       r2 = 0.01;
    L  = 5*pi;       Nx = 200;

    %------------------------------------------------------------------
    % 2. Endemic equilibrium (largest I*)
    %------------------------------------------------------------------
    a2 = beta*(d+mu0);
    a1 = d*(d+mu0) + beta*b*(d+mu1) - beta*A;
    a0 = d*b*(d+mu1) - beta*A*b;
    I_roots = roots([a2,a1,a0]);
    I_pos   = I_roots(imag(I_roots)==0 & real(I_roots)>eps);
    if isempty(I_pos), error('No positive endemic equilibrium found.'); end
    I_star = max(I_pos);
    S_star = A/(d + beta*I_star);
    fprintf('Selected equilibrium: S^* = %.6f, I^* = %.6f\n', S_star, I_star);

    %------------------------------------------------------------------
    % Jacobian & eigenvalues (unchanged)
    %------------------------------------------------------------------
    fprintf('\nEigenvalues of the reaction Jacobian (k = 0):\n');
    hI_star = mu0 + (mu1-mu0)*(b/(I_star+b));
    dhI_dI  = - (mu1-mu0)*b/(I_star+b)^2;
    J_react = [ -d - beta*I_star,   -beta*S_star;
                 beta*I_star,      beta*S_star - d - (dhI_dI*I_star + hI_star) ];
    ev0 = eig(J_react);
    fprintf('  λ₁ = %.6f %+.6fi,  λ₂ = %.6f %+.6fi\n', ...
            real(ev0(1)), imag(ev0(1)), real(ev0(2)), imag(ev0(2)));
    fprintf('\nEigenvalues of the RD model including diffusion:\n');
    fprintf(' Mode   (nπ/L)^2        λ₁                λ₂\n');
    fprintf('---------------------------------------------------------\n');
    for n = 0:(Nx-1)
        k2 = (n*pi/L)^2;
        M = J_react - diag([r1*k2, r2*k2]);
        ev = eig(M);
        fprintf('  %3d    %10.4e    % .6f%+.6fi    % .6f%+.6fi\n', ...
                n, k2, real(ev(1)), imag(ev(1)), real(ev(2)), imag(ev(2)));
    end

    %------------------------------------------------------------------
    % 3. Spatial grid & initial perturbation
    %------------------------------------------------------------------
    x  = linspace(0, L, Nx);
    dx = x(2)-x(1);
    perturb_scale = 0.01;
    cos_perturb   = cos(0.4*x);
    S0 = S_star + perturb_scale * cos_perturb;
    I0 = I_star + perturb_scale * cos_perturb;
    U0 = [S0, I0]';

    %------------------------------------------------------------------
    % 4. Time integration setup
    %------------------------------------------------------------------
    tspan = [0 15000];
    nt    = 5000;
    t_eval = linspace(tspan(1),tspan(2),nt);
    opts = odeset('RelTol',1e-14,'AbsTol',1e-14);
    sol  = ode15s(@(t,U) pde_rhs(t,U,Nx,dx,beta,d,mu0,mu1,b,A,r1,r2),tspan,U0,opts);
    U    = deval(sol,t_eval);
    S_sol = U(1:Nx,:);
    I_sol = U(Nx+1:end,:);

    %------------------------------------------------------------------
    % 5. Time-window indices
    %------------------------------------------------------------------
    idx = find(t_eval>=t_min_plot & t_eval<=t_max_plot);
    if isempty(idx), error('No points in [%g, %g].',t_min_plot,t_max_plot); end
    t_plot = t_eval(idx);

    %------------------------------------------------------------------
    % 6. Prompt for variable choice and plot
    %------------------------------------------------------------------
    varChoice = input('Enter variable to plot (S or I): ','s');
    switch upper(varChoice)
        case 'S'
            % get 2-colors from rainbow
            cmap2 = jet(2);

            % 3D surface with rainbow
            figure
            surf(x, t_plot, S_sol(:,idx)','EdgeColor','none')
            xlabel('x'), ylabel('t'), zlabel('S(x,t)')
            title(sprintf('S(x,t), t ∈ [%.1f, %.1f]',t_min_plot,t_max_plot))
            shading interp
            colormap(jet)
            colorbar
            view(40,30)

            % 2D heat-map with rainbow
            figure
            imagesc(x, t_plot, S_sol(:,idx)')
            set(gca,'YDir','normal')
            xlabel('x'), ylabel('t')
            title('Heat-map of S(x,t)')
            colormap(jet)
            colorbar

            % initial condition check (rainbow lines)
            figure
            plot(x, S_sol(:,1), 'Color',cmap2(1,:), 'LineWidth',1.5, ...
                 x, S_star*ones(size(x)), 'Color',cmap2(2,:), 'LineStyle','--','LineWidth',1.5)
            xlabel('x'), ylabel('S(x,0)')
            legend('Perturbed S','S^*','Location','best')
            title('Initial S vs. Equilibrium')
            grid on

        case 'I'
            cmap2 = jet(2);

            % 3D surface with rainbow
            figure
            surf(x, t_plot, I_sol(:,idx)','EdgeColor','none')
            xlabel('x'), ylabel('t'), zlabel('I(x,t)')
            title(sprintf('I(x,t), t ∈ [%.1f, %.1f]',t_min_plot,t_max_plot))
            shading interp
            colormap(jet)
            colorbar
            view(40,30)

            % 2D heat-map with rainbow
            figure
            imagesc(x, t_plot, I_sol(:,idx)')
            set(gca,'YDir','normal')
            xlabel('x'), ylabel('t')
            title('Heat-map of I(x,t)')
            colormap(jet)
            colorbar

            % initial condition check (rainbow lines)
            figure
            plot(x, I_sol(:,1), 'Color',cmap2(1,:), 'LineWidth',1.5, ...
                 x, I_star*ones(size(x)), 'Color',cmap2(2,:), 'LineStyle','--','LineWidth',1.5)
            xlabel('x'), ylabel('I(x,0)')
            legend('Perturbed I','I^*','Location','best')
            title('Initial I vs. Equilibrium')
            grid on

        otherwise
            error('Invalid choice. Enter S or I.');
    end
end

% -------------------------------------------------------------------------
function dUdt = pde_rhs(~,U,Nx,dx,beta,d,mu0,mu1,b,A,r1,r2)
    S = U(1:Nx);  I = U(Nx+1:end);
    S_xx = laplacian_neumann(S,dx);
    I_xx = laplacian_neumann(I,dx);
    hI   = mu0 + (mu1-mu0).*(b./(I+b));
    dSdt = A - d.*S - beta.*S.*I + r1*S_xx;
    dIdt = beta.*S.*I - d.*I - hI.*I + r2*I_xx;
    dUdt = [dSdt; dIdt];
end

% -------------------------------------------------------------------------
function u_xx = laplacian_neumann(u,dx)
    N = numel(u);
    u_xx = zeros(size(u));
    u_xx(2:N-1) = (u(1:N-2) - 2*u(2:N-1) + u(3:N)) / dx^2;
    u_xx(1)     = (u(2)     - 2*u(1)     + u(1))   / dx^2;
    u_xx(N)     = (u(N-1)   - 2*u(N)     + u(N))   / dx^2;
end
