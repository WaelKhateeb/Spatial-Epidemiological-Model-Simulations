clear;
close all;

ddebiftoolpath = '../../';
addpath(fullfile(ddebiftoolpath, 'ddebiftool'), ...
    fullfile(ddebiftoolpath, 'ddebiftool_extra_psol'), ...
    fullfile(ddebiftoolpath, 'ddebiftool_extra_nmfm'), ...
    fullfile(ddebiftoolpath, 'ddebiftool_utilities'));

format compact;
set(groot, 'defaultTextInterpreter', 'LaTeX');

parnames = {'k', 'tau', 'delay'};
cind = [parnames; num2cell(1:length(parnames))];
ind = struct(cind{:});

bounds = { ...
    'max_bound', [ind.k 150; ind.delay 22], ...
    'max_step', [1, 2], ...
    'min_bound', [ind.k 0.5; ind.delay 0]};

funcs = set_symfuncs(@sym_congestionControlModel, 'sys_tau', @() ind.tau);

k = 7;
c = 1;
d = 0.1;
a = 5;
m = 1;
r = 15;
tau = 1;
delay = 12.083;

initial_w = d*a/(c*exp(-d*delay)*m - d);
initial_q = r*c*exp(-d*delay)*a*(k*c*m*exp(-d*delay) - k*d - d*a) ...
    /(k*(c*exp(-d*delay)*m - d)^2);

stst = dde_stst_create('x', [initial_w; initial_q]);
stst.parameter(ind.k) = k;
stst.parameter(ind.tau) = tau;
stst.parameter(ind.delay) = delay;

method = df_mthod(funcs, 'stst');
[stst, success] = p_correc(funcs, stst, [], [], method.point);

method_stst = df_mthod(funcs, 'stst');
method_stst.stability.minimal_real_part = -30;
stst.stability = p_stabil(funcs, stst, method_stst.stability);

figure(1); clf;
plot(stst.stability.l1, '*');
title('Stability plot of steady state');
xlabel('$\Re(\lambda)$');
ylabel('$\Im(\lambda)$');

contpar = ind.delay;
steadystate_br = SetupStst(funcs, ...
    'x', stst.x, ...
    'parameter', stst.parameter, ...
    'step', 0.000001, ...
    'contpar', contpar, ...
    'max_step', [contpar, 0.5], ...
    bounds{:});

figure(2); clf;
ax2 = gca;
n_steps = 50;
steadystate_br = br_contn(funcs, steadystate_br, n_steps, 'plotaxis', ax2);
xlabel('$\tau$');
ylabel('$x(1)$');

[steadystate_br, ~, ind_hopf, bif1types] = LocateSpecialPoints(funcs, steadystate_br);
nunst_eqs = GetStability(steadystate_br);
fprintf('Hopf bifurcation near point %d\n', ind_hopf);

stst.parameter(ind.delay) = 12;
stst.x = [11; 74];
method_stst = df_mthod(funcs, 'stst');
stst.stability = p_stabil(funcs, stst, method_stst.stability);

hopf = p_tohopf(funcs, stst);
hopf = nmfm_hopf(funcs, hopf);
hopf.nmfm

disp('Branch off at Hopf bifurcation');
fprintf('Initial correction of periodic orbits at Hopf:\n');

[per1, suc] = SetupPsol(funcs, steadystate_br, ind_hopf, ...
    'print_residual_info', 1, ...
    'intervals', 50, ...
    'degree', 7, ...
    'max_bound', [contpar, 35], ...
    'max_step', [contpar, 0.2], ...
    'matrix', 'full');

figure(3); clf;
ax3 = gca;
xlabel('$\tau$');
ylabel('Amplitude');
per1 = br_contn(funcs, per1, 200, 'plotaxis', ax3);

norms = zeros(length(per1.point), 1);
for i = 1:length(per1.point)
    point = per1.point(i);
    stability = p_stabil(funcs, point, method_stst.stability);
    norms(i) = max(abs(stability.mu));
end

tau_values = zeros(length(per1.point), 1);
for i = 1:length(per1.point)
    tau_values(i) = per1.point(i).parameter(ind.delay);
end

figure(4); clf;
plot(tau_values, norms, '-r');
xlabel('$\tau$');
ylabel('Maximum Floquet multiplier magnitude');
title('Periodic-orbit stability along the Hopf branch');
