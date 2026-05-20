clear;

ddebiftoolpath = '../../';
addpath(fullfile(ddebiftoolpath, 'ddebiftool'), ...
    fullfile(ddebiftoolpath, 'ddebiftool_extra_symbolic'));

if dde_isoctave()
    pkg load symbolic;
end

parnames = {'k', 'tau', 'delay'};
syms(parnames{:});
par = cell2sym(parnames);

syms w wtau q qtau;

dw_dt = delay*15*w*(1 - w/k) - delay*w*q/(5 + w);
dq_dt = delay*exp(-0.1*delay)*wtau*qtau/(5 + wtau) - delay*0.1*q;

[fstr, derivs] = dde_sym2funcs( ...
    [dw_dt; dq_dt], ...
    [w, wtau; q, qtau], ...
    par, ...
    'filename', 'sym_congestionControlModel');
