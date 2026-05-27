function x_ss = gfl13_find_ss(X_net, S_ref, p)
%GFL13_FIND_SS  Find the co-rotating steady state of the GFL13 system.
%
%  Builds an initial guess from the 11-state steady state (gfl11_find_ss)
%  and augments with I_gd_g, I_gq_g using the algebraic quasi-static
%  relations (which are exact at SS).  Then runs Newton-Raphson on
%  gfl13_ode to converge to the true 13-state fixed point.
%
%  Inputs
%    X_net  – N×N reactance matrix (X_net = omega0 * L_net)
%    S_ref  – N×1 signed reference powers (>0 CBR, <0 SE)
%    p      – control params from gfl_control_params()
%
%  Output
%    x_ss   – 13N×1 steady-state state vector

N   = numel(S_ref);
nSt = 13*N;

% ── Get 11-state SS as initial guess ─────────────────────────────────────
x11 = gfl11_find_ss(X_net, S_ref, p);

% ── Augment with algebraic grid-side currents ─────────────────────────────
%  At SS: I_gd_g = X_net\vCq_g,  I_gq_g = X_net\(V_inf - vCd_g)
phi  = x11(1:11:end);
vCd  = x11(10:11:end);
vCq  = x11(11:11:end);

vCd_g    = vCd.*cos(phi) - vCq.*sin(phi);
vCq_g    = vCd.*sin(phi) + vCq.*cos(phi);
I_gd_g0  = X_net \ vCq_g;
I_gq_g0  = X_net \ (1.0*ones(N,1) - vCd_g);

x0 = zeros(nSt, 1);
for k = 1:N
    b11 = (k-1)*11;
    b13 = (k-1)*13;
    x0(b13+1:b13+11) = x11(b11+1:b11+11);
    x0(b13+12) = I_gd_g0(k);
    x0(b13+13) = I_gq_g0(k);
end

% ── Newton-Raphson on gfl13_ode ───────────────────────────────────────────
f_ode = @(x) gfl13_ode(x, p, X_net, S_ref, 1.0);
x_ss  = newton_raphson(f_ode, x0, 1e-9, 80);

res = max(abs(f_ode(x_ss)));
if res > 1e-6
    warning('gfl13_find_ss: Newton residual %.1e > 1e-6. Trying short ODE+Newton.', res);
    ode_f    = @(t, x) gfl13_ode(x, p, X_net, S_ref, 1.0);
    opts_ode = odeset('RelTol', 1e-5, 'AbsTol', 1e-4, 'MaxStep', 5e-4);
    try
        [~, xt] = ode15s(ode_f, [0 0.05], x0, opts_ode);
        x_ss    = newton_raphson(f_ode, xt(end,:)', 1e-9, 40);
        res2    = max(abs(f_ode(x_ss)));
        if res2 > 1e-5
            warning('gfl13_find_ss: ODE+Newton residual still %.1e.', res2);
        end
    catch
        warning('gfl13_find_ss: ODE fallback failed. Using best Newton result.');
    end
end

% Final residual check — throw so callers can catch and return NaN
res_final = max(abs(f_ode(x_ss)));
if res_final > 0.1
    error('gfl13_find_ss: SS not found (residual=%.2e). System may be strongly unstable.', res_final);
end
end


% ── Newton-Raphson with finite-difference Jacobian ────────────────────────
function x = newton_raphson(f, x0, tol, max_iter)
x     = x0(:);
n     = numel(x);
eps_j = 1e-7;

for iter = 1:max_iter
    fval = f(x);
    if max(abs(fval)) <= tol, return; end

    J = zeros(n);
    for j = 1:n
        xp     = x;
        xp(j)  = xp(j) + eps_j;
        J(:,j) = (f(xp) - fval) / eps_j;
    end

    dx = -J \ fval;
    step = 1.0;
    f0n  = norm(fval);
    for ls = 1:8
        xnew = x + step * dx;
        if norm(f(xnew)) < f0n, break; end
        step = step * 0.5;
    end
    x = x + step * dx;
end
end
