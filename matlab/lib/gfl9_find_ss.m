function x_ss = gfl9_find_ss(X_net, S_ref, p)
%GFL9_FIND_SS  Find the co-rotating steady state of the GFL9 system.
%
%  Solves gfl9_ode(x, p, X_net, S_ref) = 0 via Newton-Raphson with a
%  finite-difference Jacobian.  Does NOT require the Optimization Toolbox.
%
%  Inputs
%    X_net  – N×N reactance matrix
%    S_ref  – N×1 signed reference powers (>0 CBR injection, <0 SE absorption)
%    p      – control params from gfl_control_params()
%
%  Output
%    x_ss   – 9N×1 steady-state state vector

N   = numel(S_ref);
nSt = 9*N;

% ── Algebraic initial guess (quasi-static, iq ≈ 0) ───────────────────────
id_ss0 = S_ref(:);
Iq_g0  = zeros(N, 1);
Id_g0  = id_ss0;
Vd_g0  = 1.0 - X_net * Iq_g0;
Vq_g0  =       X_net * Id_g0;
phi0   = atan2(Vq_g0, Vd_g0);

x0 = zeros(nSt, 1);
for k = 1:N
    b      = (k-1)*9;
    phi_k  = phi0(k);
    x0(b+1) = phi_k;
    x0(b+3) = id_ss0(k) / p.Ki_p;             % eps_p:  id = Ki_p * eps_p at SS
    x0(b+4) = p.Rf * id_ss0(k) / p.Ki_i;      % eps_id: vd residual
    x0(b+6) = id_ss0(k);                       % id
    x0(b+8) =  Vd_g0(k)*cos(phi_k) + Vq_g0(k)*sin(phi_k);  % Vd_ff = Vd_l
    x0(b+9) = -Vd_g0(k)*sin(phi_k) + Vq_g0(k)*cos(phi_k);  % Vq_ff = Vq_l
end

% ── Newton-Raphson (no Optimization Toolbox required) ────────────────────
f_ode = @(x) gfl9_ode(x, p, X_net, S_ref, 1.0);
x_ss  = newton_raphson(f_ode, x0, 1e-9, 80);

res = max(abs(f_ode(x_ss)));
if res > 1e-6
    % Fallback: run ODE to reduce residual, then Newton from there
    warning('gfl9_find_ss: Newton residual %.1e > 1e-6. Trying ODE+Newton.', res);
    ode_f    = @(t, x) gfl9_ode(x, p, X_net, S_ref, 1.0);
    opts_ode = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'MaxStep', 1e-2);
    [~, xt]  = ode15s(ode_f, [0 1], x0, opts_ode);   % short run: damp fast transients only
    x_ss     = newton_raphson(f_ode, xt(end,:)', 1e-9, 80);
    res2 = max(abs(f_ode(x_ss)));
    if res2 > 1e-5
        warning('gfl9_find_ss: ODE+Newton residual still %.1e. Using ODE endpoint.', res2);
        x_ss = xt(end,:)';
    end
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

    % Build numerical Jacobian (n evaluations of f)
    J = zeros(n);
    for j = 1:n
        xp     = x;
        xp(j)  = xp(j) + eps_j;
        J(:,j) = (f(xp) - fval) / eps_j;
    end

    % Newton step with backtracking line search
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
