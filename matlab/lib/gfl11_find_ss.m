function x_ss = gfl11_find_ss(X_net, S_ref, p)
%GFL11_FIND_SS  Find the co-rotating steady state of the GFL11 system.
%
%  Solves gfl11_ode(x, p, X_net, S_ref) = 0 via Newton-Raphson with a
%  finite-difference Jacobian.  Does NOT require the Optimization Toolbox.
%
%  Inputs
%    X_net  – N×N reactance matrix
%    S_ref  – N×1 signed reference powers (>0 CBR, <0 SE)
%    p      – control params from gfl_control_params()
%
%  Output
%    x_ss   – 11N×1 steady-state state vector

N   = numel(S_ref);
nSt = 11*N;
w0  = p.omega0;

% ── Algebraic initial guess ───────────────────────────────────────────────
% At steady state with vCq=0 (PLL lock):
%   dvCd=0 → id = i_gd   (d-axis Cf balance)
%   dvCq=0 → iq = i_gq + w0*Cf*vCd  (q-axis Cf balance; Cf absorbs reactive)
%
% Cf (per-unit physical: Cf_pu/omega0) provides reactive current at omega0:
%   i_cap = w0 * Cf = omega0 * (0.05/omega0) = 0.05 p.u.
% The capacitor supplies this Q, so grid q-current = -i_cap (leading).

iq_cap = w0 * p.Cf;             % = 0.05 p.u. reactive current from Cf
id_ss0 = S_ref(:);              % d-axis grid current ≈ active power (vCd≈1)
Iq_g0  = -iq_cap * ones(N,1);   % grid q-current: Cf compensates
Id_g0  = id_ss0;

% Capacitor bus voltage in global frame
Vd_g0  = 1.0 - X_net * Iq_g0;  % d-axis: boosted by Cf reactive current
Vq_g0  =       X_net * Id_g0;   % q-axis: from active power flow
phi0   = atan2(Vq_g0, Vd_g0);   % PLL locks vCq=0

x0 = zeros(nSt, 1);
for k = 1:N
    b     = (k-1)*11;
    ph    = phi0(k);
    % Capacitor voltage in local dq frame (rotate global vC by -phi)
    vCd_k =  Vd_g0(k)*cos(ph) + Vq_g0(k)*sin(ph);

    x0(b+1)  = ph;
    % eps_pll = 0 at SS (vCq=0 → integrator at rest)
    x0(b+3)  = id_ss0(k) / p.Ki_p;            % eps_p
    x0(b+4)  = p.Rf * id_ss0(k) / p.Ki_i;     % eps_id
    x0(b+5)  = (p.Kp_i + p.Rf) * iq_cap / p.Ki_i;  % eps_iq
    x0(b+6)  = id_ss0(k);                     % id = i_gd at SS
    x0(b+7)  = iq_cap;                        % iq = iq_cap (Cf current)
    x0(b+8)  = vCd_k;                         % Vd_ff → vCd at SS
    % Vq_ff = 0 at SS (vCq=0)
    x0(b+10) = vCd_k;                         % vCd
    % vCq = 0 at SS
end

% ── Newton-Raphson (no Optimization Toolbox required) ─────────────────────
f_ode = @(x) gfl11_ode(x, p, X_net, S_ref, 1.0);
x_ss  = newton_raphson(f_ode, x0, 1e-9, 80);

res = max(abs(f_ode(x_ss)));
if res > 1e-6
    warning('gfl11_find_ss: Newton residual %.1e > 1e-6. Trying short ODE+Newton.', res);
    ode_f    = @(t, x) gfl11_ode(x, p, X_net, S_ref, 1.0);
    % Short integration (0.05 s) so unstable systems fail fast, not hang
    opts_ode = odeset('RelTol', 1e-5, 'AbsTol', 1e-4, 'MaxStep', 5e-4);
    try
        [~, xt] = ode15s(ode_f, [0 0.05], x0, opts_ode);
        x_ss    = newton_raphson(f_ode, xt(end,:)', 1e-9, 40);
        res2    = max(abs(f_ode(x_ss)));
        if res2 > 1e-5
            warning('gfl11_find_ss: ODE+Newton residual still %.1e.', res2);
        end
    catch
        warning('gfl11_find_ss: ODE fallback failed. Using best Newton result.');
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
