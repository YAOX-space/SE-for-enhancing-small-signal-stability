function tbl = compute_first_principles_eigenvalues(gscr_values, varargin)
%COMPUTE_FIRST_PRINCIPLES_EIGENVALUES  Compute dominant eigenvalue at each
%  target gSCR using the 9-state GFL ODE (with GFF feedforward filter).
%
%  Does NOT use paper Table V.  For each target gSCR the 39-node network
%  is scaled so that lambda_max(diag(S)*Z) = 1/gSCR, the co-rotating
%  steady state is found via fsolve, the 9-state Jacobian is linearised
%  numerically, and the eigenvalue nearest 88 rad/s is returned.
%
%  Usage
%    tbl = compute_first_principles_eigenvalues()
%        Sweeps gSCR = [2.650, 2.926, 3.256, 3.666] (Table V check)
%
%    tbl = compute_first_principles_eigenvalues(gscr_values)
%        User-specified row vector of gSCR targets.
%
%    tbl = compute_first_principles_eigenvalues(gscr_values, 'Plot', true)
%        Also plots sigma vs gSCR and compares with Table V.
%
%  Output  tbl  – table with columns: gSCR, sigma, omega, damping
%
%  Dependency: build_network_from_edges, max_positive_eig (lib/)

rootDir = fileparts(mfilename('fullpath'));
addpath(fullfile(rootDir, 'lib'));
dataDir = fullfile(rootDir, 'data');

% ── parse options ────────────────────────────────────────────────────────
doPlot = false;
if nargin >= 2
    for k = 1:2:numel(varargin)
        if strcmpi(varargin{k}, 'Plot'), doPlot = varargin{k+1}; end
    end
end

if nargin < 1 || isempty(gscr_values)
    gscr_values = [2.650, 2.926, 3.256, 3.666];   % Table V check points
end
gscr_values = gscr_values(:)';

% ── control parameters (Table IX) ────────────────────────────────────────
p.omega0  = 2*pi*50;
p.Lf      = 0.05 / p.omega0;
p.Rf      = 0.002 / p.omega0;
p.Kp_pll  = 26;    p.Ki_pll = 7800;
p.Kp_p    = 0.5;   p.Ki_p   = 5;
p.Kp_i    = 1;     p.Ki_i   = 10;
p.tau_ff  = 0.01;  % GFF filter time constant (Table IX)

% ── IEEE 39-node network, calibrated so gSCR = 2.650 with 9 CBRs ─────────
lines39 = readtable(fullfile(dataDir, 'table_xi_ieee39_lines.csv'));
nbus39  = max([lines39.from; lines39.to]);
[~, Z39] = build_network_from_edges(lines39, nbus39, 36);
keep39   = setdiff(1:nbus39, 36);
b2r39    = zeros(1, nbus39);  b2r39(keep39) = 1:38;
cbr_rows = b2r39([30 31 32 33 34 35 37 38 39]);
Z_raw    = Z39(cbr_rows, cbr_rows);   % 9×9 CBR submatrix
S_CBR    = ones(9, 1);

% Calibration: scale so lambda_max(diag(S)*Z) = 1/2.650
k_base = (1/2.650) / max_positive_eig(diag(S_CBR) * Z_raw);
Z_base = k_base * Z_raw;   % gSCR = 2.650 at this scale
fprintf('Baseline gSCR check: %.4f  (should be 2.650)\n', ...
    1/max_positive_eig(diag(S_CBR)*Z_base));

% ── paper Table V for comparison ─────────────────────────────────────────
tbl_v_gscr  = [2.650, 2.926, 3.256, 3.666];
tbl_v_sigma = [0.090, -1.758, -3.485, -5.104];
tbl_v_omega = [88.342, 88.949, 89.346, 89.564];

% ── sweep ────────────────────────────────────────────────────────────────
N   = numel(S_CBR);
nG  = numel(gscr_values);
sigma_out   = zeros(1, nG);
omega_out   = zeros(1, nG);
gscr_actual = zeros(1, nG);

for gi = 1:nG
    g_tgt = gscr_values(gi);

    % Scale network: gSCR = gSCR_base / alpha  →  alpha = gSCR_base / gSCR_target
    alpha  = 2.650 / g_tgt;
    X_net  = alpha * Z_base;
    gscr_actual(gi) = 1/max_positive_eig(diag(S_CBR)*X_net);

    fprintf('\n[gSCR target=%.4f  actual=%.4f]\n', g_tgt, gscr_actual(gi));

    % ── steady state via fsolve ───────────────────────────────────────────
    x_ss = find_steady_state_9state(X_net, S_CBR, p);
    res  = gfl9_rotating_ode(x_ss, p, X_net, S_CBR);
    fprintf('  Max |residual| at SS: %.2e\n', max(abs(res)));

    % ── numerical Jacobian ───────────────────────────────────────────────
    nSt  = 9*N;
    f0   = gfl9_rotating_ode(x_ss, p, X_net, S_CBR);
    J    = zeros(nSt);
    eps_fd = 1e-7;
    for j = 1:nSt
        xp = x_ss;  xp(j) = xp(j) + eps_fd;
        J(:,j) = (gfl9_rotating_ode(xp, p, X_net, S_CBR) - f0) / eps_fd;
    end

    ev = eig(J);
    [~, bi] = min(abs(abs(imag(ev)) - 88));
    lam = ev(bi);
    sigma_out(gi) = real(lam);
    omega_out(gi) = abs(imag(lam));

    fprintf('  First-principles:  sigma = %+.4f   omega = %.4f\n', ...
        sigma_out(gi), omega_out(gi));

    % Compare with Table V if this gSCR is in the table
    [~, ti] = min(abs(tbl_v_gscr - g_tgt));
    if abs(tbl_v_gscr(ti) - g_tgt) < 0.01
        fprintf('  Paper Table V:     sigma = %+.4f   omega = %.4f\n', ...
            tbl_v_sigma(ti), tbl_v_omega(ti));
        fprintf('  Gap:  %.4f\n', sigma_out(gi) - tbl_v_sigma(ti));
    end
end

% ── build output table ───────────────────────────────────────────────────
damping = -sigma_out ./ sqrt(sigma_out.^2 + omega_out.^2);
tbl = table(gscr_values(:), sigma_out(:), omega_out(:), damping(:), ...
    'VariableNames', {'gSCR','sigma','omega','damping'});

fprintf('\n═══ First-Principles Eigenvalue Table ═══\n');
fprintf('%-8s  %-10s  %-10s  %-10s\n', 'gSCR','sigma','omega','damping');
for gi = 1:nG
    fprintf('%-8.4f  %+10.4f  %10.4f  %+10.4f\n', ...
        gscr_values(gi), sigma_out(gi), omega_out(gi), damping(gi));
end

if doPlot
    figure('Color','w','Name','First-Principles vs Table V','Position',[100 100 700 350]);
    subplot(1,2,1);
    plot(gscr_values, sigma_out, 'b-o', 'LineWidth',1.6, 'DisplayName','First-principles'); hold on;
    plot(tbl_v_gscr, tbl_v_sigma, 'r--s', 'LineWidth',1.4, 'DisplayName','Paper Table V');
    yline(0,'--k');
    xlabel('gSCR');  ylabel('\sigma (real part of dominant eigenvalue)');
    title('\sigma vs gSCR');  legend;  grid on;

    subplot(1,2,2);
    plot(gscr_values, omega_out, 'b-o', 'LineWidth',1.6, 'DisplayName','First-principles'); hold on;
    plot(tbl_v_gscr, tbl_v_omega, 'r--s', 'LineWidth',1.4, 'DisplayName','Paper Table V');
    xlabel('gSCR');  ylabel('\omega (rad/s)');
    title('\omega vs gSCR');  legend;  grid on;

    resultsDir = fullfile(rootDir, 'results');
    if ~exist(resultsDir,'dir'), mkdir(resultsDir); end
    exportgraphics(gcf, fullfile(resultsDir, 'first_principles_vs_tableV.png'), ...
        'Resolution',180,'BackgroundColor','white');
    fprintf('\n  → results/first_principles_vs_tableV.png\n');
end
end


% ═══════════════════════════════════════════════════════════════════════
%  find_steady_state_9state  –  find co-rotating steady state via fsolve
%
%  Works for both stable and unstable operating points because it solves
%  dxdt(x_ss) = 0 directly (Newton's method), rather than integrating
%  forward in time which diverges when the steady state is unstable.
% ═══════════════════════════════════════════════════════════════════════
function x_ss = find_steady_state_9state(X_net, S_CBR, p)
N   = numel(S_CBR);
nSt = 9*N;

% Algebraic initial guess (quasi-static, iq=0)
id_ss0 = S_CBR(:);
Iq_g0  = zeros(N,1);
Id_g0  = id_ss0;
Vd_g0  = 1.0 - X_net * Iq_g0;
Vq_g0  =       X_net * Id_g0;
phi0   = atan2(Vq_g0, Vd_g0);

x0 = zeros(nSt,1);
for k = 1:N
    b = (k-1)*9;
    phi_k = phi0(k);
    x0(b+1) = phi_k;
    % eps_pll (b+2) = 0: Vq_l=0 at SS → dphi=0 requires eps_pll=0
    x0(b+3) = id_ss0(k) / p.Ki_p;            % eps_p:  id_ref=Ki_p*eps_p at SS
    x0(b+4) = p.Rf * id_ss0(k) / p.Ki_i;     % eps_id: vd_vsc residual
    % eps_iq (b+5) = 0: iq_ref=0 and iq=0
    x0(b+6) = id_ss0(k);                      % id
    % iq (b+7) = 0
    x0(b+8) = Vd_g0(k)*cos(phi_k) + Vq_g0(k)*sin(phi_k);   % Vd_ff = Vd_l
    x0(b+9) = -Vd_g0(k)*sin(phi_k) + Vq_g0(k)*cos(phi_k);  % Vq_ff = Vq_l
end

% Primary: fsolve finds dxdt=0 directly (stable or unstable fixed point)
f    = @(x) gfl9_rotating_ode(x, p, X_net, S_CBR);
opts = optimoptions('fsolve','Display','off','TolFun',1e-10,'TolX',1e-10,...
    'MaxIterations',500,'MaxFunctionEvaluations',50000);
[x_ss, fval, exitflag] = fsolve(f, x0, opts);
fprintf('    fsolve: exitflag=%d, max|f|=%.2e\n', exitflag, max(abs(fval)));

if exitflag <= 0 || max(abs(fval)) > 1e-8
    % Fallback: run ODE to reduce residual then fsolve from there
    warning('fsolve did not converge (exitflag=%d), trying ODE+fsolve ...', exitflag);
    ode_f    = @(t,x) gfl9_rotating_ode(x, p, X_net, S_CBR);
    opts_ode = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',5e-3);
    [~, xt]  = ode15s(ode_f, [0 15], x0, opts_ode);
    [x_ss, fval2, ef2] = fsolve(f, xt(end,:)', opts);
    fprintf('    ODE+fsolve: exitflag=%d, max|f|=%.2e\n', ef2, max(abs(fval2)));
end
end


% ═══════════════════════════════════════════════════════════════════════
%  gfl9_rotating_ode  –  9-state GFL in co-rotating frame (with GFF)
%
%  States per converter k:
%    [phi, eps_pll, eps_p, eps_id, eps_iq, id, iq, Vd_ff, Vq_ff]
%
%  phi = theta_k - omega0*t  (angle relative to global rotating frame)
%  Vd_ff, Vq_ff = GFF-filtered terminal voltages (local frame)
% ═══════════════════════════════════════════════════════════════════════
function dxdt = gfl9_rotating_ode(x, p, X_net, S_CBR)
N   = numel(S_CBR);
w0  = p.omega0;  Lf = p.Lf;  Rf = p.Rf;
tau = p.tau_ff;

phi     = x(1:9:end);
eps_pll = x(2:9:end);
eps_p   = x(3:9:end);
eps_id  = x(4:9:end);
eps_iq  = x(5:9:end);
id_loc  = x(6:9:end);
iq_loc  = x(7:9:end);
Vd_ff   = x(8:9:end);
Vq_ff   = x(9:9:end);

% Global-frame currents
Id_g = id_loc.*cos(phi) - iq_loc.*sin(phi);
Iq_g = id_loc.*sin(phi) + iq_loc.*cos(phi);

% Network voltages (global frame, quasi-static purely-inductive network)
Vd_g = 1.0 - X_net * Iq_g;
Vq_g =       X_net * Id_g;

% Terminal voltage in converter local dq frame
Vd_l =  Vd_g.*cos(phi) + Vq_g.*sin(phi);
Vq_l = -Vd_g.*sin(phi) + Vq_g.*cos(phi);

% Measured active power
P_meas = Vd_l.*id_loc + Vq_l.*iq_loc;
Pref   = S_CBR(:);

% Outer power control loop → id reference
id_ref = p.Kp_p*(Pref - P_meas) + p.Ki_p*eps_p;
iq_ref = zeros(N,1);

% Inner current control using GFF-filtered voltages as feedforward
%   GFF 1/(1+tau*s) at 88 rad/s: |G|=0.752, phase=-41° → uncancelled
%   fraction ~0.66, which is the key coupling that creates instability.
vd_vsc = p.Kp_i*(id_ref - id_loc) + p.Ki_i*eps_id - w0*Lf*iq_loc + Vd_ff;
vq_vsc = p.Kp_i*(iq_ref - iq_loc) + p.Ki_i*eps_iq + w0*Lf*id_loc + Vq_ff;

% PLL (co-rotating frame: dphi/dt = dtheta/dt - w0)
dphi     = p.Kp_pll*Vq_l + p.Ki_pll*eps_pll;
deps_pll = Vq_l;
deps_p   = Pref - P_meas;
deps_id  = id_ref - id_loc;
deps_iq  = iq_ref - iq_loc;

% LCL filter (simplified: only Lf, no Cf state)
did = (vd_vsc - Vd_l - Rf*id_loc)./Lf + w0*iq_loc;
diq = (vq_vsc - Vq_l - Rf*iq_loc)./Lf - w0*id_loc;

% GFF filter dynamics
dVd_ff = (Vd_l - Vd_ff) / tau;
dVq_ff = (Vq_l - Vq_ff) / tau;

dxdt = zeros(9*N, 1);
dxdt(1:9:end) = dphi;
dxdt(2:9:end) = deps_pll;
dxdt(3:9:end) = deps_p;
dxdt(4:9:end) = deps_id;
dxdt(5:9:end) = deps_iq;
dxdt(6:9:end) = did;
dxdt(7:9:end) = diq;
dxdt(8:9:end) = dVd_ff;
dxdt(9:9:end) = dVq_ff;
end
