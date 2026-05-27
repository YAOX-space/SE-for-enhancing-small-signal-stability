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
p.Cf      = 0.05  / p.omega0;  % LCL capacitor (11-state model)
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

    % ── steady state (Newton-Raphson, no Optimization Toolbox) ───────────
    x_ss = gfl13_find_ss(X_net, S_CBR, p);
    res  = gfl13_ode(x_ss, p, X_net, S_CBR, 1.0);
    fprintf('  Max |residual| at SS: %.2e\n', max(abs(res)));

    % ── numerical Jacobian ───────────────────────────────────────────────
    nSt  = 13*N;
    f0   = gfl13_ode(x_ss, p, X_net, S_CBR, 1.0);
    J    = zeros(nSt);
    eps_fd = 1e-7;
    for j = 1:nSt
        xp = x_ss;  xp(j) = xp(j) + eps_fd;
        J(:,j) = (gfl13_ode(xp, p, X_net, S_CBR, 1.0) - f0) / eps_fd;
    end

    ev = eig(J);
    % Select dominant PLL mode: max real part among |omega| in [70,130] rad/s
    % (avoids spurious fast LCL modes and selects the most dangerous slow mode)
    mask = abs(imag(ev)) >= 70 & abs(imag(ev)) <= 130;
    ev_pll = ev(mask);
    if isempty(ev_pll)
        [~, bi] = min(abs(abs(imag(ev)) - 88));
        lam = ev(bi);
    else
        [~, bi] = max(real(ev_pll));
        lam = ev_pll(bi);
    end
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


% Local fsolve-based functions removed.
% Steady-state computation now uses gfl9_find_ss (lib/), which requires no
% Optimization Toolbox.  ODE evaluations use gfl9_ode (lib/) with V_inf=1.0.
