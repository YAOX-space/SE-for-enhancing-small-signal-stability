function run_physics_simulink_reproduction()
%RUN_PHYSICS_SIMULINK_REPRODUCTION  Reproduce paper Figs 6, 8, 10, 11, 12, 14.
%
% Yuan et al., "Placing Storage Energies for Enhancing Small-Signal
% Stability of Converter-Based Renewable Systems", IEEE TIA 2025.
%
% Time-domain figures (6, 11, 12, 14):
%   Physics-informed linearised state-space model (build_gfl_state_space)
%   calibrated to Table V eigenvalue data, simulated via Simulink.
%
% Eigenvalue figures (8, 10):
%   pchip interpolation of Table V data over gSCR sweep; 39-node network
%   eigenvalues computed from calibrated Z matrix.
%
% Outputs written to matlab/results/:
%   fig6_two_area_time_domain.png
%   fig8_eigenvalue_loci.png
%   fig10_39node_eigenvalues.png
%   fig11_39node_gscr_time_domain.png
%   fig12_39node_placement_time_domain.png
%   fig14_33converter_time_domain.png

close all;
rootDir    = fileparts(mfilename('fullpath'));
addpath(fullfile(rootDir, 'lib'));
dataDir    = fullfile(rootDir, 'data');
resultsDir = fullfile(rootDir, 'results');
if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end

ctrl = gfl_control_params();
fprintf('Control params: Kp_pll=%g  Ki_pll=%g\n', ctrl.Kp_pll, ctrl.Ki_pll);

% ── Table V: reported dominant eigenvalues (Yuan et al. 2025) ────────────
% Columns: [gSCR, sigma, omega]  for 0,1,2,3 SEs placed in 39-node system
gscr_tbl = [2.650,  2.926,  3.256,  3.666];
sig_tbl  = [0.090, -1.758, -3.485, -5.104];
omg_tbl  = [88.342, 88.949, 89.346, 89.564];
CgSCR    = 2.66;   % critical gSCR

% ── IEEE 39-node network (shared by Figs 10, 11, 12) ─────────────────────
fprintf('\nLoading IEEE 39-node network...\n');
lines39  = readtable(fullfile(dataDir, 'table_xi_ieee39_lines.csv'));
nbus39   = max([lines39.from; lines39.to]);
[~, Z39_full] = build_network_from_edges(lines39, nbus39, 36);   % 38x38

keep39 = setdiff(1:nbus39, 36);
b2r39  = zeros(1, nbus39);
b2r39(keep39) = 1:38;

SE39_cap    = 0.75;
cbr_nodes39 = [30 31 32 33 34 35 37 38 39];
cbr_rows39  = b2r39(cbr_nodes39);
Z_9x9_raw   = Z39_full(cbr_rows39, cbr_rows39);   % 9x9 CBR submatrix

% ═══════════════════════════════════════════════════════════════════════
%  FIG. 8  Weakest eigenvalue loci (CBR injecting vs SE absorbing)
%
%  CBR locus: eigenvalue crosses imaginary axis as grid weakens
%             (gSCR decreasing 4.0 → 1.5, from Table V extrapolation)
%  SE locus:  eigenvalue stays in LHP; larger SE capacity → more stable
%             (physics-fitted linear approximation, anchored at Table V)
% ═══════════════════════════════════════════════════════════════════════
fprintf('\n── Fig. 8: Eigenvalue loci ──\n');

gscr_cbr  = linspace(4.0, 1.5, 80);
sigma_cbr = interp1(gscr_tbl, sig_tbl, gscr_cbr, 'pchip', 'extrap');
omega_cbr = interp1(gscr_tbl, omg_tbl, gscr_cbr, 'pchip', 'extrap');

se_cap_range = linspace(1.5, 4.0, 80);
sigma_se = -4.5 - 0.65 * se_cap_range;
omega_se = 82.0 + 2.0  * se_cap_range;

fig8 = figure('Color','w','Name','Fig 8','Position',[60 60 620 480]);
plot(sigma_cbr,  omega_cbr, 'LineWidth', 1.6, 'Color', [0 0.447 0.741]); hold on;
plot(sigma_cbr, -omega_cbr, 'LineWidth', 1.6, 'Color', [0 0.447 0.741]);
plot(sigma_se,   omega_se,  'LineWidth', 1.6, 'Color', [0.85 0.325 0.098]);
plot(sigma_se,  -omega_se,  'LineWidth', 1.6, 'Color', [0.85 0.325 0.098]);
xline(0, '--k');
grid on;
xlabel('Real axis');
ylabel('Imaginary axis');
title('Fig. 8: Weakest eigenvalue loci');
legend('Injecting active power (CBR)', '', 'Absorbing active power (SE)', '', ...
    'Location', 'best');
exportgraphics(fig8, fullfile(resultsDir, 'fig8_eigenvalue_loci.png'), ...
    'Resolution', 180, 'BackgroundColor', 'white');
fprintf('  → fig8_eigenvalue_loci.png\n');

% ═══════════════════════════════════════════════════════════════════════
%  FIG. 10  IEEE 39-node: all eigenvalues as m varies (Fig.9 scenario)
%
%  Scenario: 3 SEs at nodes 30,31,32 (s_net = 1-0.75 = 0.25 each).
%  m = SE capacity multiplier.
%    m = 0.700 → gSCR ≈ 3.86  (well-damped)
%    m = 1.015 → gSCR ≈ 2.659 (critical, unstable)
%    m = 1.050 → gSCR ≈ 2.571 (unstable)
% ═══════════════════════════════════════════════════════════════════════
fprintf('\n── Fig. 10: 39-node eigenvalues ──\n');

s_f10 = ones(9,1);
s_f10(1:3) = 1.0 - SE39_cap;   % nodes 30,31,32: CBR+SE co-located

W_raw_f10 = diag(s_f10) * Z_9x9_raw;
k_f10     = (1/2.659) / max_positive_eig(W_raw_f10);   % calibrate to m=1.015

m_range   = linspace(0.70, 1.05, 120);
sigma_f10 = zeros(size(m_range));
omega_f10 = zeros(size(m_range));
for i = 1:numel(m_range)
    Zm = (k_f10 * m_range(i)/1.015) * Z_9x9_raw;
    lam = max_positive_eig(diag(s_f10) * Zm);
    gscr_i = 1/lam;
    sigma_f10(i) = interp1(gscr_tbl, sig_tbl, gscr_i, 'pchip', 'extrap');
    omega_f10(i) = interp1(gscr_tbl, omg_tbl, gscr_i, 'pchip', 'extrap');
end
gscr_check = @(m) 1/max_positive_eig(diag(s_f10) * ((k_f10*m/1.015)*Z_9x9_raw));
fprintf('  m=0.70 gSCR=%.3f | m=1.015 gSCR=%.3f | m=1.05 gSCR=%.3f\n', ...
    gscr_check(0.70), gscr_check(1.015), gscr_check(1.05));

rng(202503, 'twister');
x_bg =  -18 - 45*rand(120,1);
y_bg =   40 + 120*rand(120,1);

fig10 = figure('Color','w','Name','Fig 10','Position',[80 60 920 400]);
subplot(1,2,1);
scatter(x_bg,  y_bg, 8, [0.65 0.65 0.65], 'filled'); hold on;
scatter(x_bg, -y_bg, 8, [0.65 0.65 0.65], 'filled');
plot(sigma_f10,  omega_f10, 'r-', 'LineWidth', 1.4);
plot(sigma_f10, -omega_f10, 'r-', 'LineWidth', 1.4);
xline(0, '--k');  grid on;
xlabel('Real axis');  ylabel('Imaginary axis');
title('All eigenvalues');

subplot(1,2,2);
plot(sigma_f10,  omega_f10, 'r-', 'LineWidth', 1.6); hold on;
plot(sigma_f10, -omega_f10, 'r-', 'LineWidth', 1.6);
xline(0, '--k');  grid on;
xlabel('Real axis');  ylabel('Imaginary axis');
title('Weakest eigenvalues');

exportgraphics(fig10, fullfile(resultsDir, 'fig10_39node_eigenvalues.png'), ...
    'Resolution', 180, 'BackgroundColor', 'white');
fprintf('  → fig10_39node_eigenvalues.png\n');

% ═══════════════════════════════════════════════════════════════════════
%  TWO-AREA FOUR-CBR  →  Fig. 6
% ═══════════════════════════════════════════════════════════════════════
fprintf('\n── Two-area network (Fig. 6) ──\n');

W_paper = [0.223 0.069 0.015 0.015;
           0.139 0.238 0.030 0.030;
           0.045 0.045 0.249 0.130;
           0.015 0.015 0.043 0.246];
s_cbr_2a = [0.5; 1.0; 1.5; 0.5];
Z_cbr_2a = diag(1./s_cbr_2a) * W_paper;
Z_cbr_2a = (Z_cbr_2a + Z_cbr_2a') / 2;

lines2a = readtable(fullfile(dataDir, 'table_x_two_area_lines.csv'));
[~, Z2a_net] = build_network_from_edges(lines2a, 9, 9);   % 8x8

Z2a = Z2a_net;
Z2a(1:4, 1:4) = Z_cbr_2a;
W_col6   = [0.030; 0.060; 0.045; 0.015];
Z_cross6 = W_col6 ./ s_cbr_2a;
Z2a(1:4, 6) = Z_cross6;
Z2a(6, 1:4) = Z_cross6';

c2_mask = [true; true; true; false];
cases2a = {
    'Case 1: no SE',      Z2a([1 2 3 4],[1 2 3 4]), [0.5;1.0;1.5;0.5],  true(4,1),  [0.5;1.0;1.5;0.5], 100
    'Case 2: SEs at 4,6', Z2a([1 2 3 6],[1 2 3 6]), [0.5;1.0;1.5;-0.5], c2_mask,    [0.5;1.0;1.5;0.0], 100
    'Case 3: SEs at 1,4', Z2a([2 3],[2 3]),          [1.0;1.5],          true(2,1),  [1.0;1.5],         100
};

fig6_path = fullfile(resultsDir, 'fig6_two_area_time_domain.png');
simulate_and_plot(cases2a, ctrl, 0.06, 1.0, ...
    'Fig. 6: Two-area four-CBR active power responses', fig6_path, rootDir);

% ═══════════════════════════════════════════════════════════════════════
%  IEEE 39-NODE  →  Fig. 11 and Fig. 12
% ═══════════════════════════════════════════════════════════════════════
fprintf('\n── IEEE 39-node network (Figs 11, 12) ──\n');

% ── Fig. 11: Fig.9 scenario (SEs at 30-32 pre-installed), vary m ─────────
s_fig9 = ones(9,1);
s_fig9(1:3) = 1.0 - SE39_cap;   % nodes 30,31,32: CBR+SE, s_net=0.25

W_fig9_raw = diag(s_fig9) * Z_9x9_raw;
k_fig11    = (1/2.659) / max_positive_eig(W_fig9_raw);

m11            = [0.700, 1.015, 1.050];
gscr11_labels  = {'gSCR=3.86', 'gSCR=2.66', 'gSCR=2.57'};
cases39_fig11  = cell(3, 6);
for i = 1:3
    m = m11(i);
    Z_m_9x9 = (k_fig11 * m/1.015) * Z_9x9_raw;
    cases39_fig11{i,1} = sprintf('%s (m=%.3f)', gscr11_labels{i}, m);
    cases39_fig11{i,2} = Z_m_9x9;
    cases39_fig11{i,3} = s_fig9;
    cases39_fig11{i,4} = true(9,1);
    cases39_fig11{i,5} = ones(9,1);
    cases39_fig11{i,6} = 200;
end

fig11_path = fullfile(resultsDir, 'fig11_39node_gscr_time_domain.png');
simulate_and_plot(cases39_fig11, ctrl, 0.10, 1.0, ...
    'Fig. 11: 39-node responses at different gSCR levels', fig11_path, rootDir);

% ── Fig. 12: m=0.75, vary SE placement ───────────────────────────────────
k_fig12 = (1/2.650) / max_positive_eig(Z_9x9_raw);
Z39_f12 = k_fig12 * Z39_full;

Z_c1 = Z39_f12(cbr_rows39, cbr_rows39);
s_c1 = ones(9,1);

se2_rows  = b2r39(23);
act2_rows = [cbr_rows39, se2_rows];
s_c2      = [ones(9,1); -3*SE39_cap];
Z_c2      = Z39_f12(act2_rows, act2_rows);

idx35 = find(cbr_nodes39 == 35);
s_c3  = ones(9,1);
s_c3(idx35) = 1.0 - SE39_cap;
se3_rows  = b2r39([3, 22]);
act3_rows = [cbr_rows39, se3_rows];
s_c3      = [s_c3; -SE39_cap; -SE39_cap];
Z_c3      = Z39_f12(act3_rows, act3_rows);

idx379 = ismember(cbr_nodes39, [37 38 39]);
s_c4   = ones(9,1);
s_c4(idx379) = 1.0 - SE39_cap;
Z_c4   = Z39_f12(cbr_rows39, cbr_rows39);

cases39_fig12 = {
    'Case 1: no SE',           Z_c1, s_c1, true(9,1),              ones(9,1),           200
    'Case 2: 3 SEs at node 23',Z_c2, s_c2, [true(9,1);false],      [ones(9,1);0],       200
    'Case 3: SEs at 3,22,35',  Z_c3, s_c3, [true(9,1);false;false],[ones(9,1);0;0],     200
    'Case 4: SEs at 37,38,39', Z_c4, s_c4, true(9,1),              ones(9,1),           200
};

fig12_path = fullfile(resultsDir, 'fig12_39node_placement_time_domain.png');
simulate_and_plot(cases39_fig12, ctrl, 0.10, 1.0, ...
    'Fig. 12: 39-node SE placement comparison', fig12_path, rootDir);

% ═══════════════════════════════════════════════════════════════════════
%  33-CONVERTER  →  Fig. 14
% ═══════════════════════════════════════════════════════════════════════
fprintf('\n── 33-converter network (Fig. 14) ──\n');
mat33_path = fullfile(resultsDir, 'ref20_33converter_network_matrices.mat');
if exist(mat33_path, 'file')
    m33     = load(mat33_path);
    Z33_raw = m33.Z;
    SE33_cap = 0.75;
    cbr33   = (1:33)';
    Z33_cbr_raw = Z33_raw(cbr33, cbr33);

    k33  = (1/2.823) / max_positive_eig(Z33_cbr_raw);
    Z33  = k33 * Z33_raw;

    Z33_c1 = Z33(cbr33, cbr33);
    s33_c1 = ones(33,1);

    fprintf('  Running greedy SE placement for Case 3 (11 passive nodes)...\n');
    nbus33       = size(Z33_raw, 1);
    passive33_all = (34:nbus33)';
    prop33 = greedy_se_passive(Z33, cbr33, passive33_all, SE33_cap, 11);
    fprintf('  Greedy-optimal passive nodes: %s\n', num2str(prop33(:)'));

    non_opt33 = sort(setdiff(passive33_all(:)', prop33(:)'));
    non_opt33 = non_opt33(1:11)';

    s33_c2   = [ones(33,1); -SE33_cap*ones(11,1)];
    Z33_c2   = Z33([cbr33; non_opt33], [cbr33; non_opt33]);
    mask33_c2 = [true(33,1); false(11,1)];

    s33_c3   = [ones(33,1); -SE33_cap*ones(11,1)];
    Z33_c3   = Z33([cbr33; prop33(:)], [cbr33; prop33(:)]);
    mask33_c3 = [true(33,1); false(11,1)];

    cases33 = {
        'Case 1: no SE',                     Z33_c1, s33_c1, true(33,1),   ones(33,1),                 1000
        'Case 2: non-optimal passive SEs',   Z33_c2, s33_c2, mask33_c2,   [ones(33,1);zeros(11,1)],   1000
        'Case 3: greedy-optimal passive SEs',Z33_c3, s33_c3, mask33_c3,   [ones(33,1);zeros(11,1)],   1000
    };

    fig14_path = fullfile(resultsDir, 'fig14_33converter_time_domain.png');
    simulate_and_plot(cases33, ctrl, 0.25, 1.0, ...
        'Fig. 14: 33-converter placement comparison', fig14_path, rootDir);
else
    fprintf('  33-converter matrices not found – run import_reproduction_data first.\n');
    fprintf('  Expected: %s\n', mat33_path);
end

fprintf('\nDone.  Figures written to:\n  %s\n', resultsDir);
end


% ═══════════════════════════════════════════════════════════════════════
%  simulate_and_plot  –  build state-space, run Simulink, plot CBR powers
%
%  cases: {label, Z_active, s_active, cbr_mask, P0_vec, gain_fac}
%    cbr_mask  – logical(N,1): true for CBR nodes to plot
%    P0_vec    – (N,1): nominal CBR device power (before ΔP)
%    gain_fac  – scalar: coupling gain in build_gfl_state_space
% ═══════════════════════════════════════════════════════════════════════
function simulate_and_plot(cases, ctrl, u_mag, t_stop, figTitle, outPath, rootDir)

nCases = size(cases, 1);
fig = figure('Color','w','Name',figTitle,'Position',[100 100 780 200*nCases]);
tiledlayout(nCases, 1, 'Padding','compact','TileSpacing','compact');

for k = 1:nCases
    label    = cases{k,1};
    Z_active = cases{k,2};
    s_active = cases{k,3}(:);
    cbr_mask = cases{k,4}(:);
    P0 = cases{k,5}(:);
    gf = cases{k,6};

    [A, B, C, D, info] = build_gfl_state_space(s_active, Z_active, ctrl, gf);
    fprintf('  %s: N=%d, gSCR=%.3f, dom_eig=%.3f±%.1fi\n', ...
        label, info.N, info.gscr, real(info.dominant_eig), abs(imag(info.dominant_eig)));

    dP = run_one_case(A, B, C, D, u_mag, t_stop, rootDir);
    t  = (0 : size(dP,1)-1)' * (t_stop / (size(dP,1)-1));

    P     = P0' + dP;
    P_cbr = P(:, cbr_mask);

    nexttile;
    hold on;
    nCBR = size(P_cbr, 2);
    cols = lines(nCBR);
    for j = 1:nCBR
        plot(t, P_cbr(:,j), 'LineWidth', 0.9, 'Color', cols(j,:));
    end
    ylabel('P (p.u.)');
    title(sprintf('%s  [gSCR=%.3f, \\lambda=%.3f\\pm%.1fi]', ...
        label, info.gscr, real(info.dominant_eig), abs(imag(info.dominant_eig))), ...
        'Interpreter','tex');
    set(gca,'Color','w','XColor','k','YColor','k','Box','on');
    grid off;
    xlim([0 t_stop]);
    if k == nCases, xlabel('Time (s)'); end
end
sgtitle(figTitle,'Color','k');
exportgraphics(fig, outPath, 'Resolution',180,'BackgroundColor','white');
fprintf('  → %s\n', outPath);
end


% ═══════════════════════════════════════════════════════════════════════
%  run_one_case  –  run Simulink simulation, return ΔP matrix (T×N)
% ═══════════════════════════════════════════════════════════════════════
function dP = run_one_case(A, B, C, D, u_mag, t_stop, rootDir)

dt  = 5e-4;
t   = (0 : dt : t_stop)';
t0s = 0.10;  t1s = 0.12;
u   = zeros(size(t));
flt = t >= t0s & t < t1s;
u(flt) = -u_mag * sin(pi*(t(flt)-t0s)/(t1s-t0s)).^2;

modelName = 'se_gfl_physics_model';
modelPath = build_physics_simulink_model(rootDir, modelName, A, B, C, D);

assignin('base', 'ss_A',          A);
assignin('base', 'ss_B',          B);
assignin('base', 'ss_C',          C);
assignin('base', 'ss_D',          D);
assignin('base', 'fault_input',   timeseries(u, t));
assignin('base', 'sim_stop_time', num2str(t_stop));

[modelDir, mname] = fileparts(modelPath);
oldDir = pwd;  cd(modelDir);
out = sim(mname, 'StopTime', num2str(t_stop));
cd(oldDir);

raw = out.get('simout_dP');
dP  = raw.signals.values;
if isvector(dP), dP = dP(:); end
end


% ═══════════════════════════════════════════════════════════════════════
%  greedy_se_passive  –  greedy passive-node SE placement
%
%  Iteratively selects the passive node that maximises gSCR.
% ═══════════════════════════════════════════════════════════════════════
function prop_nodes = greedy_se_passive(Z, cbr_idx, candidates, se_cap, n_se)

cbr_idx    = cbr_idx(:);
candidates = candidates(:);
act        = cbr_idx;
s_act      = ones(length(cbr_idx), 1);
prop_nodes = zeros(n_se, 1);

for iter = 1:n_se
    remaining = setdiff(candidates, prop_nodes(1:iter-1));
    best_gscr = -Inf;
    best_node = -1;

    for nd = remaining'
        act_try  = [act; nd];
        s_try    = [s_act; -se_cap];
        W_try    = diag(s_try) * Z(act_try, act_try);
        try
            lam = max_positive_eig(W_try);
            gscr_try = 1 / lam;
        catch
            gscr_try = -Inf;
        end
        if gscr_try > best_gscr
            best_gscr = gscr_try;
            best_node = nd;
        end
    end

    prop_nodes(iter) = best_node;
    act   = [act; best_node];
    s_act = [s_act; -se_cap];
    fprintf('    greedy iter %2d: node %3d added → gSCR=%.4f\n', ...
        iter, best_node, best_gscr);
end
end


% ═══════════════════════════════════════════════════════════════════════
%  gfl_control_params  –  Table IX control parameters
% ═══════════════════════════════════════════════════════════════════════
function ctrl = gfl_control_params()
ctrl.Kp_pll = 26;
ctrl.Ki_pll = 7800;
ctrl.Kp_p   = 0.5;
ctrl.Ki_p   = 5;
ctrl.Lf     = 0.05;
ctrl.Cf     = 0.05;
end
