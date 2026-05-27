function run_physics_simulink_reproduction()
%RUN_PHYSICS_SIMULINK_REPRODUCTION  Reproduce paper Figs 6, 8, 10, 11, 12, 14.
%
% Yuan et al., "Placing Storage Energies for Enhancing Small-Signal
% Stability of Converter-Based Renewable Systems", IEEE TIA 2025.
%
% Time-domain figures (6, 11, 12, 14):
%   Full 9-state nonlinear GFL ODE integrated with ode15s (run_gfl_full_nonlinear).
%   No linearisation, no Simulink, no paper Table V used for waveforms.
%
% Eigenvalue figures (8, 10):
%   pchip interpolation of first-principles table over gSCR sweep; 39-node
%   network eigenvalues computed from calibrated Z matrix.
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

% ── First-principles eigenvalue table (replaces hardcoded Table V) ────────
% Sweeps gSCR = 1.5 … 5.0 using the 9-state GFL ODE with GFF filter.
% This is the gSCR -> (sigma, omega) calibration curve derived independently
% from the control model, without using paper Table V as input.
fprintf('\nComputing first-principles eigenvalue table (gSCR 1.5 to 5.0)...\n');
fprintf('(~15 fsolve+Jacobian evaluations; takes ~1 min)\n');
fp_gscr_sweep = [2.4, 2.55, linspace(2.65, 3.0, 5), linspace(3.2, 5.0, 5)];
fp_raw  = compute_first_principles_eigenvalues(fp_gscr_sweep);
eig_table.gSCR  = fp_raw.gSCR(:)';
eig_table.sigma = fp_raw.sigma(:)';
eig_table.omega = fp_raw.omega(:)';
CgSCR = interp1(eig_table.sigma, eig_table.gSCR, 0, 'pchip');   % gSCR where sigma=0
fprintf('First-principles critical gSCR (sigma=0): %.4f  (paper reports ~2.66)\n', CgSCR);

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
%  FIG. 8  Loci of weakest eigenvalues (single-converter, N=1)
%
%  Matches paper Fig. 8:
%    CBR locus: s=+1  (injecting), Z11=1/SCR, SCR: 4.0→1.5  (rightward)
%    SE  locus: s=-1  (absorbing), Z11=1/|SCR|, |SCR|: 1.5→4.0 (rightward)
%  Both plotted as discrete markers: SE=blue stars, CBR=red triangles.
% ═══════════════════════════════════════════════════════════════════════
fprintf('\n── Fig. 8: Eigenvalue loci (single-converter N=1) ──\n');

% -- CBR locus: s=+1, SCR sweeps 4.0 → 1.5 ---------------------------
scr_cbr = linspace(4.0, 1.5, 20);
sigma_cbr = nan(size(scr_cbr));
omega_cbr = nan(size(scr_cbr));
fprintf('  CBR locus (%d pts)...', numel(scr_cbr));
for i8 = 1:numel(scr_cbr)
    [ok8, sigma_cbr(i8), omega_cbr(i8)] = ...
        compute_single_conv_eig(1/scr_cbr(i8), 1.0, ctrl);
    if ~ok8, sigma_cbr(i8) = NaN; omega_cbr(i8) = NaN; end
end
fprintf(' done.\n');
valid_cbr = ~isnan(sigma_cbr);

% -- SE locus: s=-1, |SCR| sweeps 1.5 → 4.0 --------------------------
scr_se = linspace(1.5, 4.0, 20);
sigma_se = nan(size(scr_se));
omega_se = nan(size(scr_se));
fprintf('  SE  locus (%d pts)...', numel(scr_se));
for i8 = 1:numel(scr_se)
    [ok8, sigma_se(i8), omega_se(i8)] = ...
        compute_single_conv_eig(1/scr_se(i8), -1.0, ctrl);
    if ~ok8, sigma_se(i8) = NaN; omega_se(i8) = NaN; end
end
fprintf(' done.\n');
valid_se = ~isnan(sigma_se);

fig8 = figure('Color','w','Name','Fig 8','Position',[60 60 620 480],'Visible','off');
hold on;

% SE: blue filled stars (pentagram)
scatter(sigma_se(valid_se),  omega_se(valid_se),  70, 'p', ...
    'MarkerFaceColor',[0.15 0.45 0.80], 'MarkerEdgeColor',[0.10 0.30 0.65], ...
    'LineWidth',0.8, 'DisplayName','When absorbing active power from the grid');
scatter(sigma_se(valid_se), -omega_se(valid_se),  70, 'p', ...
    'MarkerFaceColor',[0.15 0.45 0.80], 'MarkerEdgeColor',[0.10 0.30 0.65], ...
    'LineWidth',0.8, 'HandleVisibility','off');

% CBR: red open triangles
scatter(sigma_cbr(valid_cbr),  omega_cbr(valid_cbr), 70, '^', ...
    'MarkerEdgeColor',[0.80 0.10 0.10], 'MarkerFaceColor','none', ...
    'LineWidth',1.0, 'DisplayName','When injecting active power into the grid');
scatter(sigma_cbr(valid_cbr), -omega_cbr(valid_cbr), 70, '^', ...
    'MarkerEdgeColor',[0.80 0.10 0.10], 'MarkerFaceColor','none', ...
    'LineWidth',1.0, 'HandleVisibility','off');

xline(0, '--k', 'LineWidth', 0.8);
ax8 = gca;
grid off;
xlabel('Real Part');
ylabel('Imaginary Part');

% Dashed boxes around each cluster
se_x = sigma_se(valid_se);  se_w = omega_se(valid_se);
if ~isempty(se_x)
    x1 = min(se_x)-2; x2 = max(se_x)+2;
    y1 = -(max(se_w)+15); y2 = max(se_w)+15;
    rectangle('Position',[x1,y1,x2-x1,y2-y1], ...
              'EdgeColor',[0.3 0.3 0.3],'LineStyle','--','LineWidth',0.9);
    % Direction arrow: left→right (SCR: -1.5→-4.0)
    annotation('textarrow',[0.18 0.30],[0.50 0.50], ...
        'String','','HeadWidth',6,'HeadLength',6,'LineWidth',0.9,'Units','normalized');
    text(ax8, mean([x1 x2]), 0, 'SCR: -1.5\rightarrow-4.0', ...
         'FontSize',8,'HorizontalAlignment','center');
end
cbr_x = sigma_cbr(valid_cbr);  cbr_w = omega_cbr(valid_cbr);
if ~isempty(cbr_x)
    x1 = min(cbr_x)-2; x2 = max(cbr_x)+2;
    y1 = -(max(cbr_w)+15); y2 = max(cbr_w)+15;
    rectangle('Position',[x1,y1,x2-x1,y2-y1], ...
              'EdgeColor',[0.3 0.3 0.3],'LineStyle','--','LineWidth',0.9);
    annotation('textarrow',[0.55 0.88],[0.50 0.50], ...
        'String','','HeadWidth',6,'HeadLength',6,'LineWidth',0.9,'Units','normalized');
    text(ax8, mean([x1 x2]), 0, 'SCR: 4.0\rightarrow1.5', ...
         'FontSize',8,'HorizontalAlignment','center');
end

legend('Location','north','Box','off','FontSize',8);
set(ax8,'Color','w','Box','on','XColor','k','YColor','k');
exportgraphics(fig8, fullfile(resultsDir,'fig8_eigenvalue_loci.png'), ...
    'Resolution',180,'BackgroundColor','white');
close(fig8);
fprintf('  → fig8_eigenvalue_loci.png\n');

% ═══════════════════════════════════════════════════════════════════════
%  FIG. 10  IEEE 39-node: eigenvalue loci as m varies
%
%  Computes the actual N=9 gfl13 Jacobian (117×117) at each m value.
%  Panel (a): all 117 eigenvalues at m=1.015 (gSCR≈2.66)
%  Panel (b): PLL-band eigenvalues [70,130 rad/s] as discrete markers
% ═══════════════════════════════════════════════════════════════════════
fprintf('\n── Fig. 10: 39-node eigenvalues ──\n');

s_f10     = ones(9,1);
s_f10(1:3) = 1.0 - SE39_cap;          % nodes 30-32: s_net=0.25
W_raw_f10  = diag(s_f10) * Z_9x9_raw;
k_f10      = (1/2.659) / max_positive_eig(W_raw_f10);
nSt_f10    = 13*9;

% Discrete m values: 10 points from well-damped to unstable
m_vals_f10 = [0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.015, 1.03, 1.05];
nM_f10     = numel(m_vals_f10);

sigma_pll_f10 = cell(nM_f10, 1);   % real parts of PLL-band eigs
omega_pll_f10 = cell(nM_f10, 1);   % imag parts of PLL-band eigs
ev_full_ref   = [];                 % all 117 eigs at m=1.015 for panel (a)

for ii = 1:nM_f10
    m_ii  = m_vals_f10(ii);
    Z_m   = (k_f10 * m_ii/1.015) * Z_9x9_raw;
    gscr_m = 1/max_positive_eig(diag(s_f10)*Z_m);
    fprintf('  m=%.3f (gSCR=%.3f) Jacobian...', m_ii, gscr_m);
    try
        x_ss = gfl13_find_ss(Z_m, s_f10, ctrl);
        f0   = gfl13_ode(x_ss, ctrl, Z_m, s_f10, 1.0);
        J_m  = zeros(nSt_f10);
        for jj = 1:nSt_f10
            xp = x_ss; xp(jj) = xp(jj) + 1e-7;
            J_m(:,jj) = (gfl13_ode(xp, ctrl, Z_m, s_f10, 1.0) - f0) / 1e-7;
        end
        ev_m = eig(J_m);
        if abs(m_ii - 1.015) < 0.001
            ev_full_ref = ev_m;          % save full spectrum for panel (a)
        end
        mask_pll = abs(imag(ev_m)) >= 70 & abs(imag(ev_m)) <= 130;
        sigma_pll_f10{ii} = real(ev_m(mask_pll));
        omega_pll_f10{ii} = imag(ev_m(mask_pll));
        fprintf(' %d PLL eigs\n', sum(mask_pll));
    catch e
        fprintf(' failed: %s\n', e.message);
    end
end

fig10 = figure('Color','w','Name','Fig 10','Position',[80 60 920 420],'Visible','off');

% ── Panel (a): all eigenvalues at m=1.015 ─────────────────────────────
subplot(1,2,1);
hold on;
if ~isempty(ev_full_ref)
    % All non-PLL eigenvalues in grey
    mask_bg = abs(imag(ev_full_ref)) < 70 | abs(imag(ev_full_ref)) > 130;
    scatter(real(ev_full_ref(mask_bg)),  imag(ev_full_ref(mask_bg)), ...
            14, [0.55 0.55 0.55], 'x', 'LineWidth', 0.8);
    % PLL eigenvalues in red
    mask_a = abs(imag(ev_full_ref)) >= 70 & abs(imag(ev_full_ref)) <= 130;
    ev_pll_a = ev_full_ref(mask_a);
    scatter(real(ev_pll_a), imag(ev_pll_a), 40, [0.8 0 0], 'x', 'LineWidth', 1.5);
end
xline(0,'--k','LineWidth',0.8);
xlabel('Real Part'); ylabel('Imaginary Part');
title('(a)');
set(gca,'Color','w','Box','on'); grid off;

% ── Panel (b): PLL eigenvalues, discrete scatter markers ──────────────
subplot(1,2,2);
hold on;
for ii = 1:nM_f10
    sig = sigma_pll_f10{ii};
    omg = omega_pll_f10{ii};
    if isempty(sig), continue; end
    % Blue-outlined red-filled squares (matching paper style)
    scatter(sig,  omg, 45, 's', 'MarkerEdgeColor',[0 0 0.75], ...
            'MarkerFaceColor',[0.85 0.08 0.08], 'LineWidth', 0.7);
    scatter(sig, -omg, 45, 's', 'MarkerEdgeColor',[0 0 0.75], ...
            'MarkerFaceColor',[0.85 0.08 0.08], 'LineWidth', 0.7);
end
xline(0,'--k','LineWidth',1.0);

% Direction arrow annotation (left to right = gSCR 3.86 → 2.57)
ax2 = gca;
xl  = xlim(ax2);
annotation('textarrow', [0.58 0.78], [0.62 0.62], 'String','', ...
    'HeadWidth',6,'HeadLength',6,'LineWidth',1.0, ...
    'Units','normalized');
text(ax2,-5.2, 98,'gSCR: 3.86\rightarrow2.57','FontSize',8,'Color','k');

% Highlight gSCR≈2.66 (m=1.015) cluster
sig_crit = sigma_pll_f10{abs(m_vals_f10 - 1.015) < 0.001};
if ~isempty(sig_crit)
    x1 = min(sig_crit)-0.25;  x2 = max(sig_crit)+0.25;
    rectangle('Position',[x1, 68, x2-x1, 130-68+4], ...
              'EdgeColor','r','LineWidth',1.2,'LineStyle','-');
    text(ax2, mean(sig_crit), -102, 'gSCR=2.66', ...
         'Color','r','FontSize',8,'HorizontalAlignment','center');
end

xlabel('Real Part'); ylabel('Imaginary Part');
title('(b)');
set(ax2,'Color','w','Box','on'); grid off;

exportgraphics(fig10, fullfile(resultsDir,'fig10_39node_eigenvalues.png'), ...
    'Resolution',180,'BackgroundColor','white');
close(fig10);
fprintf('  → fig10_39node_eigenvalues.png\n');

% ═══════════════════════════════════════════════════════════════════════
%  TWO-AREA FOUR-CBR  →  Fig. 6
%  Same 4-CBR system at three gSCR levels: 2.66 (unstable), 2.95, 3.57.
%  Power is normalised by s_ref so all converters start at 1.0 pu.
% ═══════════════════════════════════════════════════════════════════════
fprintf('\n── Two-area network (Fig. 6) ──\n');

W_paper = [0.223 0.069 0.015 0.015;
           0.139 0.238 0.030 0.030;
           0.045 0.045 0.249 0.130;
           0.015 0.015 0.043 0.246];
s_cbr_2a  = [0.5; 1.0; 1.5; 0.5];
Z_base_2a = diag(1./s_cbr_2a) * W_paper;
Z_base_2a = (Z_base_2a + Z_base_2a') / 2;   % symmetrise
W_base_2a = diag(s_cbr_2a) * Z_base_2a;
lam_base_2a  = max_positive_eig(W_base_2a);
gSCR_base_2a = 1 / lam_base_2a;
fprintf('  Base two-area gSCR = %.3f\n', gSCR_base_2a);

gSCR_f6     = [2.66, 2.95, 3.57];
panel_lbl   = {'(a)','(b)','(c)'};
results_f6  = cell(3,1);

for ii = 1:3
    k_ii   = gSCR_base_2a / gSCR_f6(ii);   % scale Z → gSCR = gSCR_f6(ii)
    Z_ii   = k_ii * Z_base_2a;
    gscr_v = 1 / max_positive_eig(diag(s_cbr_2a) * Z_ii);
    fprintf('  Simulating gSCR=%.2f (actual=%.3f)...\n', gSCR_f6(ii), gscr_v);
    [t_ii, P_ii, ~] = run_gfl_full_nonlinear(s_cbr_2a, Z_ii, ctrl, 0.10, 1.0);
    P_norm = P_ii ./ s_cbr_2a(:)';          % normalise: each column / s_ref_k
    results_f6{ii} = struct('t', t_ii, 'P', P_norm, 'gscr', gscr_v);
end

fig6_path = fullfile(resultsDir, 'fig6_two_area_time_domain.png');
f6_cols   = lines(4);
cbr_names = {'CBR1','CBR2','CBR3','CBR4'};

fig6 = figure('Color','w','Name','Fig 6','Position',[80 60 480 580],'Visible','off');
tiledlayout(fig6, 3, 1, 'Padding','compact','TileSpacing','compact');
for ii = 1:3
    r  = results_f6{ii};
    ax = nexttile;
    hold(ax,'on');
    for j = 1:4
        plot(ax, r.t, r.P(:,j), 'LineWidth', 0.9, 'Color', f6_cols(j,:));
    end
    set(ax,'Color','w','Box','on','XColor','k','YColor','k');
    grid(ax,'off');
    xlim(ax,[0 1.0]);
    ylabel(ax,'P/(p.u.)');
    text(ax,0.97,0.08,sprintf('gSCR=%.2f',r.gscr),'Units','normalized',...
         'HorizontalAlignment','right','FontSize',9);
    text(ax,0.03,0.08,panel_lbl{ii},'Units','normalized',...
         'HorizontalAlignment','left','FontSize',9,'FontWeight','bold');
    if ii == 1
        legend(ax, cbr_names{:}, 'Location','north','Orientation','horizontal',...
               'FontSize',8,'Box','off');
    end
    if ii == 3, xlabel(ax,'time/s'); end
end
exportgraphics(fig6, fig6_path, 'Resolution',180,'BackgroundColor','white');
close(fig6);
fprintf('  → fig6_two_area_time_domain.png\n');

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
end

fig11_path = fullfile(resultsDir, 'fig11_39node_gscr_time_domain.png');
simulate_and_plot(cases39_fig11, ctrl, 0.10, 1.0, ...
    'Fig. 11: 39-node responses at different gSCR levels', fig11_path);

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
    'Case 1: no SE',           Z_c1, s_c1, true(9,1)
    'Case 2: 3 SEs at node 23',Z_c2, s_c2, [true(9,1);false]
    'Case 3: SEs at 3,22,35',  Z_c3, s_c3, [true(9,1);false;false]
    'Case 4: SEs at 37,38,39', Z_c4, s_c4, true(9,1)
};

fig12_path = fullfile(resultsDir, 'fig12_39node_placement_time_domain.png');
simulate_and_plot(cases39_fig12, ctrl, 0.10, 1.0, ...
    'Fig. 12: 39-node SE placement comparison', fig12_path);

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
        'Case 1: no SE',                     Z33_c1, s33_c1, true(33,1)
        'Case 2: non-optimal passive SEs',   Z33_c2, s33_c2, mask33_c2
        'Case 3: greedy-optimal passive SEs',Z33_c3, s33_c3, mask33_c3
    };

    fig14_path = fullfile(resultsDir, 'fig14_33converter_time_domain.png');
    simulate_and_plot(cases33, ctrl, 0.25, 1.0, ...
        'Fig. 14: 33-converter placement comparison', fig14_path);
else
    fprintf('  33-converter matrices not found – run import_reproduction_data first.\n');
    fprintf('  Expected: %s\n', mat33_path);
end

fprintf('\nDone.  Figures written to:\n  %s\n', resultsDir);
end


% ═══════════════════════════════════════════════════════════════════════
%  compute_single_conv_eig  –  dominant eigenvalue for N=1 GFL converter
%
%  Z11   – scalar network impedance seen by the single converter
%  S_ref – signed capacity: +1 for CBR, -se_cap for SE
%  ctrl  – struct from gfl_control_params()
%
%  Returns ok=true and (sigma, omega) if SS converges; NaN otherwise.
% ═══════════════════════════════════════════════════════════════════════
function [ok, sigma, omega] = compute_single_conv_eig(Z11, S_ref, ctrl)
ok = false; sigma = NaN; omega = NaN;
try
    x_ss = gfl13_find_ss(Z11, S_ref, ctrl);
    f0   = gfl13_ode(x_ss, ctrl, Z11, S_ref, 1.0);
    J    = zeros(13);
    eps_fd = 1e-7;
    for j = 1:13
        xp = x_ss; xp(j) = xp(j) + eps_fd;
        J(:,j) = (gfl13_ode(xp, ctrl, Z11, S_ref, 1.0) - f0) / eps_fd;
    end
    ev = eig(J);
    mask = abs(imag(ev)) >= 70 & abs(imag(ev)) <= 130;
    ev_pll = ev(mask);
    if isempty(ev_pll)
        [~, bi] = min(abs(abs(imag(ev)) - 88));
        lam = ev(bi);
    else
        [~, bi] = max(real(ev_pll));
        lam = ev_pll(bi);
    end
    sigma = real(lam);
    omega = abs(imag(lam));
    ok = true;
catch
end
end


% ═══════════════════════════════════════════════════════════════════════
%  simulate_and_plot  –  9-state nonlinear ODE simulation, plot CBR powers
%
%  cases: {label, Z_active, s_active, cbr_mask}
%    cbr_mask – logical(N,1): true for CBR nodes to include in the plot
%
%  Two-phase design: run all ODE integrations first, then create and export
%  the figure.  Avoids batch-mode figure-handle invalidation during long
%  computations (R2025a issue with exportgraphics on long-lived figure handles).
% ═══════════════════════════════════════════════════════════════════════
function simulate_and_plot(cases, ctrl, u_mag, t_stop, figTitle, outPath)

nCases = size(cases, 1);

% ── Phase 1: run all simulations (no figure open) ────────────────────
results = cell(nCases, 1);
for k = 1:nCases
    label    = cases{k,1};
    Z_active = cases{k,2};
    s_active = cases{k,3}(:);
    cbr_mask = cases{k,4}(:);
    [t, P, info] = run_gfl_full_nonlinear(s_active, Z_active, ctrl, u_mag, t_stop);
    % Normalise by s_ref so each converter starts at 1.0 p.u.
    s_cbr    = abs(s_active(cbr_mask));   % rated capacity of plotted CBRs
    P_norm   = P(:, cbr_mask) ./ s_cbr(:)';
    results{k} = struct('label', label, 't', t, ...
                        'P_cbr', P_norm, 'gscr', info.gscr);
end

% ── Phase 2: create a fresh figure and plot ───────────────────────────
fig = figure('Color','w','Name',figTitle, ...
             'Position',[100 100 780 200*nCases], 'Visible','off');
tiledlayout(fig, nCases, 1, 'Padding','compact','TileSpacing','compact');

for k = 1:nCases
    r    = results{k};
    nCBR = size(r.P_cbr, 2);
    cols = lines(nCBR);

    ax = nexttile;
    hold(ax, 'on');
    for j = 1:nCBR
        plot(ax, r.t, r.P_cbr(:,j), 'LineWidth', 0.9, 'Color', cols(j,:));
    end
    ylabel(ax, 'P/(p.u.)');
    title(ax, sprintf('%s  [gSCR=%.3f]', r.label, r.gscr), 'Interpreter','none');
    set(ax, 'Color','w','XColor','k','YColor','k','Box','on');
    grid(ax, 'off');
    xlim(ax, [0 t_stop]);
    if k == nCases, xlabel(ax, 'time/s'); end
end
sgtitle(fig, figTitle, 'Color','k');
exportgraphics(fig, outPath, 'Resolution',180,'BackgroundColor','white');
close(fig);
fprintf('  → %s\n', outPath);
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


% gfl_control_params is defined in lib/gfl_control_params.m (added to path above).
