%% GFL 9-CBR Time-Domain Simulation
%
% Reproduces the 9-converter scenarios from Yuan et al., "Placing Storage
% Energies for Enhancing Small-Signal Stability of Converter-Based Renewable
% Systems", IEEE TIA 2025.
%
% Network: IEEE 39-bus.  9 GFL CBRs are located at buses 30-35, 37-39.
%   The remaining buses are passive (loads / network nodes).
%
% This script is the GFL analog of IEEE39BusSystemExample.m:
%   - IEEE39BusSystemExample  →  synchronous generators, rotor-speed response
%   - THIS SCRIPT              →  GFL CBRs + SEs, active-power response
%
% Part A (Fig. 11 equivalent): fixed SE placement at nodes 30/31/32,
%   three gSCR levels achieved by scaling the Z matrix.
%
% Part B (Fig. 12 equivalent): fixed gSCR baseline, four SE placement
%   strategies compared.
%
% Outputs written to matlab/results/:
%   gfl_9cbr_vary_gscr.png
%   gfl_9cbr_vary_placement.png
%
% Requires: MATLAB + Simulink (no Simscape Electrical licence needed).

%% ── Setup ────────────────────────────────────────────────────────────────
rootMATLAB = fullfile(fileparts(mfilename('fullpath')), '..', 'matlab');
addpath(fullfile(rootMATLAB, 'lib'));
dataDir    = fullfile(rootMATLAB, 'data');
resultsDir = fullfile(rootMATLAB, 'results');
if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end

ctrl = gfl_control_params();

%% ── Load IEEE 39-bus network ─────────────────────────────────────────────
fprintf('Loading IEEE 39-bus network...\n');
lines39  = readtable(fullfile(dataDir, 'table_xi_ieee39_lines.csv'));
nbus39   = max([lines39.from; lines39.to]);
[~, Z39_full] = build_network_from_edges(lines39, nbus39, 36);   % 38×38

keep39 = setdiff(1:nbus39, 36);
b2r39  = zeros(1, nbus39);
b2r39(keep39) = 1:38;

SE39_cap    = 0.75;                              % SE rated capacity (pu)
cbr_nodes39 = [30 31 32 33 34 35 37 38 39];     % 9 CBR bus numbers
cbr_rows39  = b2r39(cbr_nodes39);               % row indices in Z39_full
Z_9x9_raw   = Z39_full(cbr_rows39, cbr_rows39); % 9×9 CBR sub-matrix

%% ── Eigenvalue table ─────────────────────────────────────────────────────
% Using hardcoded Table V (paper) for speed.
% For first-principles accuracy (no paper data), replace with:
%   eig_table = compute_first_principles_eigenvalues(linspace(1.5,5.0,13));
eig_table.gSCR  = [2.650, 2.926, 3.256, 3.666];
eig_table.sigma = [0.090, -1.758, -3.485, -5.104];
eig_table.omega = [88.342, 88.949, 89.346, 89.564];
fprintf('Eigenvalue table: Table V (4 points).\n');
fprintf('  → For independent check: eig_table = compute_first_principles_eigenvalues(...);\n\n');

%% ════════════════════════════════════════════════════════════════════════
%%  PART A  –  Varying gSCR (analog of Fig. 11)
%
%  SEs pre-installed at nodes 30, 31, 32 (co-located with CBRs).
%  Net signed capacity at those nodes = 1.0 - 0.75 = 0.25 pu.
%  gSCR is varied by scaling the Z matrix via multiplier m.
%    m = 0.700  →  gSCR ≈ 3.86  (well-damped)
%    m = 1.015  →  gSCR ≈ 2.66  (critical: sigma ≈ 0)
%    m = 1.050  →  gSCR ≈ 2.57  (unstable)
% ════════════════════════════════════════════════════════════════════════
fprintf('=== Part A: Varying gSCR ===\n');

s_A     = ones(9,1);
s_A(1:3) = 1.0 - SE39_cap;   % nodes 30,31,32: s_net = 0.25

% Calibrate so that m=1.015 gives gSCR = 2.659
k_A = (1/2.659) / max_positive_eig(diag(s_A) * Z_9x9_raw);

m_vals   = [0.700,              1.015,               1.050            ];
labels_A = {'gSCR \approx 3.86 (stable)', ...
            'gSCR \approx 2.66 (critical)', ...
            'gSCR \approx 2.57 (unstable)'};

fig_A = figure('Color','w','Name','GFL 9-CBR: Varying gSCR', ...
               'Position',[100 80 800 560]);
tiledlayout(3, 1, 'Padding','compact', 'TileSpacing','compact');

for ci = 1:3
    Zm = (k_A * m_vals(ci)/1.015) * Z_9x9_raw;
    [A, B, C, D, info] = build_gfl_state_space(s_A, Zm, ctrl, 200, eig_table);
    fprintf('  %s:  gSCR=%.3f,  dom_eig=%+.3f±%.1fi\n', ...
        labels_A{ci}, info.gscr, real(info.dominant_eig), abs(imag(info.dominant_eig)));

    [t, dP] = run_gfl_case(A, B, C, D, 0.10, 1.0, rootMATLAB);
    P = ones(1, 9) + dP;   % nominal P = 1.0 pu

    nexttile;
    cols = lines(9);
    hold on;
    for k = 1:9
        plot(t, P(:,k), 'LineWidth', 0.9, 'Color', cols(k,:));
    end
    title(sprintf('%s  [gSCR=%.3f, \\lambda=%+.3f\\pm%.1fi]', ...
        labels_A{ci}, info.gscr, real(info.dominant_eig), abs(imag(info.dominant_eig))), ...
        'Interpreter','tex');
    ylabel('P (p.u.)');
    grid on;
    xlim([0 1.0]);
    if ci == 3, xlabel('Time (s)'); end
end
sgtitle('GFL 9-CBR: Active-power response to voltage sag (varying gSCR)', 'Color','k');
exportgraphics(fig_A, fullfile(resultsDir, 'gfl_9cbr_vary_gscr.png'), ...
    'Resolution', 180, 'BackgroundColor', 'white');
fprintf('  Saved → results/gfl_9cbr_vary_gscr.png\n\n');

%% ════════════════════════════════════════════════════════════════════════
%%  PART B  –  Varying SE placement (analog of Fig. 12)
%
%  Z is calibrated so that gSCR = 2.650 with no SE (9 CBRs, s = ones(9,1)).
%  Four SE placement strategies are compared:
%    Case 1: No SE                        → gSCR ≈ 2.65 (unstable)
%    Case 2: 3 SEs clustered at node 23  → sub-optimal passive placement
%    Case 3: SEs spread at nodes 3,22,35 → mixed passive + collocated
%    Case 4: SEs collocated at 37,38,39  → optimal collocated placement
% ════════════════════════════════════════════════════════════════════════
fprintf('=== Part B: Varying SE placement ===\n');

k_B     = (1/2.650) / max_positive_eig(Z_9x9_raw);
Z39_cal = k_B * Z39_full;                             % full 38×38, calibrated

% Case 1: no SE  (9 CBRs, s = 1.0 each)
Z_B1 = Z39_cal(cbr_rows39, cbr_rows39);
s_B1 = ones(9,1);

% Case 2: 3 SEs at passive node 23 (total SE = 3×0.75 = 2.25 pu)
se2_row   = b2r39(23);
act2_rows = [cbr_rows39(:)', se2_row];
Z_B2      = Z39_cal(act2_rows, act2_rows);            % 10×10
s_B2      = [ones(9,1); -3*SE39_cap];                 % node 23: big SE

% Case 3: SEs at passive nodes 3,22 + collocated SE at CBR node 35
idx35    = find(cbr_nodes39 == 35);
s_B3_cbr = ones(9,1);  s_B3_cbr(idx35) = 1.0 - SE39_cap;   % net = 0.25
se3_rows = b2r39([3, 22]);
act3_rows = [cbr_rows39(:)', se3_rows(:)'];
Z_B3     = Z39_cal(act3_rows, act3_rows);             % 11×11
s_B3     = [s_B3_cbr; -SE39_cap; -SE39_cap];

% Case 4: SEs collocated at CBR nodes 37, 38, 39
idx379 = ismember(cbr_nodes39, [37 38 39]);
s_B4   = ones(9,1);  s_B4(idx379) = 1.0 - SE39_cap;  % net = 0.25
Z_B4   = Z39_cal(cbr_rows39, cbr_rows39);             % 9×9

% cbr_mask: which rows of dP correspond to CBR nodes (exclude SE-only nodes)
cases_B = {
    'Case 1: No SE',              Z_B1, s_B1, true(9,1),               ones(9,1)
    'Case 2: 3 SEs at node 23',   Z_B2, s_B2, [true(9,1); false],      [ones(9,1); 0]
    'Case 3: SEs at 3, 22, 35',   Z_B3, s_B3, [true(9,1);false;false], [ones(9,1);0;0]
    'Case 4: SEs at 37, 38, 39',  Z_B4, s_B4, true(9,1),               ones(9,1)
};
nB = size(cases_B, 1);

fig_B = figure('Color','w','Name','GFL 9-CBR: SE placement', ...
               'Position',[150 80 800 740]);
tiledlayout(nB, 1, 'Padding','compact', 'TileSpacing','compact');

for ci = 1:nB
    label    = cases_B{ci,1};
    Z_act    = cases_B{ci,2};
    s_act    = cases_B{ci,3}(:);
    cbr_mask = cases_B{ci,4}(:);
    P0       = cases_B{ci,5}(:);

    [A, B, C, D, info] = build_gfl_state_space(s_act, Z_act, ctrl, 200, eig_table);
    fprintf('  %s:  gSCR=%.3f,  dom_eig=%+.3f±%.1fi\n', ...
        label, info.gscr, real(info.dominant_eig), abs(imag(info.dominant_eig)));

    [t, dP] = run_gfl_case(A, B, C, D, 0.10, 1.0, rootMATLAB);
    P_cbr   = (P0' + dP);
    P_cbr   = P_cbr(:, cbr_mask);

    nexttile;
    cols = lines(size(P_cbr, 2));
    hold on;
    for k = 1:size(P_cbr,2)
        plot(t, P_cbr(:,k), 'LineWidth', 0.9, 'Color', cols(k,:));
    end
    title(sprintf('%s  [gSCR=%.3f, \\lambda=%+.3f\\pm%.1fi]', ...
        label, info.gscr, real(info.dominant_eig), abs(imag(info.dominant_eig))), ...
        'Interpreter','tex');
    ylabel('P (p.u.)');
    grid on;
    xlim([0 1.0]);
    if ci == nB, xlabel('Time (s)'); end
end
sgtitle('GFL 9-CBR: Active-power response to voltage sag (varying SE placement)', 'Color','k');
exportgraphics(fig_B, fullfile(resultsDir, 'gfl_9cbr_vary_placement.png'), ...
    'Resolution', 180, 'BackgroundColor', 'white');
fprintf('  Saved → results/gfl_9cbr_vary_placement.png\n');
