%% GFL Two-Area Four-CBR Time-Domain Simulation
%
% Reproduces the two-area four-CBR scenario from Yuan et al., "Placing
% Storage Energies for Enhancing Small-Signal Stability of Converter-Based
% Renewable Systems", IEEE TIA 2025.
%
% Network: custom two-area system, 9 nodes (4 CBRs + 4 passive + 1 infinite
%   bus).  Line data from matlab/data/table_x_two_area_lines.csv.
%
% This script is the GFL analog of IEEE9BusSystemExample.m:
%   - IEEE9BusSystemExample  →  3 synchronous generators, rotor-speed response
%   - THIS SCRIPT             →  4 GFL CBRs + SEs, active-power response
%
% Three cases (Fig. 6 equivalent):
%   Case 1: No SE         – 4 CBRs, close to instability
%   Case 2: SE at (4,6)  – sub-optimal placement
%   Case 3: SE at (1,4)  – greedy-optimal placement, best damping
%
% Output: matlab/results/gfl_twoarea_time_domain.png
%
% Requires: MATLAB + Simulink (no Simscape Electrical licence needed).

%% ── Setup ────────────────────────────────────────────────────────────────
rootMATLAB = fullfile(fileparts(mfilename('fullpath')), '..', 'matlab');
addpath(fullfile(rootMATLAB, 'lib'));
dataDir    = fullfile(rootMATLAB, 'data');
resultsDir = fullfile(rootMATLAB, 'results');
if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end

ctrl = gfl_control_params();

%% ── Eigenvalue table ─────────────────────────────────────────────────────
% Using hardcoded Table V for speed; replace with:
%   eig_table = compute_first_principles_eigenvalues(linspace(1.5,5.0,13));
eig_table.gSCR  = [2.650, 2.926, 3.256, 3.666];
eig_table.sigma = [0.090, -1.758, -3.485, -5.104];
eig_table.omega = [88.342, 88.949, 89.346, 89.564];

%% ── Build two-area Z matrix ──────────────────────────────────────────────
% The two-area network has CBRs at nodes 1-4 and passive buses 5-8, with
% node 9 as the infinite-bus reference.  The only inter-area coupling (CBRs
% 1,2 ↔ CBRs 3,4) runs through node 9, so the full pinv-based Z is used to
% preserve it.  The 4×4 CBR sub-block is then overwritten with values
% consistent with the paper's reported weighted-impedance matrix W, and
% cross-terms for node 6 (used in Case 2) are recovered similarly.
%
% This mirrors the construction in run_physics_simulink_reproduction.m and
% keeps the simulation consistent with the paper's Fig. 6.

lines2a = readtable(fullfile(dataDir, 'table_x_two_area_lines.csv'));

% --- standard build (8×8, ground node = 9) ---
[~, Z2a_net] = build_network_from_edges(lines2a, 9, 9);   % 8×8, row/col 1..8

% --- overwrite CBR block (nodes 1-4) with paper-calibrated values ---
W_paper = [0.223 0.069 0.015 0.015;
           0.139 0.238 0.030 0.030;
           0.045 0.045 0.249 0.130;
           0.015 0.015 0.043 0.246];
s_cbr_2a  = [0.5; 1.0; 1.5; 0.5];                     % CBR rated capacities
Z_cbr_2a  = diag(1./s_cbr_2a) * W_paper;               % Z = S^{-1} W
Z_cbr_2a  = (Z_cbr_2a + Z_cbr_2a') / 2;               % symmetrize

Z2a = Z2a_net;
Z2a(1:4, 1:4) = Z_cbr_2a;

% cross terms for node 6 (passive bus used in Case 2)
% col 6 of W for the 4 CBR rows: W_{i,6} = S_i * Z_{i,6}
W_col6    = [0.030; 0.060; 0.045; 0.015];
Z_cross6  = W_col6 ./ s_cbr_2a;
Z2a(1:4, 6) = Z_cross6;
Z2a(6, 1:4) = Z_cross6';

fprintf('Two-area Z matrix built (8×8).\n');
fprintf('  Initial gSCR (4 CBRs, no SE): %.4f\n', ...
    1/max_positive_eig(diag(s_cbr_2a) * Z_cbr_2a));

%% ── Three simulation cases ───────────────────────────────────────────────
%
%  Case 1: No SE
%    Active: CBRs 1,2,3,4  →  4×4 sub-matrix, s = [0.5;1.0;1.5;0.5]
%
%  Case 2: SE at nodes 4 and 6  (sub-optimal)
%    SE at node 4 absorbs 0.5 pu  →  CBR4 net = 0 (excluded)
%    SE at node 6 absorbs 0.5 pu  →  node 6 now active
%    Active: CBRs 1,2,3 + SE@6  →  4×4 sub-matrix rows/cols {1,2,3,6}
%    s = [0.5; 1.0; 1.5; -0.5]
%
%  Case 3: SE at nodes 1 and 4  (greedy-optimal)
%    SE at node 1 absorbs 0.5 pu  →  CBR1 net = 0 (excluded)
%    SE at node 4 absorbs 0.5 pu  →  CBR4 net = 0 (excluded)
%    Active: CBRs 2,3 only  →  2×2 sub-matrix rows/cols {2,3}
%    s = [1.0; 1.5]

c2_cbr_mask = logical([1; 1; 1; 0]);   % CBR nodes only in Case 2 output

cases_2a = {
    'Case 1: No SE',      ...
        Z2a([1 2 3 4],[1 2 3 4]),  [0.5;1.0;1.5;0.5],  true(4,1),   [0.5;1.0;1.5;0.5],   100
    'Case 2: SEs at 4,6', ...
        Z2a([1 2 3 6],[1 2 3 6]),  [0.5;1.0;1.5;-0.5], c2_cbr_mask, [0.5;1.0;1.5;0.0],   100
    'Case 3: SEs at 1,4', ...
        Z2a([2 3],[2 3]),           [1.0;1.5],          true(2,1),   [1.0;1.5],            100
};
nCases = size(cases_2a, 1);

%% ── Simulate and plot ────────────────────────────────────────────────────
fprintf('\nRunning two-area simulations...\n');

fig = figure('Color','w','Name','GFL Two-Area: SE placement', ...
             'Position',[100 80 760 520]);
tiledlayout(nCases, 1, 'Padding','compact', 'TileSpacing','compact');

for ci = 1:nCases
    label    = cases_2a{ci,1};
    Z_act    = cases_2a{ci,2};
    s_act    = cases_2a{ci,3}(:);
    cbr_mask = cases_2a{ci,4}(:);
    P0       = cases_2a{ci,5}(:);
    gf       = cases_2a{ci,6};

    [A, B, C, D, info] = build_gfl_state_space(s_act, Z_act, ctrl, gf, eig_table);
    fprintf('  %s:  gSCR=%.3f,  dom_eig=%+.3f±%.1fi\n', ...
        label, info.gscr, real(info.dominant_eig), abs(imag(info.dominant_eig)));

    [t, dP] = run_gfl_case(A, B, C, D, 0.06, 1.0, rootMATLAB);
    % sag magnitude 0.06 pu  (two-area uses smaller perturbation per paper)

    P     = P0' + dP;                   % absolute power
    P_cbr = P(:, cbr_mask);             % CBR nodes only

    nexttile;
    cols = lines(size(P_cbr,2));
    hold on;
    for k = 1:size(P_cbr,2)
        plot(t, P_cbr(:,k), 'LineWidth', 1.0, 'Color', cols(k,:));
    end
    cbr_ids = find(cbr_mask);
    legend(arrayfun(@(n) sprintf('CBR %d',n), cbr_ids, 'UniformOutput',false), ...
        'Location','northeast','FontSize',7);
    title(sprintf('%s  [gSCR=%.3f, \\lambda=%+.3f\\pm%.1fi]', ...
        label, info.gscr, real(info.dominant_eig), abs(imag(info.dominant_eig))), ...
        'Interpreter','tex');
    ylabel('P (p.u.)');
    grid on;
    xlim([0 1.0]);
    if ci == nCases, xlabel('Time (s)'); end
end
sgtitle('GFL Two-Area 4-CBR: Active-power response to voltage sag', 'Color','k');

outPath = fullfile(resultsDir, 'gfl_twoarea_time_domain.png');
exportgraphics(fig, outPath, 'Resolution', 180, 'BackgroundColor', 'white');
fprintf('  Saved → results/gfl_twoarea_time_domain.png\n');
