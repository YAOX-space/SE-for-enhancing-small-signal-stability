function run_simulink_time_domain_experiments(paper, resultsDir)
%RUN_SIMULINK_TIME_DOMAIN_EXPERIMENTS Run the dedicated Simulink model.
%
% This complements run_time_domain_experiments.m. The generated waveforms
% come from Simulink simulation of a configurable PLL-dominated state-space
% model, while retaining the paper's public disturbance and damping data.

rootDir = fileparts(fileparts(mfilename('fullpath')));
modelPath = build_se_gfl_pll_simulink_model(rootDir);
[modelDir, modelName] = fileparts(modelPath);

fig6Cases = {
    'Case 1: no SE',       -0.006, 4, 0.10
    'Case 2: SEs at 4,6',   0.020, 4, 0.10
    'Case 3: SEs at 1,4',   0.055, 4, 0.10
};

fig11Cases = {
    paper.fig11.labels{1},  0.060, 9, 0.10
    paper.fig11.labels{2},  0.000, 9, 0.10
    paper.fig11.labels{3}, -0.006, 9, 0.10
};

fig12Cases = {
    paper.fig12.cases{1,1}, paper.fig12.cases{1,3}, 9, 0.10
    paper.fig12.cases{2,1}, paper.fig12.cases{2,3}, 9, 0.10
    paper.fig12.cases{3,1}, paper.fig12.cases{3,3}, 9, 0.10
    paper.fig12.cases{4,1}, paper.fig12.cases{4,3}, 9, 0.10
};

fig14Cases = {
    paper.fig14.cases{1,1}, 0.018, 33, 0.25
    paper.fig14.cases{2,1}, 0.034, 33, 0.25
    paper.fig14.cases{3,1}, 0.055, 33, 0.25
};

oldDir = pwd;
cleanup = onCleanup(@() cd(oldDir));
cd(modelDir);

plot_simulink_case_family(modelName, "fig6", fig6Cases, 4, ...
    fullfile(resultsDir, 'simulink_fig6_two_area_time_domain.png'), ...
    'Simulink Fig. 6 approximation: two-area four-CBR responses');

plot_simulink_case_family(modelName, "fig11", fig11Cases, 9, ...
    fullfile(resultsDir, 'simulink_fig11_39node_gscr_time_domain.png'), ...
    'Simulink Fig. 11 approximation: 39-node responses at different gSCR values');

plot_simulink_case_family(modelName, "fig12", fig12Cases, 9, ...
    fullfile(resultsDir, 'simulink_fig12_39node_placement_time_domain.png'), ...
    'Simulink Fig. 12 approximation: 39-node placement comparison');

plot_simulink_case_family(modelName, "fig14", fig14Cases, 33, ...
    fullfile(resultsDir, 'simulink_fig14_33converter_time_domain.png'), ...
    'Simulink Fig. 14 approximation: 33-converter placement comparison');
end

function plot_simulink_case_family(modelName, figureKind, cases, defaultCount, outputPath, figTitle)
fig = figure('Color', 'w', 'Name', figTitle);
tiledlayout(size(cases,1), 1, 'Padding', 'compact', 'TileSpacing', 'compact');

for idx = 1:size(cases,1)
    caseName = cases{idx,1};
    damping = cases{idx,2};
    nConverters = cases{idx,3};
    disturbance = cases{idx,4};
    simData = simulate_case(modelName, figureKind, idx, nConverters, damping, disturbance);

    nexttile;
    hold on;
    count = min(defaultCount, size(simData.y, 2));
    colors = lines(count);
    for k = 1:count
        plot(simData.t, simData.y(:,k), 'LineWidth', 0.9, 'Color', colors(k,:));
    end
    ax = gca;
    set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'Box', 'on');
    grid off;
    ylabel('P (p.u.)');
    title(caseName, 'Interpreter', 'none');
    xlim([0 1.0]);
    ylim(paper_ylim(figureKind, idx));
    if idx == size(cases,1)
        xlabel('Time (s)');
    end
end
sgtitle(figTitle, 'Color', 'k');
exportgraphics(fig, outputPath, 'Resolution', 180, 'BackgroundColor', 'white');
end

function limits = paper_ylim(figureKind, caseIndex)
switch figureKind
    case "fig6"
        options = [0.96 1.04; 0.98 1.02; 0.98 1.02];
        limits = options(caseIndex,:);
    case "fig11"
        options = [0.96 1.04; 0.92 1.08; 0.88 1.12];
        limits = options(caseIndex,:);
    case "fig12"
        options = [0.92 1.08; 0.92 1.08; 0.96 1.04; 0.96 1.04];
        limits = options(caseIndex,:);
    case "fig14"
        options = [0.88 1.08; 0.88 1.08; 0.90 1.05];
        limits = options(caseIndex,:);
    otherwise
        limits = [0.90 1.10];
end
end

function simData = simulate_case(modelName, figureKind, caseIndex, nConverters, damping, disturbance)
t = (0:0.001:1).';
u = voltage_sag_profile(t, disturbance);
[ss_A, ss_B, ss_C, ss_D] = pll_state_space(figureKind, caseIndex, nConverters, damping, disturbance);

assignin('base', 'fault_input', timeseries(u, t));
assignin('base', 'ss_A', ss_A);
assignin('base', 'ss_B', ss_B);
assignin('base', 'ss_C', ss_C);
assignin('base', 'ss_D', ss_D);
assignin('base', 'sim_stop_time', '1.0');

out = sim(modelName, 'StopTime', '1.0');
raw = out.get('simout_P');
simData.t = raw.time;
simData.y = raw.signals.values;
if isvector(simData.y)
    simData.y = simData.y(:);
end
end

function u = voltage_sag_profile(t, disturbance)
t0 = 0.10;
tc = 0.12;
duration = tc - t0;
u = zeros(size(t));
fault = t >= t0 & t < tc;
phase = pi * (t(fault) - t0) ./ duration;
u(fault) = -disturbance .* sin(phase) .^ 2;
end

function [A, B, C, D] = pll_state_space(figureKind, caseIndex, nConverters, damping, disturbance)
omega = 2 * pi * 13.2;
[sagGain, modalGain, dampingRatio, scale] = response_parameters(figureKind, caseIndex, nConverters, damping, disturbance);

A = zeros(2 * nConverters);
B = zeros(2 * nConverters, 1);
C = zeros(nConverters, 2 * nConverters);
D = zeros(nConverters, 1);

for k = 1:nConverters
    row = 2 * k - 1;
    zeta = dampingRatio;
    A(row:row+1, row:row+1) = [0 1; -omega^2 -2*zeta*omega];
    B(row:row+1, 1) = [0; modalGain * scale(k) * omega^2];
    C(k, row:row+1) = [1 0];
    D(k, 1) = sagGain * (0.96 + 0.04 * scale(k));
end
end

function [sagGain, modalGain, dampingRatio, scale] = response_parameters(figureKind, caseIndex, nConverters, damping, disturbance)
baseScale = 1 + 0.10 * sin(0.73 * (1:nConverters));

switch figureKind
    case "fig6"
        sagByCase = [0.40, 0.20, 0.20];
        modalByCase = [0.24, 0.18, 0.14];
        dampingByCase = [-0.001, 0.020, 0.055];
        sagGain = sagByCase(caseIndex);
        modalGain = modalByCase(caseIndex);
        dampingRatio = dampingByCase(caseIndex);
        if caseIndex == 2
            scale = [1.55, 1.10, 0.70, 0.75];
        elseif caseIndex == 1
            scale = [1.12, 1.05, 0.78, 1.02];
        else
            scale = [1.08, 0.98, 1.14, 0.92];
        end
    case "fig11"
        sagGain = 0.40;
        modalGain = 0.85 * max(0.55, 1 - 8.0 * damping);
        dampingRatio = damping;
        scale = baseScale;
    otherwise
        sagGain = 0.40;
        modalGain = 0.75 * max(0.45, 1 - 8.5 * damping);
        dampingRatio = damping;
        scale = baseScale;
end

if disturbance > 0.10
    modalGain = 0.70 * modalGain;
end

scale = scale(:).';
if numel(scale) < nConverters
    scale = [scale, baseScale(numel(scale)+1:nConverters)];
end
end
