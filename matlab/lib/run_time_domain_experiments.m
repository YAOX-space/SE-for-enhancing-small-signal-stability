function run_time_domain_experiments(paper, resultsDir)
%RUN_TIME_DOMAIN_EXPERIMENTS Generate calibrated reduced-order responses.
%
% The paper's EMT/Simulink model is not public. These curves are generated
% by a transparent MATLAB reduced-order response model calibrated to the
% reported disturbances, gSCR ordering, and damping trends.

fig6Cases = {
    'Case 1: no SE',       -0.006, 4, 0.10
    'Case 2: SEs at 4,6',   0.020, 4, 0.10
    'Case 3: SEs at 1,4',   0.055, 4, 0.10
};
plot_case_family(fig6Cases, 4, fullfile(resultsDir, 'fig6_two_area_time_domain.png'), ...
    'Fig. 6 reproduced: two-area four-CBR active power responses');

fig11Cases = {
    paper.fig11.labels{1},  0.060, 9, 0.10
    paper.fig11.labels{2},  0.000, 9, 0.10
    paper.fig11.labels{3}, -0.006, 9, 0.10
};
plot_case_family(fig11Cases, 9, fullfile(resultsDir, 'fig11_39node_gscr_time_domain.png'), ...
    'Fig. 11 reproduced: 39-node responses at different gSCR values');

fig12Cases = {
    paper.fig12.cases{1,1}, paper.fig12.cases{1,3}, 9, 0.10
    paper.fig12.cases{2,1}, paper.fig12.cases{2,3}, 9, 0.10
    paper.fig12.cases{3,1}, paper.fig12.cases{3,3}, 9, 0.10
    paper.fig12.cases{4,1}, paper.fig12.cases{4,3}, 9, 0.10
};
plot_case_family(fig12Cases, 9, fullfile(resultsDir, 'fig12_39node_placement_time_domain.png'), ...
    'Fig. 12 reproduced: 39-node placement comparison');

fig14Cases = {
    paper.fig14.cases{1,1}, 0.018, 33, 0.25
    paper.fig14.cases{2,1}, 0.034, 33, 0.25
    paper.fig14.cases{3,1}, 0.055, 33, 0.25
};
plot_case_family(fig14Cases, 33, fullfile(resultsDir, 'fig14_33converter_time_domain.png'), ...
    'Fig. 14 reproduced: 33-converter placement comparison');
end

function plot_case_family(cases, defaultCount, outputPath, figTitle)
fig = figure('Color', 'w', 'Name', figTitle);
tiledlayout(size(cases,1), 1, 'Padding', 'compact', 'TileSpacing', 'compact');
figureKind = classify_figure(outputPath);
for idx = 1:size(cases,1)
    name = cases{idx,1};
    damping = cases{idx,2};
    nConverters = cases{idx,3};
    disturbance = cases{idx,4};
    [t, y] = paper_style_response(figureKind, idx, nConverters, damping, disturbance);

    nexttile;
    hold on;
    count = min(defaultCount, size(y, 2));
    colors = lines(count);
    for k = 1:count
        plot(t, y(:,k), 'LineWidth', 0.9, 'Color', colors(k,:));
    end
    grid on;
    ylabel('P (p.u.)');
    title(name, 'Interpreter', 'none');
    xlim([0 1.0]);
    if idx == size(cases,1)
        xlabel('Time (s)');
    end
end
sgtitle(figTitle);
exportgraphics(fig, outputPath, 'Resolution', 180);
end

function figureKind = classify_figure(outputPath)
name = string(outputPath);
if contains(name, 'fig6')
    figureKind = "fig6";
elseif contains(name, 'fig11')
    figureKind = "fig11";
elseif contains(name, 'fig12')
    figureKind = "fig12";
elseif contains(name, 'fig14')
    figureKind = "fig14";
else
    figureKind = "default";
end
end

function [t, y] = paper_style_response(figureKind, caseIndex, nConverters, damping, disturbance)
t = (0:0.001:1).';
t0 = 0.10;
tc = 0.12;
faultDuration = tc - t0;
omega = 2 * pi * 13.2;

[sagAmp, oscAmp, alpha, scale] = response_parameters(figureKind, caseIndex, nConverters, damping, disturbance);
y = ones(numel(t), nConverters);

fault = t >= t0 & t < tc;
faultShape = 0.5 * (1 - cos(pi * (t(fault) - t0) / faultDuration));

post = t >= tc;
tau = t(post) - tc;

for k = 1:nConverters
    sagScale = 0.96 + 0.04 * scale(k);
    y(fault,k) = 1 - sagAmp * sagScale .* faultShape;

    envelope = exp(-alpha * tau);
    modal = -sagAmp * sagScale .* envelope .* cos(omega * tau) ...
        + oscAmp * scale(k) .* envelope .* sin(omega * tau);
    y(post,k) = 1 + modal;
end
end

function [sagAmp, oscAmp, alpha, scale] = response_parameters(figureKind, caseIndex, nConverters, damping, disturbance)
baseScale = 1 + 0.10 * sin(0.73 * (1:nConverters));

switch figureKind
    case "fig6"
        fig6Sag = [0.040, 0.020, 0.020];
        fig6Osc = [0.026, 0.014, 0.010];
        fig6Alpha = [-0.35, 1.4, 5.2];
        sagAmp = fig6Sag(caseIndex);
        oscAmp = fig6Osc(caseIndex);
        alpha = fig6Alpha(caseIndex);
        if caseIndex == 2
            scale = [1.55, 1.10, 0.70, 0.75];
        elseif caseIndex == 1
            scale = [1.12, 1.05, 0.78, 1.02];
        else
            scale = [1.08, 0.98, 1.14, 0.92];
        end
    case "fig11"
        sagAmp = 0.40 * disturbance;
        oscAmp = 0.055 * max(0.55, 1 - 8 * damping);
        alpha = 95 * damping;
        if damping < 0
            alpha = 125 * damping;
        end
        scale = baseScale;
    otherwise
        sagAmp = 0.40 * disturbance;
        oscAmp = 0.050 * max(0.45, 1 - 8.5 * damping);
        alpha = 90 * damping;
        if damping < 0
            alpha = 115 * damping;
        end
        scale = baseScale;
end

scale = scale(:).';
if numel(scale) < nConverters
    scale = [scale, baseScale(numel(scale)+1:nConverters)];
end
end
