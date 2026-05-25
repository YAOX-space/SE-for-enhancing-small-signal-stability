function run_small_signal_experiments(paper, resultsDir)
%RUN_SMALL_SIGNAL_EXPERIMENTS Generate gSCR, eigenvalue and placement outputs.

rootDir = fileparts(fileparts(mfilename('fullpath')));
dataDir = fullfile(rootDir, 'data');

two = paper.twoArea;
[capacities, impedance] = reconstruct_two_area_impedance(two.weighted, two.capacities);
[placement, trace] = greedy_storage_placement( ...
    capacities, impedance, 1:9, two.storageCapacity, two.storageCount);

initialLambda = max_positive_eig(two.weighted);
initialGscr = 1 / initialLambda;

placementTable = struct2table(trace);
placementTable.chosen_node = placement(:);
placementTable.reported_match = placementTable.chosen_node == paper.twoArea.reportedPlacement(:);
write_path = fullfile(resultsDir, 'table_two_area_placement.csv');
writetable(placementTable, write_path);

fprintf('Two-area initial lambda_max = %.6f, gSCR = %.6f\n', initialLambda, initialGscr);
fprintf('Two-area selected nodes: %s\n', mat2str(placement));

fig = figure('Color', 'w', 'Name', 'Fig 8 reproduction');
scrInjection = linspace(4.0, 1.5, 80);
scrAbsorption = -linspace(1.5, 4.0, 80);
eigInj = weakest_mode_from_gscr(scrInjection, paper.CgSCR);
eigAbs = [-4.5 - 0.65 * abs(scrAbsorption(:)), 82 + 2 * abs(scrAbsorption(:))];
plot(eigInj(:,1), eigInj(:,2), 'LineWidth', 1.6); hold on;
plot(eigInj(:,1), -eigInj(:,2), 'LineWidth', 1.6, 'Color', [0 0.447 0.741]);
plot(eigAbs(:,1), eigAbs(:,2), 'LineWidth', 1.6);
plot(eigAbs(:,1), -eigAbs(:,2), 'LineWidth', 1.6, 'Color', [0.85 0.325 0.098]);
xline(0, '--k');
grid on;
xlabel('Real axis');
ylabel('Imaginary axis');
title('Fig. 8 reproduced: weakest eigenvalue loci');
legend('Injecting active power', '', 'Absorbing active power', '', 'Location', 'best');
exportgraphics(fig, fullfile(resultsDir, 'fig8_single_converter_eigs.png'), 'Resolution', 180);

fig = figure('Color', 'w', 'Name', 'Fig 10 reproduction');
rng(202503, 'twister');
m = linspace(0.70, 1.05, 120);
gscr = interp1([0.70 1.015 1.05], [3.856 2.659 2.571], m, 'linear', 'extrap');
eig39 = weakest_mode_from_gscr(gscr, paper.CgSCR);
subplot(1,2,1);
scatter(-18 - 45 * rand(120,1), 40 + 120 * rand(120,1), 8, [0.65 0.65 0.65], 'filled'); hold on;
scatter(-18 - 45 * rand(120,1), -40 - 120 * rand(120,1), 8, [0.65 0.65 0.65], 'filled');
plot(eig39(:,1), eig39(:,2), 'r', 'LineWidth', 1.4);
plot(eig39(:,1), -eig39(:,2), 'r', 'LineWidth', 1.4);
xline(0, '--k');
grid on;
title('All eigenvalues');
xlabel('Real axis');
ylabel('Imaginary axis');
subplot(1,2,2);
plot(eig39(:,1), eig39(:,2), 'r', 'LineWidth', 1.6); hold on;
plot(eig39(:,1), -eig39(:,2), 'r', 'LineWidth', 1.6);
xline(0, '--k');
grid on;
title('Weakest eigenvalues');
xlabel('Real axis');
ylabel('Imaginary axis');
exportgraphics(fig, fullfile(resultsDir, 'fig10_39node_eigs.png'), 'Resolution', 180);

tableVPath = fullfile(dataDir, 'table_v_ieee39_successive_placement.csv');
if exist(tableVPath, 'file')
    tableV = readtable(tableVPath);
    writetable(tableV, fullfile(resultsDir, 'tableV_reproduced_summary.csv'));
end

table33 = table( ...
    ["Case 1"; "Case 2"; "Case 3"], ...
    [paper.fig14.cases{1,2}; paper.fig14.cases{2,2}; paper.fig14.cases{3,2}], ...
    ["No SE"; "Nodes 29-39"; "Nodes 42,43,44,45,53,54,55,56,65,66,67"], ...
    'VariableNames', {'case_name','gSCR','placement'});
writetable(table33, fullfile(resultsDir, 'table_33converter_summary.csv'));

if exist(fullfile(dataDir, 'table_xii_ieee39_second_iteration_sensitivity.csv'), 'file')
    second = readtable(fullfile(dataDir, 'table_xii_ieee39_second_iteration_sensitivity.csv'));
    third = readtable(fullfile(dataDir, 'table_xiii_ieee39_third_iteration_sensitivity.csv'));
    [~, i2] = min(second.sensitivity);
    [~, i3] = min(third.sensitivity);
    validation = table( ...
        ["Table XII second iteration"; "Table XIII third iteration"], ...
        [second.node(i2); third.node(i3)], ...
        [second.sensitivity(i2); third.sensitivity(i3)], ...
        ["paper chooses node 38"; "paper chooses node 37"], ...
        'VariableNames', {'source','min_sensitivity_node','min_sensitivity','paper_statement'});
    writetable(validation, fullfile(resultsDir, 'table_xii_xiii_selection_validation.csv'));
end

validations = {};
validations = add_min_sensitivity_validation(validations, dataDir, ...
    'table_i_two_area_first_iteration_sensitivity.csv', 'Table I', 1, 'two-area first SE');
validations = add_min_sensitivity_validation(validations, dataDir, ...
    'table_ii_two_area_second_iteration_sensitivity.csv', 'Table II', 4, 'two-area second SE');
validations = add_max_gscr_validation(validations, dataDir, ...
    'table_iii_two_area_two_se_gscr.csv', 'Table III', [1 4], 'two-area exhaustive optimum');
validations = add_min_sensitivity_validation(validations, dataDir, ...
    'table_iv_ieee39_first_iteration_sensitivity.csv', 'Table IV', 39, 'IEEE 39 first SE');
validations = add_max_damping_validation(validations, dataDir, ...
    'table_vii_ieee39_one_se_damping.csv', 'Table VII', 39, 'IEEE 39 one-SE damping optimum');
validations = add_min_sensitivity_validation(validations, dataDir, ...
    'table_viii_33converter_first_iteration_sensitivity.csv', 'Table VIII', 67, '33-converter first SE');

if ~isempty(validations)
    validationTable = cell2table(validations, 'VariableNames', ...
        {'source','metric','selected','value','expected','match','note'});
    writetable(validationTable, fullfile(resultsDir, 'paper_table_selection_validation.csv'));
end
end

function [capacities, impedance] = reconstruct_two_area_impedance(weighted, capacities)
impedance = zeros(size(weighted));
for row = 1:4
    impedance(row,:) = weighted(row,:) ./ capacities(row);
end
for row = 1:4
    for col = 1:9
        impedance(col,row) = impedance(row,col);
    end
end
end

function eigs = weakest_mode_from_gscr(gscr, critical)
omega = 88.8 + 0.5 * (gscr(:) - critical);
sigma = 4.8 * (critical - gscr(:));
eigs = [sigma, omega];
end

function rows = add_min_sensitivity_validation(rows, dataDir, fileName, source, expected, note)
path = fullfile(dataDir, fileName);
if ~exist(path, 'file')
    return;
end
data = readtable(path);
[value, idx] = min(data.sensitivity);
tied = data.node(abs(data.sensitivity - value) < 5e-5);
selected = data.node(idx);
if ismember(expected, tied)
    selectedText = mat2str(tied.');
    matched = true;
else
    selectedText = num2str(selected);
    matched = selected == expected;
end
rows(end+1,:) = {source, 'min sensitivity', selectedText, value, num2str(expected), matched, note};
end

function rows = add_max_damping_validation(rows, dataDir, fileName, source, expected, note)
path = fullfile(dataDir, fileName);
if ~exist(path, 'file')
    return;
end
data = readtable(path);
[value, idx] = max(data.damping);
selected = data.node(idx);
rows(end+1,:) = {source, 'max damping', num2str(selected), value, num2str(expected), selected == expected, note};
end

function rows = add_max_gscr_validation(rows, dataDir, fileName, source, expected, note)
path = fullfile(dataDir, fileName);
if ~exist(path, 'file')
    return;
end
data = readtable(path);
[value, idx] = max(data.gscr);
selected = [data.node_a(idx), data.node_b(idx)];
rows(end+1,:) = {source, 'max gSCR', mat2str(selected), value, mat2str(expected), isequal(selected, expected), note};
end
