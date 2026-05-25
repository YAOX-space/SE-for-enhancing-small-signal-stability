function run_simulink_reproduction()
%RUN_SIMULINK_REPRODUCTION Build and run the dedicated Simulink approximation.

close all;
rootDir = fileparts(mfilename('fullpath'));
addpath(fullfile(rootDir, 'lib'));
resultsDir = fullfile(rootDir, 'results');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

paper = paper_case_data();
run_simulink_time_domain_experiments(paper, resultsDir);

fprintf('\nSimulink model and outputs written to:\n');
fprintf('  %s\n', fullfile(rootDir, 'se_gfl_pll_small_signal_model.slx'));
fprintf('  %s\n', fullfile(resultsDir, 'simulink_fig6_two_area_time_domain.png'));
fprintf('  %s\n', fullfile(resultsDir, 'simulink_fig11_39node_gscr_time_domain.png'));
fprintf('  %s\n', fullfile(resultsDir, 'simulink_fig12_39node_placement_time_domain.png'));
fprintf('  %s\n', fullfile(resultsDir, 'simulink_fig14_33converter_time_domain.png'));
end
