function run_all_reproduction()
%RUN_ALL_REPRODUCTION Reproduce the paper's numerical experiments.
%
% This package reproduces the paper at two levels:
%   1) gSCR/eigenvalue/storage-placement calculations from the equations.
%   2) Calibrated reduced-order time-domain validation of reported cases.
%
% The time-domain response model uses the weakest PLL-induced mode plus a
% smooth fault-period power sag. It is intended to reproduce stability trends
% and reported gSCR/damping comparisons, not the authors' unpublished EMT
% implementation.

close all;
rootDir = fileparts(mfilename('fullpath'));
addpath(fullfile(rootDir, 'lib'));
resultsDir = fullfile(rootDir, 'results');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

fprintf('Running small-signal and placement reproductions...\n');
paper = paper_case_data();
dataDir = fullfile(rootDir, 'data');
data = import_reproduction_data(dataDir);
write_imported_data_summary(data, resultsDir);
run_small_signal_experiments(paper, resultsDir);

fprintf('Running calibrated reduced-order time-domain reproductions...\n');
run_time_domain_experiments(paper, resultsDir);

fprintf('\nAll reproduction artifacts written to:\n  %s\n', resultsDir);
end
