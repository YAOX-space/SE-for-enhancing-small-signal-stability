function [placement, trace] = greedy_storage_placement(capacities, impedance, candidates, storageCapacity, count)
%GREEDY_STORAGE_PLACEMENT Paper Fig. 4 placement algorithm.

capacities = capacities(:).';
candidateNodes = candidates(:).';
placement = zeros(1, count);
trace = repmat(struct( ...
    'iteration', 0, ...
    'lambda_max_before', 0, ...
    'gscr_before', 0, ...
    'chosen_node', 0, ...
    'min_sensitivity', 0), count, 1);

for iter = 1:count
    weighted = diag(capacities) * impedance;
    lambda = max_positive_eig(weighted);
    [V,D] = eig(weighted.');
    [~, idx] = min(abs(diag(D) - lambda));
    v = real(V(:,idx));
    v = v ./ norm(v);
    sensitivities = -lambda * (v(candidateNodes) .^ 2);
    [minSens, minIdx] = min(sensitivities);
    chosen = candidateNodes(minIdx);

    placement(iter) = chosen;
    trace(iter).iteration = iter;
    trace(iter).lambda_max_before = lambda;
    trace(iter).gscr_before = 1 / lambda;
    trace(iter).chosen_node = chosen;
    trace(iter).min_sensitivity = minSens;

    capacities(chosen) = capacities(chosen) - storageCapacity;
end
end
