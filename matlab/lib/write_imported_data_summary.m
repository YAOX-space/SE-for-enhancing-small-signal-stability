function write_imported_data_summary(data, resultsDir)
%WRITE_IMPORTED_DATA_SUMMARY Summarize imported tables and derived networks.

names = fieldnames(data);
rows = cell(numel(names), 4);
for k = 1:numel(names)
    rows{k,1} = names{k};
    rows{k,2} = height(data.(names{k}));
    rows{k,3} = width(data.(names{k}));
    rows{k,4} = strjoin(data.(names{k}).Properties.VariableNames, ',');
end
summary = cell2table(rows, 'VariableNames', {'table_name','rows','columns','variables'});
writetable(summary, fullfile(resultsDir, 'imported_data_summary.csv'));

if isfield(data, 'ref20_33converter_topology_inferred_from_fig7') && ...
        isfield(data, 'ref20_33converter_derived_reactances')
    topo = data.ref20_33converter_topology_inferred_from_fig7;
    rx = data.ref20_33converter_derived_reactances;
    rxElement = string(rx{:,1});
    rxValue = rx{:,2};
    if iscell(rxValue) || isstring(rxValue)
        rxValue = str2double(string(rxValue));
    end
    boost = rxValue(contains(rxElement, 'converter_boost_transformer'));
    collector = rxValue(contains(rxElement, 'collector_line_between_adjacent'));
    grid = rxValue(contains(rxElement, 'grid_transformer'));
    boost = boost(1);
    collector = collector(1);
    grid = grid(1);

    x = zeros(height(topo), 1);
    for k = 1:height(topo)
        element = string(topo{k,3});
        if contains(element, 'converter_boost_transformer')
            x(k) = boost;
        elseif contains(element, 'grid_transformer')
            x(k) = grid;
        else
            x(k) = collector;
        end
    end
    edgeTable = table(topo.from, topo.to, x, string(topo{:,3}), ...
        'VariableNames', {'from','to','x','element'});
    writetable(edgeTable, fullfile(resultsDir, 'ref20_33converter_imported_edges.csv'));

    [Q, Z] = build_network_from_edges(edgeTable(:, {'from','to','x'}), 68, 68);
    save(fullfile(resultsDir, 'ref20_33converter_network_matrices.mat'), 'Q', 'Z', 'edgeTable');

    networkSummary = table(height(edgeTable), 68, nnz(Q), cond(Q), ...
        'VariableNames', {'edge_count','node_count','grounded_laplacian_nonzeros','grounded_laplacian_condition'});
    writetable(networkSummary, fullfile(resultsDir, 'ref20_33converter_network_summary.csv'));
end
end
