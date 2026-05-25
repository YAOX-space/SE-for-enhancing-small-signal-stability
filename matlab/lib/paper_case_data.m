function paper = paper_case_data()
%PAPER_CASE_DATA Reported values and reconstructed matrices from the paper.

paper.CgSCR = 2.66;

paper.twoArea.weighted = [
    0.223 0.069 0.015 0.015 0.069 0.030 0.015 0.015 0.015
    0.139 0.238 0.030 0.030 0.139 0.060 0.030 0.030 0.030
    0.045 0.045 0.249 0.130 0.045 0.045 0.045 0.071 0.130
    0.015 0.015 0.043 0.246 0.015 0.015 0.015 0.024 0.043
    0     0     0     0     0     0     0     0     0
    0     0     0     0     0     0     0     0     0
    0     0     0     0     0     0     0     0     0
    0     0     0     0     0     0     0     0     0
    0     0     0     0     0     0     0     0     0
];
paper.twoArea.capacities = [0.5 1.0 1.5 0.5 0 0 0 0 0];
paper.twoArea.storageCapacity = 0.5;
paper.twoArea.storageCount = 2;
paper.twoArea.reportedPlacement = [1 4];
paper.twoArea.comparisonPlacement = [4 6];

paper.fig11.gscr = [3.856 2.659 2.571];
paper.fig11.labels = {'gSCR = 3.86', 'gSCR = 2.66', 'gSCR = 2.57'};

paper.fig12.cases = {
    'Case 1: no SE',              2.650, -0.001, []
    'Case 2: SEs at node 23',     2.801,  0.011, [23 23 23]
    'Case 3: SEs at 3,22,35',     3.043,  0.027, [3 22 35]
    'Case 4: SEs at 37,38,39',    3.666,  0.057, [37 38 39]
};

paper.ieee39.optimalPlacement = [37 38 39];
paper.ieee39.reportedOtherEig = {
    'SEs all at node 23', '-0.976 +/- 88.714i', 0.011
    'SEs at nodes 3,22,35', '-2.427 +/- 89.122i', 0.027
};

paper.fig14.cases = {
    'Case 1: no SE',             2.823, []
    'Case 2: SEs at 29-39',      3.059, 29:39
    'Case 3: proposed nodes',    3.308, [42 43 44 45 53 54 55 56 65 66 67]
};
paper.conv33.optimalPlacement = [42 43 44 45 53 54 55 56 65 66 67];

paper.controlTable.note = [
    "Appendix Table IX is transcribed into matlab/data/table_ix_control_parameters.csv. " + ...
    "The calibrated reduced-order reproduction still fits the weakest mode " + ...
    "and fault-period power sag to the reported CgSCR, gSCR and damping " + ...
    "ratios because the authors' EMT controller implementation is unpublished."
];
end
