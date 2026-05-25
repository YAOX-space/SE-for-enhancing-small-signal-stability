function [Q, Z] = build_network_from_edges(edgeTable, nodeCount, groundNode)
%BUILD_NETWORK_FROM_EDGES Build grounded Laplacian Q and impedance Z.

Qfull = zeros(nodeCount, nodeCount);
for k = 1:height(edgeTable)
    i = edgeTable.from(k);
    j = edgeTable.to(k);
    x = edgeTable.x(k);
    b = 1 / x;
    Qfull(i,i) = Qfull(i,i) + b;
    Qfull(j,j) = Qfull(j,j) + b;
    Qfull(i,j) = Qfull(i,j) - b;
    Qfull(j,i) = Qfull(j,i) - b;
end

keep = setdiff(1:nodeCount, groundNode);
Q = Qfull(keep, keep);
Z = inv(Q);
end
