function Z = build_network_z_pinv(edgeTable, n_nodes)
%BUILD_NETWORK_Z_PINV  Z-bus via Moore-Penrose pseudo-inverse of full Laplacian.
%
%  Unlike build_network_from_edges, this does NOT remove a reference node.
%  All inter-node coupling (including paths through the reference bus) is
%  preserved.  Required for the two-area system where the only inter-area
%  path runs through node 9 (infinite bus): removing node 9 severs it.
%
%  edgeTable  – table with columns: from (int), to (int), x (reactance)
%  n_nodes    – total node count including the reference/ground node
%
%  Z = pinv(Q_full),  where Q_full is the n_nodes×n_nodes susceptance Laplacian.

Q = zeros(n_nodes, n_nodes);
for k = 1:height(edgeTable)
    i = edgeTable.from(k);
    j = edgeTable.to(k);
    b = 1 / edgeTable.x(k);
    Q(i,i) = Q(i,i) + b;
    Q(j,j) = Q(j,j) + b;
    Q(i,j) = Q(i,j) - b;
    Q(j,i) = Q(j,i) - b;
end
Z = pinv(Q);
end
