function test_gfl13_spectrum9()
%TEST_GFL13_SPECTRUM9  Print eigenvalues near the PLL band for N=9 system.
addpath(fullfile(fileparts(mfilename('fullpath')), 'lib'));

p = gfl_control_params();
rootDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(rootDir, 'data');

lines39 = readtable(fullfile(dataDir, 'table_xi_ieee39_lines.csv'));
nbus39  = max([lines39.from; lines39.to]);
[~, Z39] = build_network_from_edges(lines39, nbus39, 36);
keep39   = setdiff(1:nbus39, 36);
b2r39    = zeros(1, nbus39);  b2r39(keep39) = 1:38;
cbr_rows = b2r39([30 31 32 33 34 35 37 38 39]);
Z_raw    = Z39(cbr_rows, cbr_rows);
S_CBR    = ones(9, 1);

k_base = (1/2.650) / max_positive_eig(diag(S_CBR) * Z_raw);
Z_base = k_base * Z_raw;

for g_tgt = [2.650, 2.926]
    alpha = 2.650 / g_tgt;
    X_net = alpha * Z_base;
    N = 9;

    fprintf('\n=== gSCR=%.4f  finding SS... ===\n', g_tgt);
    x_ss = gfl13_find_ss(X_net, S_CBR, p);
    f0   = gfl13_ode(x_ss, p, X_net, S_CBR, 1.0);
    fprintf('  Max |residual|: %.2e\n', max(abs(f0)));

    nSt = 13*N;
    J   = zeros(nSt);
    eps_fd = 1e-7;
    for j = 1:nSt
        xp = x_ss; xp(j) = xp(j) + eps_fd;
        J(:,j) = (gfl13_ode(xp, p, X_net, S_CBR, 1.0) - f0) / eps_fd;
    end
    ev = eig(J);

    % Print top-10 rightmost eigenvalues
    [~, idx] = sort(real(ev), 'descend');
    fprintf('  Top-10 rightmost eigenvalues:\n');
    for k = 1:10
        fprintf('    %+9.4f  %+10.4f j\n', real(ev(idx(k))), imag(ev(idx(k))));
    end

    % Print eigenvalues with |imag| in [70, 110]
    mask = abs(imag(ev)) >= 70 & abs(imag(ev)) <= 110;
    ev_band = ev(mask);
    [~, si] = sort(real(ev_band), 'descend');
    fprintf('  Eigenvalues with |omega| in [70,110] rad/s:\n');
    if isempty(ev_band)
        fprintf('    (none)\n');
    else
        for k = 1:numel(ev_band)
            fprintf('    %+9.4f  %+10.4f j\n', real(ev_band(si(k))), imag(ev_band(si(k))));
        end
    end
end
end
