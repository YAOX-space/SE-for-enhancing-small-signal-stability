function test_gfl13_n1_spectrum()
%TEST_GFL13_N1_SPECTRUM  Print top-5 rightmost eigenvalues for single-converter.
addpath(fullfile(fileparts(mfilename('fullpath')), 'lib'));
p = gfl_control_params();

for gscr_test = [2.57, 2.66, 3.00, 3.86]
    Z11 = 1/gscr_test;
    x_ss = gfl13_find_ss(Z11, 1.0, p);
    f0   = gfl13_ode(x_ss, p, Z11, 1.0, 1.0);
    J    = zeros(13);
    for j = 1:13
        xp = x_ss; xp(j) = xp(j) + 1e-7;
        J(:,j) = (gfl13_ode(xp, p, Z11, 1.0, 1.0) - f0) / 1e-7;
    end
    ev = eig(J);
    [~, idx] = sort(real(ev), 'descend');
    fprintf('\ngSCR=%.2f  top-5 rightmost:\n', gscr_test);
    for k = 1:5
        fprintf('  %+9.3f  %+10.3f j\n', real(ev(idx(k))), imag(ev(idx(k))));
    end
end
end
