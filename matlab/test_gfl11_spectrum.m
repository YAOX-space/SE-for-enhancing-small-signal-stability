function test_gfl11_spectrum()
addpath(fullfile(fileparts(mfilename('fullpath')), 'lib'));
p = gfl_control_params();

for gscr_test = [2.66, 3.666]
    Z11 = 1/gscr_test;
    x_ss = gfl11_find_ss(Z11, 1.0, p);
    f0   = gfl11_ode(x_ss, p, Z11, 1.0, 1.0);
    J    = zeros(11);
    for j = 1:11
        xp = x_ss; xp(j) = xp(j) + 1e-7;
        J(:,j) = (gfl11_ode(xp, p, Z11, 1.0, 1.0) - f0) / 1e-7;
    end
    ev = eig(J);
    fprintf('\n=== gSCR=%.3f  ALL eigenvalues ===\n', gscr_test);
    for k = 1:11
        fprintf('  %+8.2f  %+8.2f j\n', real(ev(k)), imag(ev(k)));
    end

    % Compare 9-state at same gSCR
    x9 = gfl9_find_ss(Z11, 1.0, p);
    f09 = gfl9_ode(x9, p, Z11, 1.0, 1.0);
    J9   = zeros(9);
    for j = 1:9
        xp = x9; xp(j) = xp(j) + 1e-7;
        J9(:,j) = (gfl9_ode(xp, p, Z11, 1.0, 1.0) - f09) / 1e-7;
    end
    ev9 = eig(J9);
    fprintf('--- gSCR=%.3f  9-state eigenvalues ---\n', gscr_test);
    for k = 1:9
        fprintf('  %+8.2f  %+8.2f j\n', real(ev9(k)), imag(ev9(k)));
    end
end
end
