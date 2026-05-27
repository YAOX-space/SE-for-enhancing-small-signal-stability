function test_gfl13_quick()
addpath(fullfile(fileparts(mfilename('fullpath')), 'lib'));
p = gfl_control_params();
fprintf('Cf = %.6f  (expected %.6f)\n', p.Cf, 0.05/(2*pi*50));

fprintf('\n13-state model eigenvalue sweep:\n');
fprintf('%-8s  %-10s  %-10s\n', 'gSCR', 'sigma', 'omega');
for gscr_test = [2.57, 2.66, 3.00, 3.86]
    Z11 = 1/gscr_test;
    try
        x_ss = gfl13_find_ss(Z11, 1.0, p);
        res  = max(abs(gfl13_ode(x_ss, p, Z11, 1.0, 1.0)));
        f0   = gfl13_ode(x_ss, p, Z11, 1.0, 1.0);
        J    = zeros(13);
        for j = 1:13
            xp = x_ss; xp(j) = xp(j) + 1e-7;
            J(:,j) = (gfl13_ode(xp, p, Z11, 1.0, 1.0) - f0) / 1e-7;
        end
        ev = eig(J);
        % Find dominant PLL mode near 88 rad/s
        [~, bi] = min(abs(abs(imag(ev)) - 88));
        sigma = real(ev(bi)); omega = abs(imag(ev(bi)));
        fprintf('gSCR=%5.2f  sigma=%+7.4f  omega=%7.2f  SS_res=%.1e\n', ...
            gscr_test, sigma, omega, res);
    catch ME
        fprintf('gSCR=%5.2f  ERROR: %s\n', gscr_test, ME.message);
    end
end

fprintf('\n13-state full spectrum at gSCR=2.66:\n');
Z11 = 1/2.66;
x_ss = gfl13_find_ss(Z11, 1.0, p);
f0 = gfl13_ode(x_ss, p, Z11, 1.0, 1.0);
J = zeros(13);
for j = 1:13
    xp = x_ss; xp(j) = xp(j) + 1e-7;
    J(:,j) = (gfl13_ode(xp, p, Z11, 1.0, 1.0) - f0) / 1e-7;
end
ev = eig(J);
[~, idx] = sort(real(ev), 'descend');
for k = 1:13
    fprintf('  %+9.3f  %+10.3f j\n', real(ev(idx(k))), imag(ev(idx(k))));
end
end
