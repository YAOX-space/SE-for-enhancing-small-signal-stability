function test_gfl11_quick()
addpath(fullfile(fileparts(mfilename('fullpath')), 'lib'));
p = gfl_control_params();
fprintf('Cf = %.6f  (expected %.6f)\n', p.Cf, 0.05/(2*pi*50));

for gscr_test = [2.57, 2.66, 3.00, 3.86]
    Z11 = 1/gscr_test;
    ok = false; sigma = NaN; omega = NaN;
    try
        x_ss = gfl11_find_ss(Z11, 1.0, p);
        res  = max(abs(gfl11_ode(x_ss, p, Z11, 1.0, 1.0)));
        f0   = gfl11_ode(x_ss, p, Z11, 1.0, 1.0);
        J    = zeros(11);
        for j = 1:11
            xp = x_ss; xp(j) = xp(j) + 1e-7;
            J(:,j) = (gfl11_ode(xp, p, Z11, 1.0, 1.0) - f0) / 1e-7;
        end
        ev = eig(J);
        [~, bi] = min(abs(abs(imag(ev)) - 88));
        sigma = real(ev(bi)); omega = abs(imag(ev(bi))); ok = true;
        fprintf('gSCR=%5.2f  sigma=%+7.4f  omega=%7.2f  SS_res=%.1e\n', ...
            gscr_test, sigma, omega, res);
    catch ME
        fprintf('gSCR=%5.2f  ERROR: %s\n', gscr_test, ME.message);
    end
end
end
