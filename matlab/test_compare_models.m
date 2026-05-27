function test_compare_models()
%TEST_COMPARE_MODELS  Compare sigma from gfl9 and gfl13 for N=1 at various gSCR.
addpath(fullfile(fileparts(mfilename('fullpath')), 'lib'));
p = gfl_control_params();

fprintf('%-8s  %-14s  %-14s\n', 'gSCR', 'gfl9_sigma', 'gfl13_sigma');
fprintf('%s\n', repmat('-', 1, 40));
for gscr = [2.0, 2.4, 2.57, 2.66, 2.82, 3.0, 3.66, 5.0]
    Z11 = 1/gscr;

    % --- gfl9 ---
    try
        x9  = gfl9_find_ss(Z11, 1.0, p);
        f9  = gfl9_ode(x9, p, Z11, 1.0, 1.0);
        J9  = zeros(9);
        for j = 1:9
            xp = x9; xp(j) = xp(j)+1e-7;
            J9(:,j) = (gfl9_ode(xp, p, Z11, 1.0, 1.0) - f9) / 1e-7;
        end
        ev9 = eig(J9);
        mask9 = abs(imag(ev9)) >= 70 & abs(imag(ev9)) <= 130;
        ev9p = ev9(mask9);
        if isempty(ev9p)
            [~,bi] = min(abs(abs(imag(ev9))-88)); lam9 = ev9(bi);
        else
            [~,bi] = max(real(ev9p)); lam9 = ev9p(bi);
        end
        s9 = real(lam9);
    catch
        s9 = NaN;
    end

    % --- gfl13 ---
    try
        x13  = gfl13_find_ss(Z11, 1.0, p);
        f13  = gfl13_ode(x13, p, Z11, 1.0, 1.0);
        J13  = zeros(13);
        for j = 1:13
            xp = x13; xp(j) = xp(j)+1e-7;
            J13(:,j) = (gfl13_ode(xp, p, Z11, 1.0, 1.0) - f13) / 1e-7;
        end
        ev13 = eig(J13);
        mask13 = abs(imag(ev13)) >= 70 & abs(imag(ev13)) <= 130;
        ev13p = ev13(mask13);
        if isempty(ev13p)
            [~,bi] = min(abs(abs(imag(ev13))-88)); lam13 = ev13(bi);
        else
            [~,bi] = max(real(ev13p)); lam13 = ev13p(bi);
        end
        s13 = real(lam13);
    catch
        s13 = NaN;
    end

    fprintf('%-8.2f  %+14.4f  %+14.4f\n', gscr, s9, s13);
end
end
