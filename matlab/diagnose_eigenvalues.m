function diagnose_eigenvalues()
%DIAGNOSE_EIGENVALUES  Numerically compute the dominant eigenvalue of the
%  full 7-state GFL system at the SS operating point and compare with Table V.
%
%  Key: work in the CO-ROTATING FRAME  φ_k = θ_k - ω₀·t
%  so the steady state is a TRUE FIXED POINT (all derivatives = 0).

close all;
rootDir = fileparts(mfilename('fullpath'));
addpath(fullfile(rootDir,'lib'));
dataDir = fullfile(rootDir,'data');

fprintf('=== GFL ODE Eigenvalue Diagnostic (co-rotating frame) ===\n\n');

p.omega0  = 2*pi*50;
p.Lf      = 0.05 / p.omega0;
p.Rf      = 0.002 / p.omega0;
p.Kp_pll  = 26;   p.Ki_pll = 7800;
p.Kp_p    = 0.5;  p.Ki_p   = 5;
p.Kp_i    = 1;    p.Ki_i   = 10;

% IEEE 39-node, 9 CBR buses, calibrated gSCR = 2.650
lines39  = readtable(fullfile(dataDir,'table_xi_ieee39_lines.csv'));
nbus39   = max([lines39.from; lines39.to]);
[~, Z39] = build_network_from_edges(lines39, nbus39, 36);
keep39   = setdiff(1:nbus39, 36);
b2r39    = zeros(1,nbus39);  b2r39(keep39) = 1:38;
cbr_rows = b2r39([30 31 32 33 34 35 37 38 39]);
Z_raw    = Z39(cbr_rows, cbr_rows);
k_cal    = (1/2.650) / max_positive_eig(Z_raw);
X_net    = k_cal * Z_raw;
S_CBR    = ones(9,1);

fprintf('gSCR = %.4f\n', 1/max_positive_eig(diag(S_CBR)*X_net));
fprintf('X_net diag: '); fprintf('%.4f ', diag(X_net)'); fprintf('\n');

% ── Integrate co-rotating ODE to find steady state ───────────────────────
% State: [phi, eps_pll, eps_p, eps_id, eps_iq, id, iq] × N
% phi = theta - omega0*t  →  at SS: phi = const, d(phi)/dt = 0
N   = numel(S_CBR);
nSt = 7*N;

% Initial guess
id_ss0 = S_CBR(:);
phi0   = atan2(X_net * id_ss0, ones(N,1));
x0     = zeros(nSt,1);
for k = 1:N
    b = (k-1)*7;
    x0(b+1) = phi0(k);
    x0(b+4) = p.Rf * id_ss0(k) / p.Ki_i;
    x0(b+6) = id_ss0(k);
end

odefun = @(t,x) gfl_rotating_ode(x, p, X_net, S_CBR);
opts_ode = odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',1e-3);
fprintf('Integrating to steady state...\n');
[~, xt] = ode15s(odefun, [0 20], x0, opts_ode);
x_ss = xt(end,:)';

fprintf('SS phi_k (rad): ');   fprintf('%.4f ', x_ss(1:7:end)'); fprintf('\n');
fprintf('SS id_k:        ');   fprintf('%.4f ', x_ss(6:7:end)'); fprintf('\n');
fprintf('SS iq_k:        ');   fprintf('%.5f ', x_ss(7:7:end)'); fprintf('\n');

% Verify residual
res = gfl_rotating_ode(x_ss, p, X_net, S_CBR);
fprintf('Max |residual| at SS: %.2e\n', max(abs(res)));

% ── Numerical Jacobian ────────────────────────────────────────────────────
fprintf('\nComputing Jacobian...\n');
eps_fd = 1e-7;
f0 = gfl_rotating_ode(x_ss, p, X_net, S_CBR);
J  = zeros(nSt, nSt);
for j = 1:nSt
    xp = x_ss;  xp(j) = xp(j) + eps_fd;
    J(:,j) = (gfl_rotating_ode(xp, p, X_net, S_CBR) - f0) / eps_fd;
end

ev = eig(J);
[~,ix] = sort(real(ev),'descend');
ev = ev(ix);

fprintf('\nAll eigenvalues near omega=88 rad/s:\n');
mask = abs(abs(imag(ev)) - 88) < 5;
for k = find(mask)'
    fprintf('  %+8.4f %+8.4fi\n', real(ev(k)), imag(ev(k)));
end

fprintf('\nTop 6 eigenvalues (rightmost real part):\n');
for k = 1:6
    fprintf('  %+8.4f %+8.4fi\n', real(ev(k)), imag(ev(k)));
end

fprintf('\n── Paper Table V vs Our Model ──────────────────────────\n');
fprintf('  Paper gSCR=2.650:  sigma = +0.0900  omega = 88.342\n');
[~, bi] = min(abs(abs(imag(ev)) - 88.342));
lam = ev(bi);
fprintf('  Our model:         sigma = %+.4f  omega = %.4f\n', real(lam), abs(imag(lam)));
fprintf('  Gap: %.4f  (>0 means our model is MORE stable than paper)\n\n', real(lam) - 0.090);

% ── Sensitivity: vary Lf from 0.001 to 2.0 pu ───────────────────────────
fprintf('── Lf sensitivity (dominant eigenvalue sigma vs Lf_pu) ──\n');
Lf_pu_vals = [0.001, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 1.0];
for Lf_pu = Lf_pu_vals
    p2 = p;
    p2.Lf = Lf_pu / p.omega0;
    p2.Rf = Lf_pu * 0.04 / p.omega0;

    % Re-solve SS
    [~, xt2] = ode15s(@(t,x) gfl_rotating_ode(x,p2,X_net,S_CBR), [0 20], x0, opts_ode);
    x2 = xt2(end,:)';
    f2 = gfl_rotating_ode(x2, p2, X_net, S_CBR);
    J2 = zeros(nSt, nSt);
    for j = 1:nSt
        xp = x2; xp(j) = xp(j) + eps_fd;
        J2(:,j) = (gfl_rotating_ode(xp, p2, X_net, S_CBR) - f2) / eps_fd;
    end
    ev2 = eig(J2);
    [~, bi2] = min(abs(abs(imag(ev2)) - 88.0));
    lam2 = ev2(bi2);
    fprintf('  Lf = %.3f pu:  sigma = %+.4f  omega = %.3f\n', Lf_pu, real(lam2), abs(imag(lam2)));
end
end


% ─────────────────────────────────────────────────────────────────────────
%  7-state GFL ODE in the CO-ROTATING FRAME
%  State: [phi_k, eps_pll_k, eps_p_k, eps_id_k, eps_iq_k, id_k, iq_k]
%  phi_k = theta_k - omega0*t  →  d(phi_k)/dt = dtheta_k/dt - omega0
% ─────────────────────────────────────────────────────────────────────────
function dxdt = gfl_rotating_ode(x, p, X_net, S_CBR)
N  = numel(S_CBR);
w0 = p.omega0;  Lf = p.Lf;  Rf = p.Rf;

phi     = x(1:7:end);   % phi = theta - w0*t
eps_pll = x(2:7:end);
eps_p   = x(3:7:end);
eps_id  = x(4:7:end);
eps_iq  = x(5:7:end);
id_loc  = x(6:7:end);
iq_loc  = x(7:7:end);

% Network equations use phi (same as theta in rotating frame context)
Id_g = id_loc.*cos(phi) - iq_loc.*sin(phi);
Iq_g = id_loc.*sin(phi) + iq_loc.*cos(phi);

Vd_g = 1.0 - X_net * Iq_g;
Vq_g =       X_net * Id_g;

Vd_l =  Vd_g.*cos(phi) + Vq_g.*sin(phi);
Vq_l = -Vd_g.*sin(phi) + Vq_g.*cos(phi);

P_meas = Vd_l.*id_loc + Vq_l.*iq_loc;
Pref   = S_CBR(:);

id_ref = p.Kp_p*(Pref - P_meas) + p.Ki_p*eps_p;
iq_ref = zeros(N,1);

vd_vsc = p.Kp_i*(id_ref - id_loc) + p.Ki_i*eps_id - w0*Lf*iq_loc + Vd_l;
vq_vsc = p.Kp_i*(iq_ref - iq_loc) + p.Ki_i*eps_iq + w0*Lf*id_loc + Vq_l;

% Co-rotating frame: dphi/dt = dtheta/dt - w0
dphi     = p.Kp_pll*Vq_l + p.Ki_pll*eps_pll;   % = dtheta/dt - w0
deps_pll = Vq_l;
deps_p   = Pref - P_meas;
deps_id  = id_ref - id_loc;
deps_iq  = iq_ref - iq_loc;
did      = (vd_vsc - Vd_l - Rf*id_loc)./Lf + w0*iq_loc;
diq      = (vq_vsc - Vq_l - Rf*iq_loc)./Lf - w0*id_loc;

dxdt = zeros(7*N,1);
dxdt(1:7:end) = dphi;
dxdt(2:7:end) = deps_pll;
dxdt(3:7:end) = deps_p;
dxdt(4:7:end) = deps_id;
dxdt(5:7:end) = deps_iq;
dxdt(6:7:end) = did;
dxdt(7:7:end) = diq;
end
