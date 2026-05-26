function diagnose_feedforward_filter()
%DIAGNOSE_FEEDFORWARD_FILTER  Test whether the voltage-feedforward filter
%  (GFF = 1/(1+0.01s), from Table IX) shifts the dominant eigenvalue from
%  the -12 of the pure 7-state model toward the paper's +0.090.
%
%  9-state per converter:
%    [phi, eps_pll, eps_p, eps_id, eps_iq, id, iq, Vd_ff, Vq_ff]
%  New states: Vd_ff, Vq_ff = filtered terminal voltages fed to current ctrl.
%
%  Paper Table IX: GFF_voltage_feedforward_filter = 1/(1+0.01s)
%  => tau_ff = 0.01 s, omega_ff = 100 rad/s (16 Hz)
%  => at 88 rad/s: |G_ff| = 0.752, phase = -41.3 deg  (significant!)

close all;
rootDir = fileparts(mfilename('fullpath'));
addpath(fullfile(rootDir,'lib'));
dataDir = fullfile(rootDir,'data');

fprintf('=== GFL 9-state (+ GFF filter) Eigenvalue Diagnostic ===\n\n');

p.omega0  = 2*pi*50;
p.Lf      = 0.05 / p.omega0;
p.Rf      = 0.002 / p.omega0;
p.Kp_pll  = 26;    p.Ki_pll = 7800;
p.Kp_p    = 0.5;   p.Ki_p   = 5;
p.Kp_i    = 1;     p.Ki_i   = 10;
p.tau_ff  = 0.01;  % GFF filter time constant (Table IX)

% Verify filter phase at 88 rad/s
w88 = 88;
Gff88 = 1 / (1 + 1i*w88*p.tau_ff);
fprintf('GFF filter at 88 rad/s:  |G|=%.4f  phase=%.2f deg\n', ...
    abs(Gff88), angle(Gff88)*180/pi);
fprintf('Uncancelled fraction: %.4f\n\n', abs(1 - Gff88));

% ── IEEE 39-node network (same calibration as diagnose_eigenvalues.m) ─────
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

% ── Steady state via ODE integration ─────────────────────────────────────
N    = numel(S_CBR);
nSt  = 9*N;   % 9 states per converter

% Initial guess (from 7-state SS, extend with Vd_ff≈Vd_l, Vq_ff≈Vq_l)
id_ss0 = S_CBR(:);
phi0   = atan2(X_net * id_ss0, ones(N,1));
x0     = zeros(nSt,1);
for k = 1:N
    b = (k-1)*9;
    x0(b+1) = phi0(k);           % phi
    x0(b+4) = p.Rf * id_ss0(k) / p.Ki_i;  % eps_id
    x0(b+6) = id_ss0(k);         % id
    % Vd_ff, Vq_ff: set to approximate SS values
    phi_k   = phi0(k);
    Iq_g0   = zeros(N,1);  % first approx: iq=0
    Id_g0   = id_ss0;
    Vd_g0   = 1.0 - X_net * Iq_g0;
    Vq_g0   =       X_net * Id_g0;
    x0(b+8) = Vd_g0(k)*cos(phi_k) + Vq_g0(k)*sin(phi_k);  % Vd_l ≈ Vd_ff
    x0(b+9) = -Vd_g0(k)*sin(phi_k) + Vq_g0(k)*cos(phi_k); % Vq_l ≈ Vq_ff
end

odefun   = @(t,x) gfl9_rotating_ode(x, p, X_net, S_CBR);
opts_ode = odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',1e-3);
fprintf('Integrating 9-state model to SS (30 s)...\n');
[~, xt] = ode15s(odefun, [0 30], x0, opts_ode);
x_ss    = xt(end,:)';

res = gfl9_rotating_ode(x_ss, p, X_net, S_CBR);
fprintf('Max |residual| at SS: %.2e\n', max(abs(res)));
fprintf('SS phi (rad):    '); fprintf('%.4f ', x_ss(1:9:end)'); fprintf('\n');
fprintf('SS id:           '); fprintf('%.4f ', x_ss(6:9:end)'); fprintf('\n');
fprintf('SS Vd_ff - Vd_l: ');
% Compute Vd_l at SS to compare with Vd_ff
phi_ss  = x_ss(1:9:end);
id_ss   = x_ss(6:9:end);
iq_ss   = x_ss(7:9:end);
Id_g_ss = id_ss.*cos(phi_ss) - iq_ss.*sin(phi_ss);
Iq_g_ss = id_ss.*sin(phi_ss) + iq_ss.*cos(phi_ss);
Vd_g_ss = 1.0 - X_net*Iq_g_ss;
Vq_g_ss =       X_net*Id_g_ss;
Vd_l_ss =  Vd_g_ss.*cos(phi_ss) + Vq_g_ss.*sin(phi_ss);
Vd_ff_ss = x_ss(8:9:end);
fprintf('%.2e ', Vd_ff_ss - Vd_l_ss); fprintf('\n');

% ── Numerical Jacobian ────────────────────────────────────────────────────
fprintf('\nComputing 81×81 Jacobian...\n');
eps_fd = 1e-7;
f0 = gfl9_rotating_ode(x_ss, p, X_net, S_CBR);
J  = zeros(nSt, nSt);
for j = 1:nSt
    xp = x_ss;  xp(j) = xp(j) + eps_fd;
    J(:,j) = (gfl9_rotating_ode(xp, p, X_net, S_CBR) - f0) / eps_fd;
end

ev = eig(J);
[~,ix] = sort(real(ev),'descend');
ev = ev(ix);

fprintf('\nAll eigenvalues near omega=88 rad/s:\n');
mask = abs(abs(imag(ev)) - 88) < 10;
for k = find(mask)'
    fprintf('  %+8.4f %+8.4fi\n', real(ev(k)), imag(ev(k)));
end

fprintf('\nTop 6 eigenvalues (rightmost real part):\n');
for k = 1:6
    fprintf('  %+8.4f %+8.4fi\n', real(ev(k)), imag(ev(k)));
end

fprintf('\n── Paper Table V vs 9-state GFF Model ───────────────────────\n');
fprintf('  Paper gSCR=2.650:        sigma = +0.0900  omega = 88.342\n');
fprintf('  7-state (no GFF):        sigma = -12.035   omega = 84.142\n');
[~, bi] = min(abs(abs(imag(ev)) - 88.3));
lam = ev(bi);
fprintf('  9-state (with GFF):      sigma = %+.4f  omega = %.4f\n', ...
    real(lam), abs(imag(lam)));
fprintf('  Residual gap: %.4f\n\n', real(lam) - 0.090);

% ── Sweep tau_ff to find which value gives sigma≈0.09 ────────────────────
fprintf('── tau_ff sensitivity (sigma vs tau_ff) ──\n');
tau_vals = [0.0, 0.001, 0.002, 0.005, 0.010, 0.020, 0.050, 0.100];
for tau = tau_vals
    p2       = p;
    p2.tau_ff = tau;
    if tau == 0
        % tau=0: instantaneous feedforward → 7-state model behaviour
        fprintf('  tau_ff = 0 (no filter): ');
    else
        fprintf('  tau_ff = %.3f s: ', tau);
    end
    [~, xt2] = ode15s(@(t,x) gfl9_rotating_ode(x,p2,X_net,S_CBR), ...
                      [0 30], x0, opts_ode);
    x2  = xt2(end,:)';
    f2  = gfl9_rotating_ode(x2, p2, X_net, S_CBR);
    J2  = zeros(nSt, nSt);
    eps2 = 1e-7;
    for j = 1:nSt
        xp = x2; xp(j) = xp(j) + eps2;
        J2(:,j) = (gfl9_rotating_ode(xp,p2,X_net,S_CBR) - f2) / eps2;
    end
    ev2 = eig(J2);
    [~, bi2] = min(abs(abs(imag(ev2)) - 88.0));
    lam2 = ev2(bi2);
    fprintf('sigma = %+.4f  omega = %.3f\n', real(lam2), abs(imag(lam2)));
end

end % function


% ─────────────────────────────────────────────────────────────────────────
%  9-state GFL ODE in co-rotating frame  (adds GFF filter states)
%  State per converter: [phi, eps_pll, eps_p, eps_id, eps_iq, id, iq,
%                         Vd_ff, Vq_ff]
%  New:  tau_ff * dVd_ff/dt = Vd_l - Vd_ff  (voltage feedforward filter)
%        tau_ff * dVq_ff/dt = Vq_l - Vq_ff
%  Modified current ctrl:  vd_vsc uses Vd_ff instead of Vd_l
%                           vq_vsc uses Vq_ff instead of Vq_l
% ─────────────────────────────────────────────────────────────────────────
function dxdt = gfl9_rotating_ode(x, p, X_net, S_CBR)
N   = numel(S_CBR);
w0  = p.omega0;  Lf = p.Lf;  Rf = p.Rf;
tau = p.tau_ff;

phi     = x(1:9:end);
eps_pll = x(2:9:end);
eps_p   = x(3:9:end);
eps_id  = x(4:9:end);
eps_iq  = x(5:9:end);
id_loc  = x(6:9:end);
iq_loc  = x(7:9:end);
Vd_ff   = x(8:9:end);   % filtered Vd_l (local frame)
Vq_ff   = x(9:9:end);   % filtered Vq_l (local frame)

% Global network currents
Id_g = id_loc.*cos(phi) - iq_loc.*sin(phi);
Iq_g = id_loc.*sin(phi) + iq_loc.*cos(phi);

% Network voltages (global frame, quasi-static)
Vd_g = 1.0 - X_net * Iq_g;
Vq_g =       X_net * Id_g;

% Terminal voltage in local frame (actual, unfiltered)
Vd_l =  Vd_g.*cos(phi) + Vq_g.*sin(phi);
Vq_l = -Vd_g.*sin(phi) + Vq_g.*cos(phi);

% Power (uses actual voltages, not filtered)
P_meas = Vd_l.*id_loc + Vq_l.*iq_loc;
Pref   = S_CBR(:);

% Control
id_ref = p.Kp_p*(Pref - P_meas) + p.Ki_p*eps_p;
iq_ref = zeros(N,1);

% Current ctrl uses FILTERED voltages (GFF)
vd_vsc = p.Kp_i*(id_ref - id_loc) + p.Ki_i*eps_id - w0*Lf*iq_loc + Vd_ff;
vq_vsc = p.Kp_i*(iq_ref - iq_loc) + p.Ki_i*eps_iq + w0*Lf*id_loc + Vq_ff;

% ODEs
dphi     = p.Kp_pll*Vq_l + p.Ki_pll*eps_pll;
deps_pll = Vq_l;
deps_p   = Pref - P_meas;
deps_id  = id_ref - id_loc;
deps_iq  = iq_ref - iq_loc;
did      = (vd_vsc - Vd_l - Rf*id_loc)./Lf + w0*iq_loc;
diq      = (vq_vsc - Vq_l - Rf*iq_loc)./Lf - w0*id_loc;

% GFF filter dynamics (local-frame voltages)
if tau > 0
    dVd_ff = (Vd_l - Vd_ff) / tau;
    dVq_ff = (Vq_l - Vq_ff) / tau;
else
    % tau=0: instantaneous — Vd_ff = Vd_l always, no dynamics needed
    % (ODE system degenerates; keep dVd_ff≈0 by using large effective rate)
    dVd_ff = (Vd_l - Vd_ff) * 1e6;
    dVq_ff = (Vq_l - Vq_ff) * 1e6;
end

dxdt = zeros(9*N, 1);
dxdt(1:9:end) = dphi;
dxdt(2:9:end) = deps_pll;
dxdt(3:9:end) = deps_p;
dxdt(4:9:end) = deps_id;
dxdt(5:9:end) = deps_iq;
dxdt(6:9:end) = did;
dxdt(7:9:end) = diq;
dxdt(8:9:end) = dVd_ff;
dxdt(9:9:end) = dVq_ff;
end
