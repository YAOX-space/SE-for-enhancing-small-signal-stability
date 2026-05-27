function dxdt = gfl11_ode(x, p, X_net, S_ref, V_inf)
%GFL11_ODE  11-state GFL multi-converter ODE with LCL filter capacitor.
%
%  Extends gfl9_ode by adding filter capacitor states vCd, vCq per
%  converter.  The quasi-static network now connects from the capacitor bus
%  (not the converter terminals) to the infinite bus.  This correctly
%  captures the LCL capacitor dynamics that determine PLL-mode stability.
%
%  State layout (11 states per converter, interleaved):
%    x(11k-10) = phi        – PLL angle deviation from global frame (rad)
%    x(11k-9 ) = eps_pll    – PLL integrator
%    x(11k-8 ) = eps_p      – power outer-loop integrator
%    x(11k-7 ) = eps_id     – d-axis current inner-loop integrator
%    x(11k-6 ) = eps_iq     – q-axis current inner-loop integrator
%    x(11k-5 ) = id         – d-axis converter current through Lf
%    x(11k-4 ) = iq         – q-axis converter current through Lf
%    x(11k-3 ) = Vd_ff      – GFF-filtered terminal voltage, d-axis
%    x(11k-2 ) = Vq_ff      – GFF-filtered terminal voltage, q-axis
%    x(11k-1 ) = vCd        – filter capacitor voltage, d-axis (local frame)
%    x(11k   ) = vCq        – filter capacitor voltage, q-axis (local frame)
%
%  Network model (quasi-static, inductive):
%    vCd_g(k) = V_inf - sum_j X_net(k,j)*Iq_g(j)   [d-axis, global]
%    vCq_g(k) =          sum_j X_net(k,j)*Id_g(j)   [q-axis, global]
%  Inverted algebraically to give grid-side injection currents:
%    Id_g = X_net \ vCq_g
%    Iq_g = X_net \ (V_inf*1 - vCd_g)

if nargin < 5 || isempty(V_inf), V_inf = 1.0; end

N   = numel(S_ref);
w0  = p.omega0;
Lf  = p.Lf;   Rf = p.Rf;   Cf = p.Cf;
tau = p.tau_ff;

% ── Unpack states ─────────────────────────────────────────────────────────
phi     = x(1:11:end);
eps_pll = x(2:11:end);
eps_p   = x(3:11:end);
eps_id  = x(4:11:end);
eps_iq  = x(5:11:end);
id_loc  = x(6:11:end);
iq_loc  = x(7:11:end);
Vd_ff   = x(8:11:end);
Vq_ff   = x(9:11:end);
vCd     = x(10:11:end);
vCq     = x(11:11:end);

% ── Capacitor voltage in global frame ─────────────────────────────────────
vCd_g = vCd.*cos(phi) - vCq.*sin(phi);
vCq_g = vCd.*sin(phi) + vCq.*cos(phi);

% ── Algebraic grid-side currents from quasi-static network ────────────────
%    Iq_g = X_net \ (V_inf - vCd_g)
%    Id_g = X_net \ vCq_g
Iq_g_gl = X_net \ (V_inf*ones(N,1) - vCd_g);
Id_g_gl = X_net \ vCq_g;

% ── Grid current in local dq frame ────────────────────────────────────────
i_gd = Id_g_gl.*cos(phi) + Iq_g_gl.*sin(phi);
i_gq = -Id_g_gl.*sin(phi) + Iq_g_gl.*cos(phi);

% ── Measured active power (PCC = capacitor bus, grid-side current) ────────
P_meas = vCd.*i_gd + vCq.*i_gq;
Pref   = S_ref(:);

% ── Outer power control loop → id reference ──────────────────────────────
id_ref = p.Kp_p*(Pref - P_meas) + p.Ki_p*eps_p;
iq_ref = zeros(N, 1);

% ── Inner current control (GFF feedforward from capacitor voltage) ────────
vd_inv = p.Kp_i*(id_ref - id_loc) + p.Ki_i*eps_id - w0*Lf*iq_loc + Vd_ff;
vq_inv = p.Kp_i*(iq_ref - iq_loc) + p.Ki_i*eps_iq + w0*Lf*id_loc + Vq_ff;

% ── PLL (locks q-axis of capacitor terminal voltage to zero) ─────────────
dphi     = p.Kp_pll*vCq + p.Ki_pll*eps_pll;
deps_pll = vCq;

% ── Power and current integrators ────────────────────────────────────────
deps_p  = Pref - P_meas;
deps_id = id_ref - id_loc;
deps_iq = iq_ref - iq_loc;

% ── Lf inductor dynamics (terminal voltage = vC) ─────────────────────────
did = (vd_inv - vCd - Rf*id_loc)./Lf + w0*iq_loc;
diq = (vq_inv - vCq - Rf*iq_loc)./Lf - w0*id_loc;

% ── GFF filter (tracks capacitor terminal voltage) ───────────────────────
dVd_ff = (vCd - Vd_ff) / tau;
dVq_ff = (vCq - Vq_ff) / tau;

% ── Cf capacitor dynamics ─────────────────────────────────────────────────
%    Cf*dvCd/dt = (id - i_gd) + w0*Cf*vCq  →  dvCd/dt = (id-i_gd)/Cf + w0*vCq
dvCd = (id_loc - i_gd)./Cf + w0*vCq;
dvCq = (iq_loc - i_gq)./Cf - w0*vCd;

% ── Assemble dxdt ─────────────────────────────────────────────────────────
dxdt = zeros(11*N, 1);
dxdt(1:11:end)  = dphi;
dxdt(2:11:end)  = deps_pll;
dxdt(3:11:end)  = deps_p;
dxdt(4:11:end)  = deps_id;
dxdt(5:11:end)  = deps_iq;
dxdt(6:11:end)  = did;
dxdt(7:11:end)  = diq;
dxdt(8:11:end)  = dVd_ff;
dxdt(9:11:end)  = dVq_ff;
dxdt(10:11:end) = dvCd;
dxdt(11:11:end) = dvCq;
end
