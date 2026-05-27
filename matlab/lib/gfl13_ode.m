function dxdt = gfl13_ode(x, p, X_net, S_ref, V_inf)
%GFL13_ODE  13-state GFL multi-converter ODE with full LCL filter dynamics.
%
%  Extends gfl11_ode by replacing the quasi-static algebraic network with
%  explicit grid-side inductor current states I_gd_g, I_gq_g (GLOBAL frame).
%  This eliminates the spurious high-frequency unstable modes that appear
%  when algebraic network inversion is combined with explicit Cf states.
%
%  State layout (13 states per converter, interleaved):
%    x(13k-12) = phi        – PLL angle deviation from global frame (rad)
%    x(13k-11) = eps_pll    – PLL integrator
%    x(13k-10) = eps_p      – power outer-loop integrator
%    x(13k-9 ) = eps_id     – d-axis current inner-loop integrator
%    x(13k-8 ) = eps_iq     – q-axis current inner-loop integrator
%    x(13k-7 ) = id         – d-axis converter current through Lf (local)
%    x(13k-6 ) = iq         – q-axis converter current through Lf (local)
%    x(13k-5 ) = Vd_ff      – GFF-filtered capacitor voltage, d-axis
%    x(13k-4 ) = Vq_ff      – GFF-filtered capacitor voltage, q-axis
%    x(13k-3 ) = vCd        – filter capacitor voltage, d-axis (local frame)
%    x(13k-2 ) = vCq        – filter capacitor voltage, q-axis (local frame)
%    x(13k-1 ) = I_gd_g     – grid-side d-axis current (GLOBAL frame)
%    x(13k   ) = I_gq_g     – grid-side q-axis current (GLOBAL frame)
%
%  Network dynamics (global rotating frame, X_net = omega0 * L_net):
%    dI_gd_g/dt = omega0*(X_net\(vCd_g - V_inf)) + omega0*I_gq_g
%    dI_gq_g/dt = omega0*(X_net\vCq_g)            - omega0*I_gd_g
%
%  At steady state these reduce to the quasi-static relations:
%    I_gd_g = X_net\vCq_g,   I_gq_g = X_net\(V_inf - vCd_g)

if nargin < 5 || isempty(V_inf), V_inf = 1.0; end

N   = numel(S_ref);
w0  = p.omega0;
Lf  = p.Lf;   Rf = p.Rf;   Cf = p.Cf;
tau = p.tau_ff;

% ── Unpack states ─────────────────────────────────────────────────────────
phi     = x(1:13:end);
eps_pll = x(2:13:end);
eps_p   = x(3:13:end);
eps_id  = x(4:13:end);
eps_iq  = x(5:13:end);
id_loc  = x(6:13:end);
iq_loc  = x(7:13:end);
Vd_ff   = x(8:13:end);
Vq_ff   = x(9:13:end);
vCd     = x(10:13:end);
vCq     = x(11:13:end);
I_gd_g  = x(12:13:end);
I_gq_g  = x(13:13:end);

% ── Capacitor voltage in global frame ─────────────────────────────────────
vCd_g = vCd.*cos(phi) - vCq.*sin(phi);
vCq_g = vCd.*sin(phi) + vCq.*cos(phi);

% ── Grid-side current in local dq frame (from STATE, not algebraic) ───────
i_gd = I_gd_g.*cos(phi) + I_gq_g.*sin(phi);
i_gq = -I_gd_g.*sin(phi) + I_gq_g.*cos(phi);

% ── Measured active power (PCC = capacitor bus) ───────────────────────────
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

% ── Lf inductor dynamics ─────────────────────────────────────────────────
did = (vd_inv - vCd - Rf*id_loc)./Lf + w0*iq_loc;
diq = (vq_inv - vCq - Rf*iq_loc)./Lf - w0*id_loc;

% ── GFF filter (tracks capacitor terminal voltage) ───────────────────────
dVd_ff = (vCd - Vd_ff) / tau;
dVq_ff = (vCq - Vq_ff) / tau;

% ── Cf capacitor dynamics (uses STATE grid current) ──────────────────────
dvCd = (id_loc - i_gd)./Cf + w0*vCq;
dvCq = (iq_loc - i_gq)./Cf - w0*vCd;

% ── Grid-side inductor dynamics (global rotating frame) ──────────────────
%  KVL: vCd_g - V_inf = (X_net/w0)*dI_gd_g/dt - X_net*I_gq_g
%       vCq_g         = (X_net/w0)*dI_gq_g/dt + X_net*I_gd_g
%  Rearranged (X_net = w0*L_net):
dI_gd_g = w0*(X_net \ (vCd_g - V_inf*ones(N,1))) + w0*I_gq_g;
dI_gq_g = w0*(X_net \ vCq_g)                      - w0*I_gd_g;

% ── Assemble dxdt ─────────────────────────────────────────────────────────
dxdt = zeros(13*N, 1);
dxdt(1:13:end)  = dphi;
dxdt(2:13:end)  = deps_pll;
dxdt(3:13:end)  = deps_p;
dxdt(4:13:end)  = deps_id;
dxdt(5:13:end)  = deps_iq;
dxdt(6:13:end)  = did;
dxdt(7:13:end)  = diq;
dxdt(8:13:end)  = dVd_ff;
dxdt(9:13:end)  = dVq_ff;
dxdt(10:13:end) = dvCd;
dxdt(11:13:end) = dvCq;
dxdt(12:13:end) = dI_gd_g;
dxdt(13:13:end) = dI_gq_g;
end
