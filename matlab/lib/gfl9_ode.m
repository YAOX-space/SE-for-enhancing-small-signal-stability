function dxdt = gfl9_ode(x, p, X_net, S_ref, V_inf)
%GFL9_ODE  9-state GFL multi-converter ODE in co-rotating frame.
%
%  Identical physics to gfl9_rotating_ode inside
%  compute_first_principles_eigenvalues.m, but extracted as a standalone
%  lib function and extended with the V_inf parameter for voltage-sag
%  injection.
%
%  State layout (9 states per converter, interleaved):
%    x( 9k-8 ) = phi_k       – PLL angle deviation from global frame (rad)
%    x( 9k-7 ) = eps_pll_k   – PLL integrator
%    x( 9k-6 ) = eps_p_k     – power outer-loop integrator
%    x( 9k-5 ) = eps_id_k    – d-axis current inner-loop integrator
%    x( 9k-4 ) = eps_iq_k    – q-axis current inner-loop integrator
%    x( 9k-3 ) = id_k        – d-axis current (converter local dq frame)
%    x( 9k-2 ) = iq_k        – q-axis current (converter local dq frame)
%    x( 9k-1 ) = Vd_ff_k     – GFF-filtered terminal voltage, d-axis
%    x( 9k   ) = Vq_ff_k     – GFF-filtered terminal voltage, q-axis
%
%  Inputs
%    x      – state vector  (9N × 1)
%    p      – control params struct from gfl_control_params()
%    X_net  – N×N reactance (impedance) matrix  [quasi-static network]
%    S_ref  – N×1 reference powers (signed: >0 CBR, <0 SE)
%    V_inf  – scalar, infinite-bus d-axis voltage (default 1.0 pu)
%             Set to 1 - delta_v(t) to inject a voltage sag.
%
%  Network model (quasi-static, purely inductive):
%    Vd_g(k) = V_inf  - sum_j X_net(k,j) * Iq_g(j)
%    Vq_g(k) =          sum_j X_net(k,j) * Id_g(j)
%  where Id_g, Iq_g are converter currents in the global rotating frame.

if nargin < 5 || isempty(V_inf)
    V_inf = 1.0;
end

N   = numel(S_ref);
w0  = p.omega0;
Lf  = p.Lf;   Rf = p.Rf;
tau = p.tau_ff;

% ── Unpack states ────────────────────────────────────────────────────────
phi     = x(1:9:end);
eps_pll = x(2:9:end);
eps_p   = x(3:9:end);
eps_id  = x(4:9:end);
eps_iq  = x(5:9:end);
id_loc  = x(6:9:end);
iq_loc  = x(7:9:end);
Vd_ff   = x(8:9:end);
Vq_ff   = x(9:9:end);

% ── Global-frame currents (dq rotation by phi_k) ─────────────────────────
Id_g = id_loc.*cos(phi) - iq_loc.*sin(phi);
Iq_g = id_loc.*sin(phi) + iq_loc.*cos(phi);

% ── Network voltages (global dq frame) ───────────────────────────────────
Vd_g = V_inf - X_net * Iq_g;
Vq_g =         X_net * Id_g;

% ── Terminal voltage in converter local dq frame (rotate by -phi_k) ──────
Vd_l =  Vd_g.*cos(phi) + Vq_g.*sin(phi);
Vq_l = -Vd_g.*sin(phi) + Vq_g.*cos(phi);

% ── Measured active power ─────────────────────────────────────────────────
P_meas = Vd_l.*id_loc + Vq_l.*iq_loc;
Pref   = S_ref(:);

% ── Outer power control loop → id reference ──────────────────────────────
id_ref = p.Kp_p*(Pref - P_meas) + p.Ki_p*eps_p;
iq_ref = zeros(N, 1);

% ── Inner current control (GFF-filtered voltages as feedforward) ──────────
% GFF at 88 rad/s: |G| = 0.752, phase = -41 deg → ~34% coupling residual
vd_vsc = p.Kp_i*(id_ref - id_loc) + p.Ki_i*eps_id - w0*Lf*iq_loc + Vd_ff;
vq_vsc = p.Kp_i*(iq_ref - iq_loc) + p.Ki_i*eps_iq + w0*Lf*id_loc + Vq_ff;

% ── PLL (co-rotating frame: dphi/dt = dtheta/dt - w0) ────────────────────
dphi     = p.Kp_pll*Vq_l + p.Ki_pll*eps_pll;
deps_pll = Vq_l;

% ── Power & current loop integrators ─────────────────────────────────────
deps_p  = Pref - P_meas;
deps_id = id_ref - id_loc;
deps_iq = iq_ref - iq_loc;

% ── LCL inductor dynamics (simplified: only Lf, no Cf state) ─────────────
did = (vd_vsc - Vd_l - Rf*id_loc)./Lf + w0*iq_loc;
diq = (vq_vsc - Vq_l - Rf*iq_loc)./Lf - w0*id_loc;

% ── GFF filter dynamics ───────────────────────────────────────────────────
dVd_ff = (Vd_l - Vd_ff) / tau;
dVq_ff = (Vq_l - Vq_ff) / tau;

% ── Assemble dxdt ─────────────────────────────────────────────────────────
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
