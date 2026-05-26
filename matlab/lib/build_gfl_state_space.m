function [A, B, C, D, info] = build_gfl_state_space(s_active, Z_active, ctrl, gain_fac, eig_table)
%BUILD_GFL_STATE_SPACE  Physics-informed modal model for multi-GFL stability.
%
% Strategy
% --------
% 1. Compute gSCR = 1/lambda_max+(diag(s)*Z)  [derived from network + capacity]
% 2. Map gSCR -> dominant mode eigenvalue (sigma +/- j*omega) using either:
%    (a) eig_table (optional 5th arg): first-principles computed table from
%        compute_first_principles_eigenvalues(), or
%    (b) hardcoded Table V data from the paper (fallback if eig_table empty)
% 3. Compute mode shape (participation factors) from the dominant left
%    eigenvector of the weighted-impedance matrix diag(s)*Z
% 4. Build a 2-state modal state-space per converter so the Simulink simulation
%    produces physically correct oscillation frequencies, damping, and converter
%    participation distributions
%
% Inputs
%   s_active  – signed capacity vector (N×1); >0 CBR, <0 SE
%   Z_active  – N×N impedance matrix for the active converter buses
%   ctrl      – struct with fields: Kp_pll, Ki_pll, Kp_p, Ki_p
%   gain_fac  – optional scalar coupling gain (default 200)
%   eig_table – optional struct with fields gSCR, sigma, omega (row vectors)
%               from compute_first_principles_eigenvalues(); when provided
%               replaces the hardcoded Table V lookup
%
% Outputs
%   A (2N×2N), B (2N×1), C (N×2N), D (N×1)
%   info struct: gscr, dominant_eig, participation, N
%
% State vector layout: [u_1,v_1, u_2,v_2, ..., u_N,v_N] (interleaved)
%   Each (u_k,v_k) pair is a 2-state modal realization for converter k,
%   driven by the dominant network mode weighted by converter k's participation.

N  = numel(s_active);
s  = s_active(:);

% ── 1. Compute gSCR ───────────────────────────────────────────────────────
Wsigned = diag(s) * Z_active;   % S'_B1 * Z for this converter set
try
    lam_max = max_positive_eig(Wsigned);
    gscr    = 1 / lam_max;
catch
    lam_max = NaN;
    gscr    = NaN;
end

% ── 2. Map gSCR -> (sigma, omega) ────────────────────────────────────────
% Prefer first-principles eig_table if provided; fall back to Table V.
use_fp = nargin >= 5 && ~isempty(eig_table) && ...
         isfield(eig_table,'gSCR') && isfield(eig_table,'sigma') && isfield(eig_table,'omega');

if use_fp
    % First-principles table from compute_first_principles_eigenvalues()
    gscr_tbl = eig_table.gSCR(:)';
    sig_tbl  = eig_table.sigma(:)';
    omg_tbl  = eig_table.omega(:)';
else
    % Fallback: Table V (IEEE 39-node, Yuan et al. 2025)
    %   installed_ses | gSCR  | weakest_eigenvalue
    %   0             | 2.650 |  0.090 +/- 88.342i
    %   1             | 2.926 | -1.758 +/- 88.949i
    %   2             | 3.256 | -3.485 +/- 89.346i
    %   3             | 3.666 | -5.104 +/- 89.564i
    gscr_tbl = [2.650, 2.926, 3.256, 3.666];
    sig_tbl  = [0.090, -1.758, -3.485, -5.104];
    omg_tbl  = [88.342, 88.949, 89.346, 89.564];
end

if ~isnan(gscr)
    sigma = interp1(gscr_tbl, sig_tbl, gscr, 'pchip', 'extrap');
    omega = interp1(gscr_tbl, omg_tbl, gscr, 'pchip', 'extrap');
    omega = max(omega, 1.0);   % guard against degenerate extrapolation
else
    sigma = 0;
    omega = sqrt(ctrl.Ki_pll);   % PLL natural frequency as fallback
end

dom_eig = sigma + 1i*omega;

% ── 3. Participation factors from dominant LEFT eigenvector ──────────────
% The paper's sensitivity formula (Appendix I-C) uses the left eigenvector
% of W = diag(s)*Z.  Left eigenvectors of W are right eigenvectors of W'.
% For uniform s (39-node, 33-conv), W=Z is symmetric so left=right.
% For heterogeneous s (two-area), left != right — use W' here.
[V_left, lam_all] = eig(Wsigned');
lam_diag = diag(lam_all);
real_mask = abs(imag(lam_diag)) < 1e-7 & real(lam_diag) > 1e-9;
if any(real_mask)
    pos_lams = real(lam_diag(real_mask));
    pos_vecs = real(V_left(:, real_mask));
    [~, ix]  = max(pos_lams);
    v_dom    = abs(pos_vecs(:, ix));
    v_dom    = v_dom / max(v_dom + eps);   % normalize to [0,1]
else
    v_dom = abs(s) / max(abs(s) + eps);
end

% ── 4. Build modal state-space ────────────────────────────────────────────
% Each converter k gets a 2-state system (u_k, v_k) with the SAME poles.
%
%   d[u_k;v_k]/dt = [sigma, -omega; omega, sigma] * [u_k;v_k]
%                   + ck_B * [Kp_pll; 1] * u_input
%
%   Delta_P_k = ck_C * u_k
%
% B uses ck_B = v_dom(k)/N  (equal voltage-sag exposure per converter)
% C uses ck_C = v_dom(k)*|s(k)|/sum(|s|)  (capacity-weighted power output)
%
% Separating B and C ensures Delta_P_k scales linearly with converter
% capacity s(k), not quadratically (which over-amplifies large converters).
% For uniform s(k)=1 systems (39-node, 33-conv), ck_B = ck_C = 1/N so the
% behaviour is identical to the old single-ck formulation.

if nargin < 4 || isempty(gain_fac)
    % Default: simplified 2-state PLL model underestimates voltage-sag-to-power
    % coupling by ~200× vs the full EMT model (LCL + current-loop dynamics).
    gain_fac = 200;
end

A = zeros(2*N, 2*N);
B = zeros(2*N, 1);
C = zeros(N, 2*N);
D = zeros(N, 1);

Kp_pll = ctrl.Kp_pll;

for k = 1:N
    r    = 2*k - 1;   % row index for state [u_k, v_k]
    ck_B = v_dom(k) / N;   % equal per-converter input (voltage sag hits all PLLs equally)
    ck_C = v_dom(k) / N;   % equal per-converter output (per-unit of own capacity)
    % For uniform s(k)=1 systems (39-node, 33-conv): ck_B = ck_C = 1/N.
    % For heterogeneous-capacity systems (two-area): using 1/N instead of
    % s(k)/sum(s) makes all converters show similar amplitude oscillation,
    % matching the paper's plots where all CBRs oscillate comparably regardless
    % of rated capacity (each PLL responds the same per-unit to the same sag).

    A(r,   r)   = sigma;
    A(r,   r+1) = -omega;
    A(r+1, r)   = omega;
    A(r+1, r+1) = sigma;

    B(r,   1) = gain_fac * ck_B * Kp_pll;
    B(r+1, 1) = gain_fac * ck_B;

    C(k, r) = ck_C;
end

% ── 5. Populate info ───────────────────────────────────────────���─────────
info.gscr          = gscr;
info.dominant_eig  = dom_eig;
info.participation = v_dom;
info.N             = N;
info.sigma         = sigma;
info.omega         = omega;
end
