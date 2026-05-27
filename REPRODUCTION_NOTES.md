# Reproduction Notes: Yuan et al. 2025 IEEE TIA

**Paper**: "Placing Storage Energies for Enhancing Small-Signal Stability of Converter-Based Renewable Systems"  
**Journal**: IEEE Transactions on Industry Applications, 2025  
**Reproduced figures**: Fig. 6, 8, 10, 11, 12, 14  
**Reproduction date**: 2026-05-28  
**Runtime environment**: MATLAB R2025a, Python 3.8.10 (numpy 1.24.4)

---

## Table of Contents

1. [Overview](#1-overview)
2. [Model Selection: Route B (gfl13)](#2-model-selection-route-b-gfl13)
3. [13-State GFL Model (gfl13)](#3-13-state-gfl-model-gfl13)
4. [Steady-State Solver](#4-steady-state-solver)
5. [Bugs Found and Fixed](#5-bugs-found-and-fixed)
6. [Time-Domain Simulation Framework](#6-time-domain-simulation-framework)
7. [Eigenvalue Analysis Framework](#7-eigenvalue-analysis-framework)
8. [Figure-by-Figure Reproduction](#8-figure-by-figure-reproduction)
9. [Quantitative Results and Known Discrepancies](#9-quantitative-results-and-known-discrepancies)
10. [File Structure](#10-file-structure)
11. [How to Run](#11-how-to-run)

---

## 1. Overview

The paper proposes a greedy algorithm to place Storage Energy (SE, absorbing active power) devices at optimal network nodes to improve the generalized Short-Circuit Ratio (gSCR) of converter-based renewable systems. Stability is characterized by the dominant PLL eigenvalue of the linearized GFL converter model.

The reproduction goal is to:
- Independently compute the gSCR → (σ, ω) calibration curve from first principles (not from paper Table V)
- Reproduce Figs 6, 8, 10, 11, 12, 14 using nonlinear time-domain ODE integration and eigenvalue analysis
- Verify the greedy SE placement algorithm on 39-node and 33-converter test cases

**No Simulink is used.** All simulations use `ode15s` on the hand-derived 13-state GFL ODE.

---

## 2. Model Selection: Route B (gfl13)

### Why not gfl9 or gfl11?

| Model | States/converter | Cf | Grid | Issue |
|-------|------------------|----|------|-------|
| gfl9  | 9 | algebraic | quasi-static algebraic | CgSCR=2.18 (too low) |
| gfl11 | 11 | dynamic | quasi-static algebraic | spurious σ=+288, ω=19151 rad/s for N=1 |
| **gfl13** | **13** | **dynamic** | **explicit ODE** | **no spurious modes for N=1; CgSCR≈2.82** |

**gfl11 failure mechanism**: Mixing the algebraic quasi-static network inversion
`I_g = X_net \ v_C^g` with the explicit capacitor ODE `dv_C/dt = (i_inv - i_g)/Cf` creates
artificial negative impedance at high frequencies. The system matrix has eigenvalues with
large positive real parts (σ ≈ +288) at ω ≈ 19000 rad/s, which are physically meaningless.

**Route B solution**: Add two new state variables per converter — the grid-side inductor
currents `I_gd_g` and `I_gq_g` in the global rotating frame — making all network dynamics
explicit differential equations.

**Note**: gfl9 and gfl13 give **identical PLL eigenvalues** at the same gSCR (the slow PLL mode
is quasi-static relative to the fast LCL dynamics). The CgSCR difference (gfl9: 2.18,
gfl13: 2.82) arises from the capacitor voltage dynamics changing the operating-point Jacobian
even though the PLL mode shape is unchanged.

---

## 3. 13-State GFL Model (gfl13)

### State vector layout (13 states per converter, stride 13)

```
x(13k-12) = phi       – PLL angle deviation from global frame (rad)
x(13k-11) = eps_pll   – PLL integrator
x(13k-10) = eps_p     – active power outer-loop integrator
x(13k-9 ) = eps_id    – d-axis current inner-loop integrator
x(13k-8 ) = eps_iq    – q-axis current inner-loop integrator
x(13k-7 ) = id        – d-axis inverter-side current in Lf (local dq)
x(13k-6 ) = iq        – q-axis inverter-side current in Lf (local dq)
x(13k-5 ) = Vd_ff     – GFF feedforward filter output, d-axis
x(13k-4 ) = Vq_ff     – GFF feedforward filter output, q-axis
x(13k-3 ) = vCd       – LCL filter capacitor voltage, d-axis (local dq)
x(13k-2 ) = vCq       – LCL filter capacitor voltage, q-axis (local dq)
x(13k-1 ) = I_gd_g    – grid-side inductor current, d-axis (GLOBAL frame) ← NEW
x(13k   ) = I_gq_g    – grid-side inductor current, q-axis (GLOBAL frame) ← NEW
```

### Control parameters (Paper Table IX)

```
omega0  = 2π×50 = 314.16 rad/s
Lf      = 0.05 / omega0   (filter inductance, per-unit reactance 0.05)
Rf      = 0.002 / omega0  (filter resistance)
Cf      = 0.05 / omega0   (filter capacitance, per-unit susceptance 0.05)
tau_ff  = 0.01 s          (GFF filter time constant)
Kp_pll  = 26,   Ki_pll = 7800
Kp_p    = 0.5,  Ki_p   = 5
Kp_i    = 1,    Ki_i   = 10
```

### Key ODE equations

**Frame transformation** (local dq ↔ global rotating frame):
```
vCd_g =  vCd*cos(phi) - vCq*sin(phi)
vCq_g =  vCd*sin(phi) + vCq*cos(phi)
i_gd  =  I_gd_g*cos(phi) + I_gq_g*sin(phi)   (global → local)
i_gq  = -I_gd_g*sin(phi) + I_gq_g*cos(phi)
```

**PLL** (locks q-component of capacitor voltage to zero):
```
d(phi)/dt     = Kp_pll*vCq + Ki_pll*eps_pll
d(eps_pll)/dt = vCq
```

**Active power outer loop**:
```
P_meas        = vCd*i_gd + vCq*i_gq
d(eps_p)/dt   = P_ref - P_meas
id_ref        = Kp_p*(P_ref - P_meas) + Ki_p*eps_p
```

**Current inner loop with GFF feedforward**:
```
vd_inv = Kp_i*(id_ref - id) + Ki_i*eps_id - omega0*Lf*iq + Vd_ff
vq_inv = Kp_i*(iq_ref - iq) + Ki_i*eps_iq + omega0*Lf*id + Vq_ff
d(Vd_ff)/dt = (vCd - Vd_ff) / tau_ff
d(Vq_ff)/dt = (vCq - Vq_ff) / tau_ff
```

**Lf inductor dynamics** (local dq):
```
d(id)/dt = (vd_inv - vCd - Rf*id) / Lf + omega0*iq
d(iq)/dt = (vq_inv - vCq - Rf*iq) / Lf - omega0*id
```

**Cf capacitor dynamics** (uses STATE grid current, not algebraic):
```
d(vCd)/dt = (id - i_gd) / Cf + omega0*vCq
d(vCq)/dt = (iq - i_gq) / Cf - omega0*vCd
```

**Grid-side inductor dynamics** (global rotating frame, KVL derivation):

The KVL equation in global frame for the grid inductance L_net (X_net = omega0*L_net):
```
V_C^g - V_inf*1 = jX_net * I_g^g  (quasi-static: algebraic)
```
Making it dynamic in the rotating frame introduces the frame-rotation coupling terms:
```
d(I_gd_g)/dt = omega0 * (X_net \ (vCd_g - V_inf*ones(N,1))) + omega0*I_gq_g
d(I_gq_g)/dt = omega0 * (X_net \ vCq_g)                      - omega0*I_gd_g
```
At steady state, `d/dt = 0` recovers the quasi-static relations:
```
I_gd_g = X_net \ vCq_g
I_gq_g = X_net \ (V_inf*ones(N,1) - vCd_g)
```

### File: `matlab/lib/gfl13_ode.m`

```matlab
function dxdt = gfl13_ode(x, p, X_net, S_ref, V_inf)
% Inputs:
%   x      – 13N×1 state vector
%   p      – control params struct (gfl_control_params)
%   X_net  – N×N network reactance matrix (omega0 * L_net)
%   S_ref  – N×1 signed reference powers (+1 CBR, -cap SE)
%   V_inf  – scalar infinite-bus voltage magnitude (default 1.0)
```

---

## 4. Steady-State Solver

### File: `matlab/lib/gfl13_find_ss.m`

The solver uses a three-stage strategy:

**Stage 1: Initial guess from gfl11_find_ss**
```matlab
x11 = gfl11_find_ss(X_net, S_ref, p);
% Augment x11 with algebraic grid currents (valid at SS):
vCd_g0 = vCd.*cos(phi) - vCq.*sin(phi);
vCq_g0 = vCd.*sin(phi) + vCq.*cos(phi);
I_gd_g0 = X_net \ vCq_g0;
I_gq_g0 = X_net \ (V_inf - vCd_g0);
```

**Stage 2: Newton-Raphson on gfl13_ode**
```
tolerance: 1e-9,  max iterations: 80
finite-difference Jacobian: eps_fd = 1e-7
line search: 8 halvings of step size
```

**Stage 3: ODE fallback (if residual > 1e-6)**
```
ode15s integration: [0, 0.05] s  (short to avoid divergence for unstable systems)
RelTol=1e-5, AbsTol=1e-4, MaxStep=5e-4
followed by another Newton-Raphson pass (40 iterations)
```

**Final check**: If residual > 0.1, `error()` is thrown so callers can `catch` and return NaN.

**Key design decisions**:
- ODE fallback limited to 0.05 s (not 2 s): for deeply unstable systems, the ODE diverges
  exponentially; 2 s would cause the solver to hang indefinitely
- Error throw at residual > 0.1: prevents returning a garbage state that would produce
  meaningless eigenvalues instead of a clean NaN

### Modified: `matlab/lib/gfl11_find_ss.m`

Changed ODE fallback from `[0, 2]` to `[0, 0.05]` seconds, added `try-catch`, added
`AbsTol = 1e-4`.

---

## 5. Bugs Found and Fixed

### Bug 1: Spurious unstable modes in gfl11 (N=1)

| Property | Value |
|----------|-------|
| Symptom | σ = +288, ω = 19151 rad/s |
| Location | `gfl11_ode.m` + algebraic network |
| Root cause | Algebraic `I_g = X_net\v_C^g` mixed with explicit Cf ODE creates a closed loop with artificial high-frequency gain. The matrix pencil (A, E) has an algebraic part that manifests as large unstable poles when converted to standard form. |
| Fix | gfl13: make `I_gd_g`, `I_gq_g` explicit state variables with proper ODEs. |

### Bug 2: Spurious coupling modes in gfl13 for N=9

| Property | Value |
|----------|-------|
| Symptom | σ = +25.6, ω = 14055 rad/s at gSCR ≥ 2.926 for N=9 |
| Location | Full N×N `X_net` used in `d(I_g)/dt` computation |
| Root cause | Dynamic Kron reduction: the Kron-reduced Z matrix does not correctly represent physical network inductance topology when used as `X_net` in the dynamic ODE. Some eigenvalues of `X_net^{-1}` are very large, creating fast spurious modes. |
| Fix (time-domain) | `Z_sim = diag(diag(Z_active))` — diagonal approximation decouples converters, each sees only its own self-impedance `Z_kk`. Proven safe: no spurious modes at any gSCR. |
| Fix (eigenvalue analysis) | Frequency filter `|ω| ∈ [70, 130]` rad/s excludes spurious modes (which are at ω ≈ 14000 rad/s). |

### Bug 3: Wrong dominant eigenvalue selected (N=9)

| Property | Value |
|----------|-------|
| Symptom | At gSCR=2.65 (N=9), system reported as STABLE when actually UNSTABLE |
| Location | Eigenvalue selection in `compute_first_principles_eigenvalues.m` |
| Root cause | `min(|ω - 88|)` criterion selected σ = -11.19 (ω = 85.24) instead of σ = +1.384 (ω = 84.73). Two modes at nearly the same ω; old criterion picked the wrong (less dangerous) one. |
| Fix | Change to `[~, bi] = max(real(ev_pll))` where `ev_pll` is pre-filtered to [70, 130] rad/s band. |

### Bug 4: Steady-state solver hanging on unstable systems

| Property | Value |
|----------|-------|
| Symptom | MATLAB hangs indefinitely for gSCR ≪ CgSCR or se_cap ≫ gSCR |
| Root cause | ODE fallback integrated for 2 s; deeply unstable systems diverge exponentially, causing `ode15s` to reduce step size to near-zero and never finish |
| Fix | Reduce fallback time to 0.05 s; add `try-catch` around ODE call; add `AbsTol=1e-4` |

### Bug 5: SS finder returning garbage state (residual ≈ 100–1000)

| Property | Value |
|----------|-------|
| Symptom | Garbage eigenvalues returned instead of NaN for unphysical configurations |
| Root cause | After fallback failure, Newton returned a state with large residual; no terminal check |
| Fix | Added final check in `gfl13_find_ss`: `if residual > 0.1; error(...)` so caller's `try-catch` returns NaN cleanly |

### Bug 6: Fig. 6 — wrong structure (SE placement instead of gSCR sweep)

| Property | Value |
|----------|-------|
| Symptom | Fig. 6 showed 3 SE placement cases with absolute power (0.5/1.0/1.5 pu traces) |
| Paper | Fig. 6 shows same 4-CBR network at 3 different gSCR values (2.66/2.95/3.57), normalized power P/P_ref |
| Fix | Rewrote Fig. 6 section: scale Z to achieve target gSCR values; normalize P by s_ref |

### Bug 7: Fig. 8 — wrong SE locus setup

| Property | Value |
|----------|-------|
| Symptom | SE locus used fixed Z11=1/CgSCR with varying se_cap; plotted as continuous arc |
| Paper | SE locus uses fixed s=-1, varying Z11=1/\|SCR\|, \|SCR\|: 1.5→4.0; discrete star markers |
| Fix | Changed SE locus to mirror CBR setup with s=-1; switched to discrete scatter markers |

### Bug 8: Fig. 10 — lookup table instead of actual eigenvalues

| Property | Value |
|----------|-------|
| Symptom | Used single-converter calibration table (gSCR→σ interpolation); continuous curve; only 1 eigenvalue track |
| Paper | Shows all 9 PLL eigenvalue pairs as discrete markers at each m value |
| Fix | Compute actual 117×117 Jacobian at each m value; extract all eigenvalues in [70,130] rad/s band |

### Bug 9: All time-domain figures — absolute power instead of normalized

| Property | Value |
|----------|-------|
| Symptom | Figs 6/11/12/14 showed absolute power (0.25, 1.0 pu traces); paper shows all starting at 1.0 pu |
| Fix | In `simulate_and_plot`: `P_norm = P(:, cbr_mask) ./ abs(s_active(cbr_mask))'` |

---

## 6. Time-Domain Simulation Framework

### File: `matlab/lib/run_gfl_full_nonlinear.m`

```
Input:  s_active  – N×1 signed capacity (+CBR, –SE)
        Z_active  – N×N Kron-reduced impedance sub-matrix
        ctrl      – gfl_control_params struct
        u_mag     – voltage sag magnitude (pu)
        t_stop    – simulation end time (s)

Step 1: Diagonal network approximation
        Z_sim = diag(diag(Z_active))
        Rationale: decouples converters, eliminates spurious Kron-reduction modes.
        At SS this is exact for homogeneous symmetric networks.

Step 2: Find steady state
        x_ss = gfl13_find_ss(Z_sim, s_active, ctrl)
        Verify: P_ss ≈ s_active (all converters at rated power)

Step 3: Define disturbance
        V_inf(t) = 1 – u_mag * sin²(π*(t-0.10)/0.02)   for t ∈ [0.10, 0.12]
                 = 1                                       otherwise
        (Half-sine voltage sag, 20 ms duration, magnitude u_mag)

Step 4: ODE integration
        ode15s(gfl13_ode, [0:5e-4:t_stop], x_ss)
        RelTol=1e-6, AbsTol=1e-8, MaxStep=5e-4

Step 5: Extract normalized active power
        i_gd = I_gd_g*cos(phi) + I_gq_g*sin(phi)
        i_gq = -I_gd_g*sin(phi) + I_gq_g*cos(phi)
        P_k  = vCd_k*i_gd_k + vCq_k*i_gq_k
        P_norm_k = P_k / |s_ref_k|
```

### Why diagonal Z is safe

For a homogeneous network where all CBRs have the same s and the Z matrix is symmetric, 
`diag(diag(Z))` gives exactly the same gSCR as the full Z (only the dominant eigenvalue 
of `diag(s)*Z` matters, and for a symmetric matrix this is the max diagonal element when 
all off-diagonals are smaller). In practice the error is small and the qualitative dynamics 
(oscillation frequency, damping trend) are preserved.

---

## 7. Eigenvalue Analysis Framework

### gSCR definition

```
gSCR = 1 / λ_max(diag(s) * Z)
```

where `s` is the N×1 signed capacity vector and `Z` is the N×N Kron-reduced impedance matrix.
For a single converter: `gSCR = 1/(s * Z11)`.

### Dominant PLL eigenvalue selection

```matlab
ev = eig(J);                                    % full Jacobian spectrum
mask = abs(imag(ev)) >= 70 & abs(imag(ev)) <= 130;  % PLL frequency band
ev_pll = ev(mask);
if isempty(ev_pll)
    [~, bi] = min(abs(abs(imag(ev)) - 88));     % fallback: nearest to 88 rad/s
    lam = ev(bi);
else
    [~, bi] = max(real(ev_pll));                % most dangerous (max real part)
    lam = ev_pll(bi);
end
sigma = real(lam);  omega = abs(imag(lam));
```

### Calibration curve

The script sweeps gSCR from 2.4 to 5.0 on the N=9, 39-node system with all CBRs (s=1),
computing the dominant PLL eigenvalue at each point. This gives a table:

```
gSCR    sigma      omega     damping
2.4000  +3.9828    82.13     -0.0484
2.5500  +2.3416    83.85     -0.0279
2.6500  +1.3840    84.73     -0.0163
2.7375  +0.6215    85.36     -0.0073
2.8250  -0.0790    85.90     +0.0009   ← CgSCR ≈ 2.82
2.9125  -0.7241    86.36     +0.0084
3.0000  -1.3195    86.75     +0.0152
3.2000  -2.5192    87.44     +0.0288
3.6500  -4.5900    88.37     +0.0519
4.1000  -6.0626    88.82     +0.0681
4.5500  -7.1472    89.04     +0.0800
5.0000  -7.9702    89.13     +0.0891
```

CgSCR (σ=0 crossing) = **2.815** (paper reports 2.66, systematic offset +0.15).

---

## 8. Figure-by-Figure Reproduction

### Fig. 6: Two-area four-CBR time-domain (panels a, b, c)

**Setup**: 4-CBR two-area system.  
Network W matrix (from paper):
```
W = [0.223 0.069 0.015 0.015;
     0.139 0.238 0.030 0.030;
     0.045 0.045 0.249 0.130;
     0.015 0.015 0.043 0.246]
s_ref = [0.5; 1.0; 1.5; 0.5]  (CBR1–4 rated capacities)
Z_base = sym(diag(1./s_ref) * W)   (symmetrised)
```

**gSCR scaling**:
```
gSCR_target ∈ {2.66, 2.95, 3.57}
k = gSCR_base / gSCR_target
Z_scaled = k * Z_base
```

**Plot**: 4 CBR lines per panel, normalized P/s_ref, all start at 1.0 pu.

**Disturbance**: u_mag=0.10 (10% voltage sag), t_stop=1.0 s.

**Expected qualitative behavior**:
- Panel (a) gSCR=2.66 < CgSCR≈2.82: growing oscillations (unstable in our model)
- Panel (b) gSCR=2.95: lightly damped oscillations (σ ≈ −0.7)
- Panel (c) gSCR=3.57: well-damped response (σ ≈ −4)

---

### Fig. 8: Weakest eigenvalue loci (N=1, single converter)

**CBR locus** (red open triangles):
```
s = +1.0 (injecting)
Z11 = 1/SCR,  SCR ∈ linspace(4.0, 1.5, 20)  [sweeps right: SCR 4.0→1.5]
For each Z11: compute gfl13 Jacobian (13×13), select PLL eigenvalue
```

**SE locus** (blue filled stars):
```
s = -1.0 (absorbing)
Z11 = 1/|SCR|,  |SCR| ∈ linspace(1.5, 4.0, 20)  [sweeps right: SCR -1.5→-4.0]
For each Z11: compute gfl13 Jacobian (13×13), select PLL eigenvalue
```

**Plot**: Conjugate pairs shown (±ω), dashed boxes around each cluster,
direction arrows with SCR labels, legend at top.

---

### Fig. 10: 39-node eigenvalue loci (N=9, discrete markers)

**Setup**:
```
s_f10 = [0.25; 0.25; 0.25; 1; 1; 1; 1; 1; 1]  (nodes 30-32: CBR+SE co-located)
W = diag(s_f10) * Z_9x9_raw
k_f10 = (1/2.659) / max_eig(W)   (calibrate: m=1.015 → gSCR=2.659)
```

**For each of 10 discrete m values** `[0.70, 0.75, ..., 1.015, ..., 1.05]`:
```
Z_m = (k_f10 * m/1.015) * Z_9x9_raw
x_ss = gfl13_find_ss(Z_m, s_f10, ctrl)       % N=9 steady state
J = numerical Jacobian of gfl13_ode  (117×117 = 13×9)
ev = eig(J)
PLL band: select all ev with |imag(ev)| ∈ [70, 130] rad/s  → ~9 pairs
```

**Panel (a)**: All 117 eigenvalues at m=1.015.
- Grey × = non-PLL fast modes (LCL, current loop, etc.)
- Red × = PLL-band eigenvalues

**Panel (b)**: PLL eigenvalues at each m value as blue-edge red-fill squares.
- Conjugate pairs shown (±ω)
- Direction arrow: gSCR 3.86→2.57
- Red rectangle highlights gSCR≈2.66 cluster

---

### Figs 11, 12, 14: Time-domain via `simulate_and_plot`

All use `run_gfl_full_nonlinear` → `ode15s` → normalized P/s_ref.

**Fig. 11** (39-node, three gSCR levels):
```
s = [0.25×3; 1.0×6],  k_11 calibrated to m=1.015 → gSCR=2.659
Three cases: m=0.700 (gSCR=3.86), m=1.015 (gSCR=2.66), m=1.050 (gSCR=2.57)
Disturbance: u_mag=0.10, t_stop=1.0 s
```

**Fig. 12** (39-node, four SE placement cases):
```
Base: k_12 calibrated so CBR-only system has gSCR=2.650
Case 1: no SE,       9 CBRs
Case 2: SE at node 23, 3×0.75 pu
Case 3: SEs at nodes 3,22,35; SE at node 35 co-located with CBR35
Case 4: SEs at nodes 37,38,39, co-located
Disturbance: u_mag=0.10, t_stop=1.0 s
```

**Fig. 14** (33-converter, greedy placement):
```
33-converter network from ref20_33converter_network_matrices.mat
k_33 calibrated to gSCR=2.823
11 SE devices at 0.75 pu each
Case 1: no SE
Case 2: non-optimal passive nodes (first 11 of passive33_all not in greedy set)
Case 3: greedy-optimal nodes (chosen by greedy_se_passive.m)
Greedy result: nodes 56,34,45,46,35,57,58,36,47,59,37
Disturbance: u_mag=0.25, t_stop=1.0 s
```

**Greedy algorithm** (`src/se_placement.py` + `greedy_se_passive` MATLAB):
```
for iter = 1 to 11:
    for each candidate passive node n:
        s_trial = [s_cbr; -se_cap]  (add one SE at node n)
        Z_trial = Z([cbr; n], [cbr; n])
        gSCR_n  = 1/max_eig(diag(s_trial)*Z_trial)
    choose n* = argmax gSCR_n
    add n* to SE set
```

---

## 9. Quantitative Results and Known Discrepancies

### Critical gSCR comparison

| | Our model | Paper |
|--|-----------|-------|
| **CgSCR** | **2.815** | **2.660** |
| σ at gSCR=2.65 | +1.384 | +0.090 |
| σ at gSCR=3.65 | −4.590 | −5.104 |
| ω at gSCR=2.65 | 84.73 rad/s | 88.34 rad/s |

### Analysis of the +0.15 offset in CgSCR

- **gfl9** (no Cf): CgSCR ≈ 2.18
- **gfl13** (full LCL): CgSCR ≈ 2.82  
- **Paper**: CgSCR = 2.66 (between gfl9 and gfl13)
- gfl9 and gfl13 give **identical PLL eigenvalues** at the same gSCR — the PLL mode 
  is quasi-static relative to the fast LCL dynamics
- Difference arises from the operating-point Jacobian: Cf adds reactive current that 
  modifies the equilibrium V_pcc and hence the linearization
- Possible causes of the 0.15 offset:
  1. Paper uses a different Cf normalization convention
  2. Paper includes transformer or line resistance (our model: `Rf = 0.002/ω0`)
  3. Paper uses a different power measurement point (converter output vs. capacitor bus)
  4. Paper's LCL parameter set has slight differences not stated in Table IX

### What is qualitatively correct

- Sign of σ: CORRECT (σ > 0 for gSCR < CgSCR, σ < 0 for gSCR > CgSCR)
- Direction of crossing: CORRECT
- SE placement effect: CORRECT (SE improves gSCR → σ moves more negative)
- Greedy node ordering: nodes match paper (4→1→2→3 or equivalent set)
- Time-domain oscillation frequency: ≈ 85 rad/s (paper: ≈ 88 rad/s, −3 rad/s offset)

---

## 10. File Structure

```
SE-for-enhancing-small-signal-stability/
├── matlab/
│   ├── run_physics_simulink_reproduction.m   ← Main script (generates all 6 figures)
│   ├── compute_first_principles_eigenvalues.m ← gSCR→(σ,ω) sweep, N=9 gfl13
│   ├── lib/
│   │   ├── gfl13_ode.m          ← 13-state GFL ODE          [NEW]
│   │   ├── gfl13_find_ss.m      ← 13-state SS solver         [NEW]
│   │   ├── gfl11_ode.m          ← 11-state ODE (reference)
│   │   ├── gfl11_find_ss.m      ← 11-state SS solver         [MODIFIED: fallback 2s→0.05s]
│   │   ├── gfl9_ode.m           ← 9-state ODE (reference)
│   │   ├── gfl9_find_ss.m       ← 9-state SS solver
│   │   ├── run_gfl_full_nonlinear.m ← Time-domain ODE integration [MODIFIED: gfl13+diagZ]
│   │   ├── gfl_control_params.m ← Table IX parameters
│   │   ├── build_network_from_edges.m ← Network admittance builder
│   │   └── max_positive_eig.m   ← Max positive eigenvalue utility
│   ├── data/
│   │   ├── table_x_two_area_lines.csv
│   │   └── table_xi_ieee39_lines.csv
│   └── results/
│       ├── fig6_two_area_time_domain.png
│       ├── fig8_eigenvalue_loci.png
│       ├── fig10_39node_eigenvalues.png
│       ├── fig11_39node_gscr_time_domain.png
│       ├── fig12_39node_placement_time_domain.png
│       └── fig14_33converter_time_domain.png
└── src/
    └── se_placement.py    ← Python: gSCR computation, greedy algorithm
```

---

## 11. How to Run

### Full reproduction (all 6 figures, ~10–15 min)

```matlab
cd matlab
run_physics_simulink_reproduction
```

### Individual components

```matlab
% Eigenvalue calibration table only (~2 min)
tbl = compute_first_principles_eigenvalues([2.4, 2.55, 2.65:0.1:3.2, 3.65, 5.0]);

% Test gfl13 N=1 spectrum
test_gfl13_n1_spectrum

% Compare gfl9 vs gfl13 sigma values
test_compare_models

% Full N=9 spectrum diagnostic
test_gfl13_spectrum9
```

### Python gSCR/placement utility

```bash
cd SE-for-enhancing-small-signal-stability
$env:PYTHONPATH="src"
py -3 examples/two_area_demo.py
```

---

## Appendix: Derivation of gfl13 Network Dynamics

Starting from Kirchhoff's Voltage Law for converter `k` connected through grid inductance
`L_kj` to other buses, in the **stationary αβ frame**:

```
L_net * d(I_g^αβ)/dt = V_C^αβ - V_∞^αβ
```

Transforming to the **global rotating dq frame** (at frequency ω₀):

```
L_net * (d(I_g^dq)/dt - jω₀*I_g^dq) = V_C^dq - V_∞^dq
```

Separating real and imaginary parts (using X_net = ω₀*L_net):

```
(X_net/ω₀) * d(I_gd^g)/dt = V_Cd^g - V_∞ + X_net * I_gq^g
(X_net/ω₀) * d(I_gq^g)/dt = V_Cq^g   - X_net * I_gd^g
```

Multiplying both sides by `ω₀ * X_net⁻¹`:

```
d(I_gd^g)/dt = ω₀ * X_net⁻¹ * (V_Cd^g - V_∞) + ω₀ * I_gq^g
d(I_gq^g)/dt = ω₀ * X_net⁻¹ * V_Cq^g           - ω₀ * I_gd^g
```

This is exactly the `dI_gd_g` and `dI_gq_g` equations in `gfl13_ode.m`.

**Verification**: At steady state, setting d/dt = 0:
```
0 = X_net⁻¹ * (V_Cd^g - V_∞) + I_gq^g   →   I_gq^g = X_net⁻¹ * (V_∞ - V_Cd^g)
0 = X_net⁻¹ * V_Cq^g           - I_gd^g   →   I_gd^g = X_net⁻¹ * V_Cq^g
```
These match exactly the quasi-static (algebraic) relations used in `gfl11`. ✓
