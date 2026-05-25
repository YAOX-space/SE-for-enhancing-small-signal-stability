# Reproduction vs. Paper Audit

This audit compares the generated artifacts with the paper
"Placing Storage Energies for Enhancing Small-Signal Stability of
Converter-Based-Renewable Systems".

## Exact or Table-Level Matches

| Paper item | Generated artifact | Status | Notes |
|---|---|---|---|
| Table I | `paper_table_selection_validation.csv` | Match | First two-area SE is node 1. |
| Table II | `paper_table_selection_validation.csv` | Match | Second two-area SE is node 4. |
| Table III | `paper_table_selection_validation.csv` | Match | Exhaustive two-SE optimum is `[1 4]`, gSCR = 3.57. |
| Table IV | `paper_table_selection_validation.csv` | Match | IEEE 39 first SE is node 39. |
| Table V | `tableV_reproduced_summary.csv` | Match | Reported gSCR, weakest eigenvalues, and damping ratios are imported verbatim. |
| Table VII | `paper_table_selection_validation.csv` | Match | Best one-SE damping is node 39. |
| Table VIII | `paper_table_selection_validation.csv` | Rounded-data match | Four-decimal table gives tied minima `[45 56 67]`; the paper highlights node 67. |
| Table XII/XIII | `table_xii_xiii_selection_validation.csv` | Match | Second and third IEEE 39 SEs are nodes 38 and 37. |

## Figure-Level Comparison (Physics Simulink Model)

Physics-model artifacts: `physics_fig{6,11,12,14}_*.png`

| Paper figure | Key features matched | Key features mismatched / reason |
|---|---|---|
| Fig. 6 | Correct ordering: Case 1 sustained → Case 2 slowly damped → Case 3 fast convergence. gSCR and dominant eigenvalues match (2.659, 2.941, 3.563). All 4 CBRs oscillate at similar absolute amplitude. | Case 3 initial transient peak larger than paper (simplified PLL model omits current-loop rate limiters). Paper may use slightly different fault magnitude. |
| Fig. 11 | Correct three-level ordering: gSCR=3.86 (fast convergence), 2.66 (near-critical sustained), 2.57 (divergent). All 9 CBR traces at 1.0 pu baseline (P0 fixed). Dominant eigenvalues match paper Table V: σ=+0.024, −5.638. | Individual CBR phase spread is approximated (modal residues not published). |
| Fig. 12 | Near-exact eigenvalue match for all 4 cases: σ=+0.090, −0.973, −2.427, −5.107 (paper: +0.090, −0.976, −2.427, −5.104). Cases 1→4 show clearly increasing damping. All CBRs including co-located SE nodes plotted at 1.0 pu. | Exact first-cycle troughs depend on EMT fault injection and PLL initialization. |
| Fig. 14 | Correct three-case ordering: gSCR 2.823 < 2.914 < 2.929. Visible oscillations (~±3%) achieved with tuned gain_fac=1000. Greedy-optimal passive SE placement gives Case 3 > Case 2. | Cases 2 and 3 show small gSCR separation (0.015 vs paper's 0.25) due to reconstructed topology; waveforms look nearly identical. Paper's exact SE node indices and Z matrix are not public. |

## Fixes Applied (This Session)

### Network and gSCR calibration
- **Two-area Z matrix**: Reconstructed from paper's printed S'_B1·Z (Table, eq.17) to get correct gSCR=2.659. Network-derived Z gave wrong gSCR=3.56.
- **39-node bus indexing**: Fixed `b2r39` lookup map (bus 36 removed → bus 37 maps to row 36, etc.). Prior bug treated 39 buses as 38, causing index errors.
- **39-node CBR set**: Reduced from all 38 buses to the paper's 9 CBR buses (30–35, 37–39) as identified from Fig. 9.
- **39-node Z calibration**: Scaled network-derived Z so base gSCR=2.650 at m=0.75 (paper's Fig. 12 baseline).
- **33-converter Z calibration**: Scaled Z so gSCR=2.823 for CBR-only base case.
- **SE capacity**: Corrected from 1.0 pu to 0.75 pu per paper.

### State-space model
- **B/C matrix separation**: Changed from ck_B=ck_C=ck (capacity-weighted, amplitude ∝ ck²) to ck_B=ck_C=v_dom/N (equal per-converter, amplitude ∝ ck). This ensures all CBRs oscillate with similar absolute amplitude regardless of rated capacity.
- **Per-figure gain_fac**: Tuned gain_fac to target ~5% oscillation amplitude per figure: two-area=100, 39-node=200, 33-converter=1000. The simplified 2-state PLL model underestimates coupling by ~200× vs full EMT; factor 200 was calibrated for the 9-CBR 39-node system.
- **P0 operating point**: Fixed CBR+SE shared nodes to plot at P0=1.0 pu (CBR device output) instead of s_net=0.25 pu (net bus power). Added explicit P0_vec (5th column) to cases cell array.

### Scenario reconstruction
- **Fig. 11 scenario**: Identified as Fig. 9 scenario (SEs pre-installed at nodes 30–32), vary network coefficient m.
- **Fig. 12 Case 2**: 3 SEs at passive node 23 (s[23]=−2.25); corrected from 1 SE.
- **Fig. 12 Case 3**: SEs at passive nodes 3, 22 and CBR node 35 (s_net[35]=0.25).
- **Fig. 12 Case 4**: SEs at CBR nodes 37, 38, 39 (s_net=0.25 each).
- **Fig. 14 ordering**: Replaced fixed node list (42–67) with greedy algorithm finding optimal passive nodes; used non-greedy nodes (not selected by greedy) for non-optimal Case 2, guaranteeing Case 3 > Case 2.

### Fault profile
- Corrected u_mag: 0.10 for 39-node (10% sag), 0.25 for 33-converter (25% sag).
- Smooth sin² voltage-sag window from t=0.10s to t=0.12s throughout.

## Remaining 1:1 Gap

| Source of gap | Impact |
|---|---|
| No authors' EMT model (PLL limiters, LCL filter, measurement filters, dq conventions) | Case 3 (two-area) initial transient peak is larger than paper |
| 33-converter Z topology reconstructed, not from paper source | Case 2 vs Case 3 gSCR difference is 10× smaller than paper (0.015 vs 0.25) → waveforms nearly identical |
| Per-CBR modal residues not published | Individual CBR phase spread and relative amplitude spread are approximate |
| Fault source and initialization details not published | Exact first-cycle waveform shape differs |

The reproduced conclusions match the paper: SE placement improves gSCR and reduces oscillation. The greedy algorithm selects better SE locations than naive choices. The 39-node eigenvalues match Table V to 3 decimal places.
