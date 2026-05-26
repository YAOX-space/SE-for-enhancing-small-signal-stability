"""Reproduce the paper's two-area four-CBR storage placement example.

Tables I & II (greedy placement) and Table III (exhaustive 2-SE search) are
both reproduced from the network line data in
matlab/data/table_x_two_area_lines.csv.

Network construction: the susceptance Laplacian is built from all 9 nodes and
inverted via Moore-Penrose pseudo-inverse (build_network_z_pinv).  This
preserves the inter-area coupling that runs through node 9 (the infinite-bus
reference).  Removing node 9 as in build_network_from_edges.m severs the only
path between CBRs 1/2 and CBRs 3/4 (both areas connect only through node 9).

Calibration: Z is scaled so that λmax(diag(S_base) Z) = λmax(INITIAL_WEIGHTED),
where S_base = [0.5, 1.0, 1.5, 0.5, 0, …, 0].  This anchors the per-unit
scale to the paper without back-computing Z from the printed W (which would
be circular: W → Z → same W).

INITIAL_WEIGHTED is kept only as a calibration anchor and for direct
comparison with the paper's reported values; it is NOT used as algorithm input.
"""

from __future__ import annotations

import csv
import os
from itertools import combinations_with_replacement
from pathlib import Path

import numpy as np

from se_placement import greedy_place_storage, gscr_from_weighted_impedance
from se_placement import weighted_impedance, lambda_max_positive


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / ".matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(RESULTS / ".cache"))
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt


INITIAL_WEIGHTED = np.array(
    [
        [0.223, 0.069, 0.015, 0.015, 0.069, 0.030, 0.015, 0.015, 0.015],
        [0.139, 0.238, 0.030, 0.030, 0.139, 0.060, 0.030, 0.030, 0.030],
        [0.045, 0.045, 0.249, 0.130, 0.045, 0.045, 0.045, 0.071, 0.130],
        [0.015, 0.015, 0.043, 0.246, 0.015, 0.015, 0.015, 0.024, 0.043],
        [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
        [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
        [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
        [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
        [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
    ],
    dtype=float,
)


def build_laplacian(lines_path: Path, n_nodes: int) -> np.ndarray:
    """Build the full n_nodes × n_nodes susceptance Laplacian from a CSV."""
    Q = np.zeros((n_nodes, n_nodes))
    with open(lines_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            i = int(row["from"]) - 1   # 0-indexed
            j = int(row["to"]) - 1
            b = 1.0 / float(row["x"])
            Q[i, i] += b
            Q[j, j] += b
            Q[i, j] -= b
            Q[j, i] -= b
    return Q


def build_network_z(lines_path: Path, n_nodes: int, ground_node: int) -> np.ndarray:
    """Python mirror of build_network_from_edges.m.

    Builds the susceptance Laplacian, removes the ground node row/column,
    and returns Z = inv(Q_reduced).  Note: for the two-area system this severs
    the inter-area coupling between CBRs 1/2 and CBRs 3/4 because the only path
    between them passes through the removed reference node 9.  Use
    build_network_z_pinv() for a topology-preserving Z matrix.
    """
    Q_full = build_laplacian(lines_path, n_nodes)
    keep = [k for k in range(n_nodes) if k != ground_node - 1]
    Q = Q_full[np.ix_(keep, keep)]
    return np.linalg.inv(Q)


def build_network_z_pinv(lines_path: Path, n_nodes: int) -> np.ndarray:
    """Z matrix via Moore-Penrose pseudo-inverse of the full Laplacian.

    For the two-area system the inter-area coupling (CBRs 1/2 ↔ CBRs 3/4)
    passes through the reference node 9.  Removing node 9 (as in
    build_network_z) breaks this coupling.  Using the pseudo-inverse of the
    full n_nodes × n_nodes Laplacian preserves all inter-node coupling while
    still producing a valid impedance matrix: current injections that sum to
    zero (as required by the gSCR definition with S'_B1 ⊇ CBR + SE) yield
    consistent voltage responses.
    """
    Q_full = build_laplacian(lines_path, n_nodes)
    return np.linalg.pinv(Q_full)


def reconstruct_impedance() -> tuple[np.ndarray, np.ndarray]:
    """Reconstruct a symmetric Z consistent with the printed matrix.

    Scope limitation: only the four CBR rows of Z are recovered from the
    printed W via Z[i,:] = W[i,:] / S_i.  Symmetry then fills columns 0-3
    for all rows, but the 5x5 passive-to-passive sub-block (rows 4-8,
    cols 4-8) remains zero because those entries do not appear in the paper.

    This is sufficient for the two greedy steps reproduced here, because
    the optimal nodes (1 and 4) are both CBR nodes.  If the algorithm were
    run to place more SEs at passive nodes the missing diagonal entries
    Z[j,j] for j >= 4 would cause incorrect sensitivity values.  For a
    full 9-node greedy run, build Z from the line data in
    matlab/data/table_x_two_area_lines.csv via build_network_from_edges.
    """

    capacities = np.array([0.5, 1.0, 1.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0])
    impedance = np.zeros_like(INITIAL_WEIGHTED)
    for row, capacity in enumerate(capacities[:4]):
        impedance[row, :] = INITIAL_WEIGHTED[row, :] / capacity

    for row in range(4):
        for col in range(9):
            impedance[col, row] = impedance[row, col]

    return capacities, impedance


def compute_table_iii(
    lines_path: Path,
    se_capacity: float = 0.5,
) -> tuple[list[tuple[int, int, float]], float]:
    """Enumerate gSCR for every possible 2-SE placement (Table III).

    Uses the full 8-node Z-bus built from the line CSV (groundNode=9 matches
    the MATLAB build_network_from_edges convention).  Node 9 is the infinite-bus
    reference and is not included as a candidate.

    The raw Z from the CSV is in the physical per-unit system of the line data.
    To match the paper's Table III scale, Z is calibrated so that the initial
    gSCR (no SEs, CBR capacities [0.5, 1.0, 1.5, 0.5]) equals the value given
    by the paper's printed INITIAL_WEIGHTED matrix.  This mirrors the calibration
    done for the 39-node system: k = lam_paper / lam_csv.

    Returns (results, scale_factor) where results is a list of
    (node_a, node_b, gscr) tuples in upper-triangular order (1-indexed).
    """
    # Use pseudo-inverse of full 9-node Laplacian so that inter-area
    # coupling (nodes 1/2 ↔ nodes 3/4 through node 9) is preserved.
    Z_raw = build_network_z_pinv(lines_path, n_nodes=9)   # 9x9

    # CBR capacities for all 9 nodes (paper nodes 1-9)
    s_base = np.array([0.5, 1.0, 1.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0])

    # Calibrate: scale Z so λmax(diag(s_base) @ Z) = λmax(INITIAL_WEIGHTED)
    # Both expressions produce the same eigenvalue because passive rows are zero
    # in both W matrices, so only the 4x4 CBR block matters.
    lam_csv   = lambda_max_positive(np.diag(s_base) @ Z_raw)
    lam_paper = lambda_max_positive(INITIAL_WEIGHTED)
    k = lam_paper / lam_csv
    Z = k * Z_raw

    results: list[tuple[int, int, float]] = []
    n = len(s_base)
    for i, j in combinations_with_replacement(range(n), 2):
        s = s_base.copy()
        s[i] -= se_capacity
        s[j] -= se_capacity
        W = np.diag(s) @ Z
        try:
            gscr = gscr_from_weighted_impedance(W)
        except ValueError:
            gscr = float("nan")
        results.append((i + 1, j + 1, gscr))   # 1-indexed

    return results, k


def main() -> None:
    RESULTS.mkdir(exist_ok=True)
    lines_path = ROOT / "matlab" / "data" / "table_x_two_area_lines.csv"

    # ── greedy placement (Tables I & II) ─────────────────────────────────
    # Build Z from network data via pseudo-inverse so that inter-area coupling
    # through node 9 is preserved, then calibrate to the paper's initial W scale.
    Z_raw = build_network_z_pinv(lines_path, n_nodes=9)
    s_base = np.array([0.5, 1.0, 1.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0])
    lam_csv   = lambda_max_positive(np.diag(s_base) @ Z_raw)
    lam_paper = lambda_max_positive(INITIAL_WEIGHTED)
    k_greedy  = lam_paper / lam_csv
    Z_calibrated = k_greedy * Z_raw

    initial_gscr = gscr_from_weighted_impedance(weighted_impedance(s_base, Z_calibrated))
    final_capacities, steps = greedy_place_storage(
        s_base,
        Z_calibrated,
        candidate_nodes=list(range(9)),
        storage_capacity=0.5,
        count=2,
    )
    final_weighted = weighted_impedance(final_capacities, Z_calibrated)
    final_gscr = gscr_from_weighted_impedance(final_weighted)

    greedy_lines = [
        "iteration,chosen_node,lambda_max_before,gscr_before",
        *[
            f"{step.iteration},{step.chosen_node + 1},{step.lambda_max:.6f},{step.gscr:.6f}"
            for step in steps
        ],
    ]
    (RESULTS / "two_area_placement.csv").write_text("\n".join(greedy_lines) + "\n")

    print(f"Initial gSCR: {initial_gscr:.4f}")
    for step in steps:
        print(
            f"Iteration {step.iteration}: choose node {step.chosen_node + 1}, "
            f"lambda_max={step.lambda_max:.4f}, gSCR={step.gscr:.4f}"
        )
    print(f"Final selected nodes: {[step.chosen_node + 1 for step in steps]}")
    print(f"Final gSCR after placement: {final_gscr:.4f}")

    # ── Table III: exhaustive search ─────────────────────────────────────
    print("\n── Table III: gSCR for all 2-SE placements (nodes 1-8) ──")
    table3, z_scale = compute_table_iii(lines_path)
    print(f"Z calibration scale factor k = {z_scale:.6f}")

    # Find optimal and compare with paper
    best = max(table3, key=lambda r: r[2])
    print(f"Computed optimal: nodes ({best[0]}, {best[1]})  gSCR = {best[2]:.4f}")
    print("Paper reports global optimum at nodes (1, 4).")

    # Load paper reference and compare
    ref_path = ROOT / "matlab" / "data" / "table_iii_two_area_two_se_gscr.csv"
    ref: dict[tuple[int, int], float] = {}
    if ref_path.exists():
        with open(ref_path, newline="") as fh:
            for row in csv.DictReader(fh):
                a, b = int(row["node_a"]), int(row["node_b"])
                ref[(a, b)] = float(row["gscr"])

    print(f"\n{'Nodes':>10}  {'Computed':>10}  {'Paper':>10}  {'Gap':>8}")
    print("-" * 45)
    max_gap = 0.0
    for a, b, g_comp in sorted(table3):
        g_ref = ref.get((a, b), float("nan"))
        gap = g_comp - g_ref if not np.isnan(g_ref) else float("nan")
        if not np.isnan(gap):
            max_gap = max(max_gap, abs(gap))
        marker = " ← greedy" if (a, b) == (1, 4) else ""
        print(f"  ({a},{b}):     {g_comp:10.4f}  {g_ref:10.4f}  {gap:8.4f}{marker}")

    print(f"\nMax |gap| vs paper Table III: {max_gap:.4f}")

    # Save computed Table III
    t3_lines = ["node_a,node_b,gscr_computed"]
    for a, b, g in table3:
        t3_lines.append(f"{a},{b},{g:.6f}")
    (RESULTS / "table_iii_computed.csv").write_text("\n".join(t3_lines) + "\n")
    print("Saved → results/table_iii_computed.csv")

    # ── bar chart ─────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 3.6))
    ax.bar(
        [step.chosen_node + 1 for step in steps],
        [step.gscr for step in steps],
        color=["#315f72", "#c77b30"],
        width=0.55,
    )
    ax.set_xlabel("Chosen SE node")
    ax.set_ylabel("gSCR before placement")
    ax.set_title("Two-area four-CBR greedy SE placement")
    ax.set_xticks([step.chosen_node + 1 for step in steps])
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(RESULTS / "two_area_gscr.png", dpi=180)


if __name__ == "__main__":
    main()
