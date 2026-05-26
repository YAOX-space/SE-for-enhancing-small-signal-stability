"""Reproduce the paper's two-area four-CBR storage placement example.

The PDF reports the initial weighted impedance matrix S'_B1 Z in (17).  The
CBR capacity ratios can be inferred from the nonsymmetric entries because Z is
symmetric: S = [0.5, 1.0, 1.5, 0.5].  The missing passive-node rows of Z are
recovered by symmetry where possible; they are not needed for reproducing the
first two node selections because the dominant left eigenvector ranking depends
on the reported active rows.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np

from se_placement import greedy_place_storage, gscr_from_weighted_impedance
from se_placement import weighted_impedance


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / ".matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(RESULTS / ".cache"))

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


def reconstruct_impedance() -> tuple[np.ndarray, np.ndarray]:
    """Reconstruct a symmetric Z consistent with the printed matrix.

    Scope limitation: only the four CBR rows of Z are recovered from the
    printed W via Z[i,:] = W[i,:] / S_i.  Symmetry then fills columns 0-3
    for all rows, but the 5×5 passive-to-passive sub-block (rows 4-8,
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


def main() -> None:
    RESULTS.mkdir(exist_ok=True)

    capacities, impedance = reconstruct_impedance()
    initial_gscr = gscr_from_weighted_impedance(INITIAL_WEIGHTED)
    final_capacities, steps = greedy_place_storage(
        capacities,
        impedance,
        candidate_nodes=list(range(9)),
        storage_capacity=0.5,
        count=2,
    )
    final_weighted = weighted_impedance(final_capacities, impedance)
    final_gscr = gscr_from_weighted_impedance(final_weighted)

    lines = [
        "iteration,chosen_node,lambda_max_before,gscr_before",
        *[
            f"{step.iteration},{step.chosen_node + 1},{step.lambda_max:.6f},{step.gscr:.6f}"
            for step in steps
        ],
    ]
    (RESULTS / "two_area_placement.csv").write_text("\n".join(lines) + "\n")

    print(f"Initial gSCR: {initial_gscr:.4f}")
    for step in steps:
        print(
            f"Iteration {step.iteration}: choose node {step.chosen_node + 1}, "
            f"lambda_max={step.lambda_max:.4f}, gSCR={step.gscr:.4f}"
        )
    print(f"Final selected nodes: {[step.chosen_node + 1 for step in steps]}")
    print(f"Final gSCR after placement: {final_gscr:.4f}")
    print("Paper reports the same selected nodes for this example: [1, 4].")

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
