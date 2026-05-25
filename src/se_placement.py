"""Grid-strength based storage-energy placement.

This module implements the gSCR and greedy placement method described in
Yuan et al., "Placing Storage Energies for Enhancing Small-Signal Stability
of Converter-Based-Renewable Systems", IEEE TIA, 2025.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class PlacementStep:
    """One greedy storage-placement iteration."""

    iteration: int
    chosen_node: int
    lambda_max: float
    gscr: float
    sensitivities: dict[int, float]


def positive_eigenvalues(matrix: np.ndarray, tol: float = 1e-9) -> np.ndarray:
    """Return real positive eigenvalues of a numerically real matrix."""

    eigvals = np.linalg.eigvals(matrix)
    real_eigs = eigvals[np.abs(eigvals.imag) < 1e-7].real
    return np.sort(real_eigs[real_eigs > tol])


def lambda_max_positive(matrix: np.ndarray, tol: float = 1e-9) -> float:
    """Largest positive eigenvalue, the denominator of gSCR in the paper."""

    positive = positive_eigenvalues(matrix, tol=tol)
    if positive.size == 0:
        raise ValueError("matrix has no positive real eigenvalue")
    return float(positive[-1])


def gscr_from_weighted_impedance(weighted_impedance: np.ndarray) -> float:
    """Compute gSCR = 1 / lambda_max_positive(S'_B1 Z)."""

    return 1.0 / lambda_max_positive(weighted_impedance)


def weighted_impedance(signed_capacities: np.ndarray, impedance: np.ndarray) -> np.ndarray:
    """Build S'_B1 Z.

    `signed_capacities` is positive for CBR injection and negative for SE
    absorption, matching S_B1 = diag(S_B, -S_BE) in the paper.
    """

    signed_capacities = np.asarray(signed_capacities, dtype=float)
    impedance = np.asarray(impedance, dtype=float)
    if impedance.shape[0] != impedance.shape[1]:
        raise ValueError("impedance must be square")
    if impedance.shape[0] != signed_capacities.size:
        raise ValueError("capacity vector and impedance size do not match")
    return np.diag(signed_capacities) @ impedance


def left_eigen_sensitivities(
    weighted: np.ndarray,
    candidate_nodes: list[int],
) -> dict[int, float]:
    """Rank placement nodes using the paper's eigenvector sensitivity.

    For the largest positive eigenvalue lambda of S'_B1 Z, the paper derives
    d(lambda)/d(S_BE,j) = -a * lambda * v_j^2.  The unknown positive scaling
    `a` does not affect ranking, so this implementation reports
    -lambda * v_j^2 from the corresponding left eigenvector.

    Nodes are zero-based.
    """

    lam = lambda_max_positive(weighted)
    eigvals, eigvecs = np.linalg.eig(weighted.T)
    idx = int(np.argmin(np.abs(eigvals - lam)))
    vector = eigvecs[:, idx].real
    norm = np.linalg.norm(vector)
    if norm == 0:
        raise ValueError("degenerate eigenvector for sensitivity calculation")
    vector = vector / norm
    return {node: float(-lam * vector[node] ** 2) for node in candidate_nodes}


def greedy_place_storage(
    base_signed_capacities: np.ndarray,
    impedance: np.ndarray,
    candidate_nodes: list[int],
    storage_capacity: float,
    count: int,
) -> tuple[np.ndarray, list[PlacementStep]]:
    """Greedily place `count` identical SEs by minimizing lambda sensitivity."""

    capacities = np.array(base_signed_capacities, dtype=float)
    candidates = list(candidate_nodes)
    steps: list[PlacementStep] = []

    for iteration in range(1, count + 1):
        weighted = weighted_impedance(capacities, impedance)
        lam = lambda_max_positive(weighted)
        sensitivities = left_eigen_sensitivities(weighted, candidates)
        chosen = min(sensitivities, key=sensitivities.get)
        capacities[chosen] -= storage_capacity
        steps.append(
            PlacementStep(
                iteration=iteration,
                chosen_node=chosen,
                lambda_max=lam,
                gscr=1.0 / lam,
                sensitivities=sensitivities,
            )
        )

    return capacities, steps
