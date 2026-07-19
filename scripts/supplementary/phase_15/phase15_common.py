#!/usr/bin/env python3
"""Shared utilities and a self-contained Bayesian sampler for Phase 15.

Model
-----
Observed turnover estimates y_i are treated as noisy summaries of latent pair-level
responses theta_i. Resampling quantiles from Phase 14B provide an approximate
observation-dispersion scale s_i.

    y_i | theta_i ~ Normal(theta_i, s_i)
    theta_i | beta, u_cell, sigma, lambda_i ~ Normal(X_i beta + u_cell[i], sigma^2/lambda_i)
    lambda_i ~ Gamma(nu/2, nu/2), fixed nu=4
    u_cell ~ Normal(0, tau^2)

The Gamma scale-mixture produces a Student-t process likelihood. Regression
coefficients have zero-centered Normal priors; sigma and tau have half-Normal
priors sampled by Metropolis-within-Gibbs. The implementation uses only NumPy.
"""
from __future__ import annotations

import csv
import itertools
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np

SCRIPT_VERSION = "phase15_common_v0.1.0_2026-07-18"


def analysis_root(project_root: Path) -> Path:
    preferred = project_root / "04_analysis_USE _THIS"
    if preferred.exists():
        return preferred / "14_recent_environmental_stress"
    return project_root / "04_analysis" / "14_recent_environmental_stress"


def phase15_root(project_root: Path) -> Path:
    preferred = project_root / "04_analysis_USE _THIS"
    if preferred.exists():
        return preferred / "15_bayesian_h3_evidence"
    return project_root / "04_analysis" / "15_bayesian_h3_evidence"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    rows = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("\n", encoding="utf-8")
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")


def as_float(value: Any) -> float | None:
    try:
        result = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def zscore(values: np.ndarray, mean: float | None = None, sd: float | None = None) -> tuple[np.ndarray, float, float]:
    if mean is None:
        mean = float(np.mean(values))
    if sd is None:
        sd = float(np.std(values, ddof=1))
    if not math.isfinite(sd) or sd <= 0:
        raise ValueError("Cannot standardize a constant or invalid variable.")
    return (values - mean) / sd, float(mean), float(sd)


def period_midpoint(period_1: str, period_2: str) -> float:
    years = [int(x) for x in re.findall(r"(?:19|20)\d{2}", f"{period_1} {period_2}")]
    if len(years) < 4:
        raise ValueError(f"Cannot recover transition years from {period_1!r}, {period_2!r}")
    return (years[0] + years[-1]) / 2.0


@dataclass
class PreparedData:
    rows: list[dict[str, Any]]
    y: np.ndarray
    se: np.ndarray
    stress_raw: np.ndarray
    log_events_raw: np.ndarray
    midpoint_raw: np.ndarray
    cells: list[str]
    cell_index: np.ndarray
    unique_cells: list[str]
    pair_ids: list[str]


def load_prepared(path: Path) -> PreparedData:
    rows = read_csv(path)
    if not rows:
        raise RuntimeError(f"No rows found in {path}")
    y = np.array([float(row["y_observed"]) for row in rows], dtype=float)
    se = np.array([float(row["observation_sd_approx"]) for row in rows], dtype=float)
    stress = np.array([float(row["stress_change_raw"]) for row in rows], dtype=float)
    log_events = np.array([float(row["log_common_events_raw"]) for row in rows], dtype=float)
    midpoint = np.array([float(row["transition_midpoint_raw"]) for row in rows], dtype=float)
    cells = [row["grid_cell_id"] for row in rows]
    unique_cells = sorted(set(cells))
    cell_map = {cell: i for i, cell in enumerate(unique_cells)}
    cell_index = np.array([cell_map[cell] for cell in cells], dtype=int)
    return PreparedData(
        rows=rows, y=y, se=se, stress_raw=stress, log_events_raw=log_events,
        midpoint_raw=midpoint, cells=cells, cell_index=cell_index,
        unique_cells=unique_cells, pair_ids=[row["pair_id"] for row in rows],
    )


def build_design(data: PreparedData, include_stress: bool = True,
                 train_mask: np.ndarray | None = None) -> tuple[np.ndarray, list[str], dict[str, float]]:
    n = len(data.y)
    if train_mask is None:
        train_mask = np.ones(n, dtype=bool)
    columns = [np.ones(n)]
    names = ["intercept"]
    scaling: dict[str, float] = {}
    if include_stress:
        z, mean, sd = zscore(data.stress_raw[train_mask])
        all_z = (data.stress_raw - mean) / sd
        columns.append(all_z)
        names.append("stress_beta")
        scaling.update(stress_mean=mean, stress_sd=sd)
    for raw, base in ((data.log_events_raw, "log_events"), (data.midpoint_raw, "midpoint")):
        _, mean, sd = zscore(raw[train_mask])
        columns.append((raw - mean) / sd)
        names.append(f"{base}_beta")
        scaling.update({f"{base}_mean": mean, f"{base}_sd": sd})
    return np.column_stack(columns), names, scaling


def log_halfnormal_on_log_scale(log_value: float, scale: float) -> float:
    value = math.exp(log_value)
    return -0.5 * (value / scale) ** 2 + log_value


def sample_mvn_precision(rng: np.random.Generator, precision: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    chol = np.linalg.cholesky(precision)
    mean = np.linalg.solve(chol.T, np.linalg.solve(chol, rhs))
    z = rng.normal(size=len(rhs))
    noise = np.linalg.solve(chol.T, z)
    return mean + noise


def _sigma_target(log_sigma: float, residual: np.ndarray, lam: np.ndarray, prior_scale: float) -> float:
    sigma = math.exp(log_sigma)
    return -len(residual) * log_sigma - 0.5 * float(np.sum(lam * residual * residual)) / (sigma * sigma) + log_halfnormal_on_log_scale(log_sigma, prior_scale)


def _tau_target(log_tau: float, theta: np.ndarray, xb: np.ndarray, z_cell: np.ndarray,
                cell_index: np.ndarray, lam: np.ndarray, sigma: float, prior_scale: float) -> float:
    # Non-centered cell effect u_j = tau * z_j. The standard-normal prior on z
    # does not depend on tau; tau is informed through the process likelihood.
    tau = math.exp(log_tau)
    residual = theta - xb - tau * z_cell[cell_index]
    return -0.5 * float(np.sum(lam * residual * residual)) / (sigma * sigma) + log_halfnormal_on_log_scale(log_tau, prior_scale)



def slice_sample_log(rng: np.random.Generator, current: float, logpdf, width: float = 0.5, max_steps: int = 100) -> float:
    """Univariate stepping-out slice sampler on an unconstrained log scale."""
    log_height = float(logpdf(current)) - rng.exponential(1.0)
    left = current - rng.uniform() * width
    right = left + width
    j = int(rng.integers(0, max_steps + 1))
    k = max_steps - j
    while j > 0 and float(logpdf(left)) > log_height:
        left -= width
        j -= 1
    while k > 0 and float(logpdf(right)) > log_height:
        right += width
        k -= 1
    for _ in range(1000):
        proposal = rng.uniform(left, right)
        if float(logpdf(proposal)) >= log_height:
            return float(proposal)
        if proposal < current:
            left = proposal
        else:
            right = proposal
    raise RuntimeError("Slice sampler failed to find an acceptable point.")


def run_chain(
    y: np.ndarray, se: np.ndarray, X: np.ndarray, cell_index: np.ndarray,
    prior_sd: np.ndarray, draws: int, warmup: int, thin: int, seed: int,
    nu: float = 4.0, sigma_prior_scale: float = 0.25,
    tau_prior_scale: float = 0.15,
) -> dict[str, np.ndarray]:
    rng = np.random.default_rng(seed)
    n, p = X.shape
    j_count = int(cell_index.max()) + 1
    beta = np.zeros(p)
    try:
        beta = np.linalg.lstsq(X, y, rcond=None)[0]
    except np.linalg.LinAlgError:
        pass
    z_cell = np.zeros(j_count)
    u = np.zeros(j_count)
    theta = y.copy()
    sigma = max(0.05, float(np.std(y - X @ beta, ddof=max(1, p))))
    tau = min(0.15, max(0.03, sigma / 2))
    lam = np.ones(n)
    log_sigma = math.log(sigma)
    log_tau = math.log(tau)
    total = warmup + draws * thin
    stored_beta: list[np.ndarray] = []
    stored_sigma: list[float] = []
    stored_tau: list[float] = []
    stored_u: list[np.ndarray] = []
    se2 = np.maximum(se * se, 1e-8)
    prior_precision = np.diag(1.0 / np.maximum(prior_sd * prior_sd, 1e-12))

    for it in range(total):
        mu = X @ beta + u[cell_index]
        precision_theta = 1.0 / se2 + lam / (sigma * sigma)
        mean_theta = (y / se2 + lam * mu / (sigma * sigma)) / precision_theta
        theta = mean_theta + rng.normal(size=n) / np.sqrt(precision_theta)

        weights = lam / (sigma * sigma)
        precision_beta = X.T @ (weights[:, None] * X) + prior_precision
        rhs_beta = X.T @ (weights * (theta - u[cell_index]))
        beta = sample_mvn_precision(rng, precision_beta, rhs_beta)

        xb = X @ beta
        # Non-centered cell effects: u_j = tau * z_j, z_j ~ Normal(0, 1).
        for j in range(j_count):
            idx = np.where(cell_index == j)[0]
            prec = 1.0 + (tau * tau) * float(np.sum(weights[idx]))
            mean_z = tau * float(np.sum(weights[idx] * (theta[idx] - xb[idx]))) / prec
            z_cell[j] = mean_z + rng.normal() / math.sqrt(prec)
        u = tau * z_cell

        residual = theta - xb - u[cell_index]
        shape = (nu + 1.0) / 2.0
        rate = (nu + (residual / sigma) ** 2) / 2.0
        lam = rng.gamma(shape=shape, scale=1.0 / rate)

        log_sigma = slice_sample_log(
            rng, log_sigma,
            lambda value: _sigma_target(value, residual, lam, sigma_prior_scale),
            width=0.45,
        )
        sigma = math.exp(log_sigma)

        log_tau = slice_sample_log(
            rng, log_tau,
            lambda value: _tau_target(value, theta, xb, z_cell, cell_index, lam, sigma, tau_prior_scale),
            width=0.55,
        )
        tau = math.exp(log_tau)
        u = tau * z_cell

        if it >= warmup and (it - warmup) % thin == 0:
            stored_beta.append(beta.copy())
            stored_sigma.append(sigma)
            stored_tau.append(tau)
            stored_u.append(u.copy())

    return {
        "beta": np.asarray(stored_beta),
        "sigma": np.asarray(stored_sigma),
        "tau": np.asarray(stored_tau),
        "u": np.asarray(stored_u),
    }


def split_rhat(chains: np.ndarray) -> float:
    # chains shape: chain x draw
    m, n = chains.shape
    half = n // 2
    if m < 2 or half < 10:
        return float("nan")
    split = np.concatenate([chains[:, :half], chains[:, -half:]], axis=0)
    n = half
    chain_means = split.mean(axis=1)
    W = float(np.mean(np.var(split, axis=1, ddof=1)))
    B = float(n * np.var(chain_means, ddof=1))
    if W <= 0:
        return 1.0 if B <= 0 else float("inf")
    var_hat = (n - 1) / n * W + B / n
    return math.sqrt(var_hat / W)


def effective_sample_size(chains: np.ndarray) -> float:
    m, n = chains.shape
    if m < 2 or n < 20:
        return float("nan")
    chain_means = chains.mean(axis=1)
    W = float(np.mean(np.var(chains, axis=1, ddof=1)))
    B = float(n * np.var(chain_means, ddof=1))
    var_plus = (n - 1) / n * W + B / n
    if var_plus <= 0:
        return float(m * n)
    centered = chains - chain_means[:, None]
    rhos = []
    for lag in range(1, min(n - 1, 2000)):
        autocov = float(np.mean(np.sum(centered[:, :-lag] * centered[:, lag:], axis=1) / (n - lag)))
        rho = 1.0 - (W - autocov) / var_plus
        rhos.append(rho)
        if lag >= 2 and lag % 2 == 0 and rhos[-1] + rhos[-2] < 0:
            rhos = rhos[:-2]
            break
    # initial monotone paired sequence
    paired = []
    for i in range(0, len(rhos) - 1, 2):
        value = rhos[i] + rhos[i + 1]
        if paired:
            value = min(value, paired[-1])
        if value < 0:
            break
        paired.append(value)
    tau_hat = max(1.0, -1.0 + 2.0 * sum(1.0 + value for value in paired)) if paired else 1.0
    return min(float(m * n), float(m * n) / tau_hat)


def summarize_chains(chain_outputs: list[dict[str, np.ndarray]], names: list[str]) -> tuple[list[dict[str, Any]], dict[str, np.ndarray]]:
    arrays: dict[str, np.ndarray] = {}
    for k, name in enumerate(names):
        arrays[name] = np.stack([chain["beta"][:, k] for chain in chain_outputs])
    arrays["process_sd"] = np.stack([chain["sigma"] for chain in chain_outputs])
    arrays["cell_sd"] = np.stack([chain["tau"] for chain in chain_outputs])
    rows = []
    for name, arr in arrays.items():
        flat = arr.reshape(-1)
        rows.append({
            "parameter": name,
            "mean": float(np.mean(flat)),
            "sd": float(np.std(flat, ddof=1)),
            "median": float(np.median(flat)),
            "q025": float(np.quantile(flat, 0.025)),
            "q055": float(np.quantile(flat, 0.055)),
            "q945": float(np.quantile(flat, 0.945)),
            "q975": float(np.quantile(flat, 0.975)),
            "rhat": split_rhat(arr),
            "ess_bulk_approx": effective_sample_size(arr),
        })
    return rows, arrays


def posterior_predictive(
    rng: np.random.Generator, chain_outputs: list[dict[str, np.ndarray]], X: np.ndarray,
    cell_index: np.ndarray, se: np.ndarray, nu: float = 4.0, max_draws: int = 3000,
    new_cell: bool = False,
) -> np.ndarray:
    beta = np.concatenate([c["beta"] for c in chain_outputs], axis=0)
    sigma = np.concatenate([c["sigma"] for c in chain_outputs])
    tau = np.concatenate([c["tau"] for c in chain_outputs])
    u = np.concatenate([c["u"] for c in chain_outputs], axis=0)
    if len(beta) > max_draws:
        idx = rng.choice(len(beta), max_draws, replace=False)
        beta, sigma, tau, u = beta[idx], sigma[idx], tau[idx], u[idx]
    result = np.empty((len(beta), X.shape[0]), dtype=float)
    for d in range(len(beta)):
        if new_cell:
            cell_effects = rng.normal(0.0, tau[d], size=X.shape[0])
        else:
            cell_effects = u[d, cell_index]
        mu = X @ beta[d] + cell_effects
        process = rng.standard_t(nu, size=X.shape[0]) * sigma[d]
        obs = rng.normal(0.0, se, size=X.shape[0])
        result[d] = mu + process + obs
    return result


def empirical_crps(samples: np.ndarray, observed: float) -> float:
    samples = np.sort(np.asarray(samples, dtype=float))
    n = len(samples)
    term1 = float(np.mean(np.abs(samples - observed)))
    weights = 2.0 * np.arange(1, n + 1) - n - 1.0
    pair_term = float(np.sum(weights * samples)) * 2.0 / (n * n)
    return term1 - 0.5 * pair_term


def exact_sign_flip_pvalue(values: np.ndarray) -> tuple[float, int]:
    values = np.asarray(values, dtype=float)
    observed = abs(float(np.mean(values)))
    n = len(values)
    extreme = 0
    total = 0
    for signs in itertools.product((-1.0, 1.0), repeat=n):
        statistic = abs(float(np.mean(values * np.asarray(signs))))
        total += 1
        if statistic >= observed - 1e-12:
            extreme += 1
    return extreme / total, total


def prior_sd_vector(names: list[str], beta_scale: float, intercept_scale: float = 0.25, control_scale: float = 0.20) -> np.ndarray:
    values = []
    for name in names:
        if name == "intercept":
            values.append(intercept_scale)
        elif name == "stress_beta":
            values.append(beta_scale)
        else:
            values.append(control_scale)
    return np.asarray(values, dtype=float)
