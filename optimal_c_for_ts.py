#!/usr/bin/env python3
"""
Compute and plot c*(t) = argmax_c x(t,c) for the ODE:
    x(0)=0, x'(0)=0, x''(t) = |cos(x(t)-c)|,   with t>=0 and |c|<=pi/2.

This script produces:
  - an interactive Plotly HTML (thin blue line + red optimized-point markers)
  - a static matplotlib PNG

Hyperparameters are grouped in CONFIG below.

Dependencies:
  numpy, scipy, plotly, matplotlib
"""

import math, time, random
import numpy as np
from scipy import special, optimize
import plotly.graph_objects as go
import matplotlib.pyplot as plt


# -----------------------
# CONFIG (edit these)
# -----------------------
CONFIG = dict(
    # t-values to solve for (you can replace this with any numpy array you want)
    t_max=100.0,
    n_t_points=20,  # used only if you use the default quadratic spacing below
    t_spacing="quadratic",  # "quadratic" or "linear"
    # If you want to auto-pick n_t_points based on a timing run, set this:
    auto_resolution=True,
    time_budget_seconds=120.0,
    n_t_min=80,
    n_t_max=350,
    timing_t=100.0,  # time a single argmax at this t (usually the max t)
    # Argmax hyperparameters
    n_grid=300,  # grid scan resolution over c in [-pi/2, pi/2]
    n_rand=100,  # number of random multi-start refinement intervals
    xatol=1e-9,  # tolerance for local bounded optimizer (minimize_scalar)
    # Speed/accuracy tradeoff in x(t,c) evaluation:
    # segment locator uses exact ellipk for k < k_switch, series for k>=k_switch
    k_switch=10,
    # Output files
    out_html="optimal_c_vs_t_100.html",
    out_png="optimal_c_vs_t_100.png",
)


# -----------------------
# Fast elliptic helpers
# -----------------------
PI = math.pi
PI_OVER_2 = math.pi / 2.0


def ellipk_series(m: float) -> float:
    """Series approximation of K(m) for small m (very accurate for m<=0.1)."""
    m2 = m * m
    m3 = m2 * m
    m4 = m2 * m2
    m5 = m4 * m
    return PI_OVER_2 * (
        1.0
        + 0.25 * m
        + (9.0 / 64.0) * m2
        + (25.0 / 256.0) * m3
        + (1225.0 / 16384.0) * m4
        + (3969.0 / 65536.0) * m5
    )


def x_exact_piecewise_fast(t: float, c: float, k_switch: int = 10) -> float:
    """
    Closed-form x(t,c) using the piecewise Jacobi-elliptic solution,
    with a fast segment locator for large t.
    """
    if abs(c) >= PI_OVER_2 - 1e-14:
        return 0.0

    beta = 0.5 * (1.0 + math.sin(c))  # in (0,1)
    K0 = float(special.ellipk(beta))

    # Segment 0 (u in [-pi/2, pi/2])
    if t <= K0:
        u_arg = K0 - t
        sn = float(special.ellipj(u_arg, beta)[0])
        s = math.sqrt(beta) * sn
        s = 1.0 if s > 1.0 else (-1.0 if s < -1.0 else s)
        u = PI_OVER_2 - 2.0 * math.asin(s)
        return u + c

    # Later segments k>=1
    t_rem = t - K0
    acc_time = 0.0
    k = 1
    while True:
        denom = math.sqrt(k + beta)
        m = 1.0 / (k + beta)
        if k < k_switch:
            Kk = float(special.ellipk(m))
        else:
            Kk = ellipk_series(m)
        dt_seg = Kk / denom

        if acc_time + dt_seg >= t_rem - 1e-15:
            tau = t_rem - acc_time
            mk = 1.0 / (k + beta)
            K_exact = float(special.ellipk(mk))
            u_arg = K_exact - math.sqrt(k + beta) * tau
            sn = float(special.ellipj(u_arg, mk)[0])
            sn = 1.0 if sn > 1.0 else (-1.0 if sn < -1.0 else sn)
            u = k * PI + PI_OVER_2 - 2.0 * math.asin(sn)
            return u + c

        acc_time += dt_seg
        k += 1


def argmax_c_for_t(
    t: float, *, n_grid: int, n_rand: int, xatol: float, seed: int, k_switch: int
):
    """
    Independent argmax over c in [-pi/2, pi/2]:
      - grid scan
      - refine around local peaks on the grid
      - many random refinement intervals (multi-start)
      - bounded 1D local optimization in each interval
    """
    if t <= 1e-14:
        return 0.0, 0.0

    rng = random.Random(seed)
    lo, hi = -PI_OVER_2, PI_OVER_2

    cs_grid = np.linspace(lo, hi, n_grid)
    fs_grid = np.array(
        [
            x_exact_piecewise_fast(float(t), float(c), k_switch=k_switch)
            for c in cs_grid
        ],
        dtype=float,
    )

    intervals = [(lo, float(cs_grid[1])), (float(cs_grid[-2]), hi)]
    for i in range(1, n_grid - 1):
        if fs_grid[i] >= fs_grid[i - 1] and fs_grid[i] >= fs_grid[i + 1]:
            intervals.append((float(cs_grid[i - 1]), float(cs_grid[i + 1])))

    step = (hi - lo) / (n_grid - 1)
    w = 2.5 * step
    for _ in range(n_rand):
        c0 = rng.uniform(lo, hi)
        a = max(lo, c0 - w)
        b = min(hi, c0 + w)
        if b - a > 1e-12:
            intervals.append((a, b))

    best_c = None
    best_f = -1e300

    def neg_f(cc):
        return -x_exact_piecewise_fast(float(t), float(cc), k_switch=k_switch)

    for a, b in intervals:
        if b - a < 1e-12:
            continue
        res = optimize.minimize_scalar(
            neg_f, bounds=(a, b), method="bounded", options={"xatol": xatol}
        )
        c_hat = float(res.x)
        f_hat = x_exact_piecewise_fast(float(t), c_hat, k_switch=k_switch)
        if f_hat > best_f:
            best_f, best_c = f_hat, c_hat

    return best_c, best_f


def make_t_values(cfg):
    """Define the set of t values to optimize at."""
    t_max = float(cfg["t_max"])
    n = int(cfg["n_t_points"])
    if cfg["t_spacing"] == "quadratic":
        ts = t_max * (np.linspace(0.0, 1.0, n) ** 2)
    elif cfg["t_spacing"] == "linear":
        ts = np.linspace(0.0, t_max, n)
    else:
        raise ValueError("t_spacing must be 'quadratic' or 'linear'")
    ts[0] = 0.0
    ts[-1] = t_max
    return ts


def main():
    cfg = CONFIG.copy()

    # Optionally auto-select number of t points based on a timing run
    if cfg["auto_resolution"]:
        t_timing = float(cfg["timing_t"])
        seed = 123456
        t0 = time.perf_counter()
        _ = argmax_c_for_t(
            t_timing,
            n_grid=cfg["n_grid"],
            n_rand=cfg["n_rand"],
            xatol=cfg["xatol"],
            seed=seed,
            k_switch=cfg["k_switch"],
        )
        dt = time.perf_counter() - t0
        n_est = int(float(cfg["time_budget_seconds"]) / max(dt, 1e-6))
        n_est = max(int(cfg["n_t_min"]), min(int(cfg["n_t_max"]), n_est))
        cfg["n_t_points"] = n_est
        print(
            f"[auto_resolution] timing at t={t_timing}: {dt:.4f} s -> n_t_points={n_est}"
        )

    ts = make_t_values(cfg)

    cs_star = np.empty_like(ts)
    xs_star = np.empty_like(ts)

    for i, t in enumerate(ts):
        seed_i = 20_000_000 + i * 10007  # independent seed per t
        c_star, x_star = argmax_c_for_t(
            float(t),
            n_grid=cfg["n_grid"],
            n_rand=cfg["n_rand"],
            xatol=cfg["xatol"],
            seed=seed_i,
            k_switch=cfg["k_switch"],
        )
        cs_star[i] = c_star
        xs_star[i] = x_star
        if (i + 1) % 50 == 0:
            print(f"Optimized {i + 1}/{len(ts)} points...")

    c_min = float(np.min(cs_star))
    c_max = float(np.max(cs_star))
    rng = c_max - c_min
    margin = 0.1 if rng < 1e-12 else 0.10 * rng
    y_lo = max(-PI_OVER_2, c_min - margin)
    y_hi = min(PI_OVER_2, c_max + margin)

    # Plotly interactive
    hover_text = [
        f"t={t:.6f}<br>c*={c:.9f} rad<br>x(t,c*)={xv:.6f}"
        for t, c, xv in zip(ts, cs_star, xs_star)
    ]

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=ts,
            y=cs_star,
            mode="lines",
            line=dict(color="blue", width=1),
            name="c*(t) (connected)",
            hoverinfo="skip",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=ts,
            y=cs_star,
            mode="markers",
            marker=dict(color="red", size=5),
            name="Optimized points",
            text=hover_text,
            hoverinfo="text",
        )
    )
    fig.update_layout(
        title="argmax_c x(t,c): optimal c* vs t",
        xaxis_title="t",
        yaxis_title="c* (radians)",
        yaxis=dict(range=[y_lo, y_hi]),
        template="plotly_white",
    )
    fig.write_html(cfg["out_html"], include_plotlyjs="cdn")
    print(f"Wrote interactive plot: {cfg['out_html']}")

    # Static matplotlib
    plt.figure(figsize=(9.5, 4.8))
    plt.plot(ts, cs_star, linewidth=0.75, color="blue")
    plt.scatter(ts, cs_star, s=8, color="red")
    plt.ylim([y_lo, y_hi])
    plt.xlabel("t")
    plt.ylabel("c* (radians)")
    plt.title("argmax_c x(t,c): optimal c* vs t")
    plt.tight_layout()
    plt.savefig(cfg["out_png"], dpi=160)
    plt.close()
    print(f"Wrote static plot: {cfg['out_png']}")


if __name__ == "__main__":
    main()
