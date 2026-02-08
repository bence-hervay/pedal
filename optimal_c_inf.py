import mpmath as mp


def a_coeff(n):
    # a_n = (C(2n,n)^2) / 16^n
    return mp.binomial(2 * n, n) ** 2 / (16**n)


def F_limit(c, J=80, N=20, a=None):
    """
    Limit objective F(c) whose minimizer is c_infty.

    J: number of 'exact' terms in the finite sum
    N: number of K-series coefficients used in the tail
    """
    if a is None:
        a = [a_coeff(n) for n in range(N + 1)]

    beta = (1 + mp.sin(c)) / 2

    # T0 = K(beta)
    F = mp.ellipk(beta)

    # Finite exact part: sum_{k=1..J} [ΔT_k - π/(2√k)]
    for k in range(1, J + 1):
        F += mp.ellipk(1 / (k + beta)) / mp.sqrt(k + beta) - mp.pi / (2 * mp.sqrt(k))

    # Tail (k > J), accelerated with Hurwitz zeta
    # n=0 piece corresponds to (k+beta)^(-1/2); we need the difference to k^(-1/2)
    tail = mp.zeta(mp.mpf("0.5"), J + 1 + beta) - mp.zeta(mp.mpf("0.5"), J + 1)

    # n>=1 pieces from K-series: a_n (k+beta)^(-(n+1/2))
    for n in range(1, N + 1):
        s = n + mp.mpf("0.5")
        tail += a[n] * mp.zeta(s, J + 1 + beta)

    F += (mp.pi / 2) * tail
    return F


def c_infty(dps=80, J=120, N=20):
    mp.mp.dps = dps
    a = [a_coeff(n) for n in range(N + 1)]
    F = lambda cc: F_limit(cc, J=J, N=N, a=a)

    # Coarse scan to get a robust initial guess
    lo = mp.mpf("0.2")
    hi = mp.mpf("0.8")
    grid = [lo + (hi - lo) * i / 200 for i in range(201)]
    vals = [F(c) for c in grid]
    i0 = min(range(len(vals)), key=lambda i: vals[i])
    c0 = grid[i0]

    # Solve F'(c)=0 near the minimizer (use a bracketed pair)
    dF = lambda cc: mp.diff(F, cc)
    return mp.findroot(dF, (c0 - mp.mpf("0.05"), c0 + mp.mpf("0.05")))


def main():
    c = c_infty(dps=100, J=200, N=30)
    print("c_infty =", c)
    print("deg     =", c * 180 / mp.pi)


if __name__ == "__main__":
    main()
