import mpmath as mp


# ---------- Elliptic K-series coefficient a_n ----------
def a_coeff(n: int) -> mp.mpf:
    # a_n = binomial(2n,n)^2 / 16^n
    return mp.binomial(2 * n, n) ** 2 / mp.power(16, n)


def precompute_a(N: int):
    return [a_coeff(n) for n in range(N)]


# ---------- dK/dm for complete elliptic K(m) (parameter m) ----------
def dK_dm(m: mp.mpf) -> mp.mpf:
    # d/dm ellipk(m) = (E(m) - (1-m)K(m)) / (2 m (1-m))
    K = mp.ellipk(m)
    E = mp.ellipe(m)
    return (E - (1 - m) * K) / (2 * m * (1 - m))


# ---------- derivative wrt s of g(s) = K(1/s)/sqrt(s) ----------
def gprime_term(s: mp.mpf) -> mp.mpf:
    m = 1 / s
    K = mp.ellipk(m)
    dK = dK_dm(m)
    return (-mp.mpf("0.5")) * s ** (-mp.mpf("1.5")) * K - s ** (-mp.mpf("2.5")) * dK


# ---------- H'(beta) with split + Hurwitz-zeta tail, plus an explicit tail bound ----------
def Hprime_beta(beta: mp.mpf, M: int, L: int, A) -> tuple[mp.mpf, mp.mpf]:
    """
    Returns (approx_value, rigorous_tail_bound) for H'(beta).

    H(beta) = K(beta) + sum_{k>=1} [ K(1/(k+beta))/sqrt(k+beta) - pi/(2sqrt(k)) ].
    The subtraction term vanishes under d/d beta, so:

      H'(beta) = dK/dm(beta) + sum_{k>=1} d/d beta [ K(1/(k+beta))/sqrt(k+beta) ].

    We split k=1..M directly; k>M via K-series and Hurwitz zeta:
      tail = -(pi/2) * sum_{n>=0} a_n (n+1/2) zeta(n+3/2, M+1+beta).

    We truncate n=0..L-1 and bound the remainder by a geometric bound.
    """
    pi = mp.pi
    # finite part
    val = dK_dm(beta)
    for k in range(1, M + 1):
        val += gprime_term(k + beta)

    # tail part (truncated)
    q = M + 1 + beta
    tail = mp.mpf("0")
    for n in range(L):
        s = n + mp.mpf("1.5")  # n + 3/2
        tail += A[n] * (n + mp.mpf("0.5")) * mp.zeta(s, q)
    val += (-pi / 2) * tail

    # explicit remainder bound for n >= L
    # Use a_n <= a_L for n>=L (monotone decreasing), and
    # zeta(s,q) <= q^{-s} + ∫_0^∞ (x+q)^{-s} dx = q^{-s} + q^{1-s}/(s-1)
    # with s=n+3/2, s-1 = n+1/2
    aL = A[L]
    r = 1 / q

    # sum_{n=L}∞ (n+1/2) q^{-n-3/2}  = q^{-L-3/2} [ (L+1/2)/(1-r) + r/(1-r)^2 ]
    part1 = q ** (-L - mp.mpf("1.5")) * (
        (L + mp.mpf("0.5")) / (1 - r) + r / (1 - r) ** 2
    )

    # sum_{n=L}∞ q^{-n-1/2} = q^{-L-1/2}/(1-r)
    part2 = q ** (-L - mp.mpf("0.5")) / (1 - r)

    err_tail = (pi / 2) * aL * (part1 + part2)

    # be conservative: add a tiny extra for roundoff
    err_round = mp.mpf(10) ** (-(mp.mp.dps - 10))
    return val, err_tail + err_round


def sign_interval(val: mp.mpf, err: mp.mpf) -> int:
    """Sign of an interval [val-err, val+err]; returns -1, +1, or 0 if it straddles 0."""
    if val + err < 0:
        return -1
    if val - err > 0:
        return +1
    return 0


# ---------- stable decimal digits in an interval ----------
def stable_decimals(lo: mp.mpf, hi: mp.mpf, maxd: int) -> int:
    """Largest d <= maxd such that floor(lo*10^d) == floor(hi*10^d)."""
    left, right = 0, maxd
    while left < right:
        mid = (left + right + 1) // 2
        factor = mp.power(10, mid)
        if mp.floor(lo * factor) == mp.floor(hi * factor):
            left = mid
        else:
            right = mid - 1
    return left


def format_fixed_from_int(n: int, d: int) -> str:
    s = str(n)
    if d == 0:
        return s
    if len(s) <= d:
        s = "0" * (d - len(s) + 1) + s
    return s[:-d] + "." + s[-d:]


# ---------- main computation ----------
def compute_c_infty(target_digits: int = 200, M: int = 50):
    """
    target_digits: how many digits you *want to aim for* in c. The script will print only
                  the digits it can certify as stable.
    M: split point for the k-sum.
    """
    # Working precision. Increase this if you push target_digits large.
    mp.mp.dps = target_digits + 120

    # Tail length L: rule of thumb + margin
    # L ~ (target_digits / log10(M+2))
    L = int(mp.ceil((target_digits + 40) / mp.log10(M + 2))) + 10

    A = precompute_a(L + 5)

    def Hp(b):
        return Hprime_beta(b, M, L, A)[0]

    # root near 0.7569 (your hint corresponds to beta ~ (1+sin(0.5396))/2 ~ 0.7569)
    beta_star = mp.findroot(Hp, (mp.mpf("0.755"), mp.mpf("0.76")))

    # Build a certified bracket beta_star ± delta by checking signs
    delta = mp.mpf(10) ** (-(target_digits))
    vL, eL = Hprime_beta(beta_star - delta, M, L, A)
    vR, eR = Hprime_beta(beta_star + delta, M, L, A)

    sL = sign_interval(vL, eL)
    sR = sign_interval(vR, eR)

    # If sign is not certified, tighten tail or increase precision / use larger delta
    if not (sL == -1 and sR == +1):
        raise RuntimeError(
            "Could not certify bracket with current settings. "
            "Try increasing mp.mp.dps, increasing L, or using a larger delta."
        )

    beta_lo = beta_star - delta
    beta_hi = beta_star + delta

    # Map to c = asin(2beta - 1)
    c_lo = mp.asin(2 * beta_lo - 1)
    c_hi = mp.asin(2 * beta_hi - 1)

    width = c_hi - c_lo
    half = width / 2
    rel_half = half / abs((c_lo + c_hi) / 2)

    # Determine how many decimals are *stable*
    d_stable = stable_decimals(c_lo, c_hi, maxd=target_digits)

    # Print only the guaranteed digits: use floor(c_lo*10^d_stable)
    factor = mp.power(10, d_stable)
    common_int = int(mp.floor(c_lo * factor))
    c_str = format_fixed_from_int(common_int, d_stable)

    print("Certified stable decimals in c:", d_stable)
    print("c_infty (stable prefix):")
    print(c_str)
    print()
    print(
        "Bracket width  |c_hi - c_lo|  =",
        mp.nstr(width, 8),
        " (~ 10^",
        mp.nstr(mp.log10(width), 6),
        ")",
    )
    print("Half-width     =", mp.nstr(half, 8))
    print("Relative half-width <=", mp.nstr(rel_half, 8))


if __name__ == "__main__":
    compute_c_infty(target_digits=260, M=75)
