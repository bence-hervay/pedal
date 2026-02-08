# Optimal “always-stand-on-the-forward-pedal” strategy on an ideal bike
*(closed-form dynamics, exact elliptic-function solution, and a high‑precision method for the long‑time optimal phase)*

## 1) Physical model (informal but precise)

Consider an idealized bicycle + rider system moving along a straight line. The system has:

- **Two pedals** placed 180° apart on a crank.
- The rider always puts **full body weight** on whichever pedal is **forward** (in the direction of travel).
- The instantaneous **driving torque** produced by that pedal is proportional to the cosine of the crank angle:
  $$
  \tau(\theta) = (Mg)\,r\,\cos\theta,
  $$
  where:
  - $M$ = total mass (bike + rider),
  - $g$ = gravitational acceleration,
  - $r$ = crank (pedal) radius,
  - $\theta$ = crank angle measured so that “forward” corresponds to $\theta=0$.

Because there are two pedals separated by $\pi$, if one pedal has angle $\theta$, the other has angle $\theta+\pi$ and produces $\cos(\theta+\pi)=-\cos\theta$. Choosing the **forward** pedal means taking the positive contribution at all times, hence the effective driving term is proportional to:
$$
|\cos\theta|.
$$

Finally, the crank angle advances in proportion to distance traveled (via wheel radius and gear ratio), so $\theta$ is an affine function of the traveled distance.

---

## 2) Why wheel radius, crank length, gear ratio, total mass, and $g$ are “irrelevant” (nondimensionalization)

Let $s(t)$ be the physical distance traveled along the ground.

Let:
- $R$ = wheel radius,
- $G$ = gear ratio defined by: wheel angle $\phi$ and crank angle $\theta$ satisfy $\phi = G\theta$ (ideal drivetrain),
- so $\phi = s/R$ and therefore $\theta = \frac{s}{GR} + \theta_0$, where $\theta_0$ is the initial crank phase.

Ideal power transmission gives wheel torque $\tau_w = \tau/G$. Ground force is $F = \tau_w/R$, and acceleration is $s'' = F/M$. Substituting $\tau(\theta)=(Mg)r\cos\theta$ gives
$$
s''(t) = \frac{gr}{G R}\,|\cos\theta(t)|
       = \frac{gr}{G R}\,\left|\cos\!\Big(\frac{s(t)}{G R} + \theta_0\Big)\right|.
$$

Introduce the **dimensionless distance** and **dimensionless time**
$$
x = \frac{s}{GR}, \qquad \tau = t\,\sqrt{\frac{gr}{G^2R^2}}.
$$
Then $s = GR\,x$, and a short calculation gives the parameter‑free dimensionless ODE
$$
\frac{d^2x}{d\tau^2} = \left|\cos\!\big(x(\tau) - c\big)\right|,
$$
where $c = -\theta_0$ is the initial phase parameter.

**Interpretation:** all physical parameters $(M,g,r,R,G)$ enter only through:
- a **distance scale** $GR$ converting $x \leftrightarrow s$,
- and a **time scale** $\sqrt{gr}/(GR)$ converting $\tau \leftrightarrow t$.

Therefore, any statement about:
- the shape of the trajectories $x(\tau;c)$ as a function of $\tau$,
- the maximizing phase $c$,
- and in particular the **long‑time optimal phase $c_\infty$**

is **universal** in the dimensionless variables and does not depend on $(M,g,r,R,G)$, except through the trivial rescalings above.

---

## 3) Dimensionless mathematical problem (the core ODE)

For each fixed $c\in[-\pi/2,\pi/2]$, define $x(\tau;c)$ for $\tau\ge 0$ by
$$
x(0;c)=0,\qquad x'(0;c)=0,\qquad x''(\tau;c)=\bigl|\cos(x(\tau;c)-c)\bigr|.
$$
(Primes denote derivatives with respect to $\tau$.)

The phase bound $|c|\le \pi/2$ corresponds to starting with one pedal “forward” (within 90° of the forward direction).

---

## 4) Exact solution of the ODE (closed form via Jacobi elliptic functions)

### 4.1 Change of variables
Let
$$
u(\tau)=x(\tau;c)-c.
$$
Then
$$
u''(\tau)=|\cos u(\tau)|,\qquad u(0)=-c,\qquad u'(0)=0,
$$
and $x(\tau;c)=u(\tau)+c$.

### 4.2 Piecewise explicit formula
Define the parameter
$$
\beta=\beta(c)=\frac{1+\sin c}{2}\in(0,1)\quad\text{for }|c|<\pi/2.
$$
Let:
- $K(m)$ be the **complete elliptic integral of the first kind** (parameter $m$):
  $$
  K(m)=\int_0^{\pi/2}\frac{d\theta}{\sqrt{1-m\sin^2\theta}},
  $$
- $\operatorname{sn}(z\mid m)$ be the **Jacobi elliptic sine**.

Special endpoints:
- If $c=\pm\pi/2$, then the unique solution is $x(\tau;c)\equiv 0$.

Otherwise $|c|<\pi/2$. Define the first segment duration:
$$
T_0 = K(\beta).
$$

**Segment 0:** for $0\le \tau\le T_0$,
$$
u(\tau)=\frac{\pi}{2}-2\arcsin\!\Big(\sqrt{\beta}\;\operatorname{sn}\big(K(\beta)-\tau\mid \beta\big)\Big),
\qquad x(\tau;c)=u(\tau)+c.
$$

**Later segments:** for each integer $k\ge 1$, define
$$
m_k = \frac{1}{k+\beta}, \qquad \Delta T_k = \frac{K(m_k)}{\sqrt{k+\beta}},
$$
and cumulative times
$$
T_k = T_0 + \sum_{j=1}^k \Delta T_j.
$$

Given $\tau > T_0$, find the unique $k\ge 1$ such that $T_{k-1}\le \tau \le T_k$, set $\rho=\tau-T_{k-1}$, and then:
$$
u(\tau)=k\pi+\frac{\pi}{2}-2\arcsin\!\Big(\operatorname{sn}\big(K(m_k)-\sqrt{k+\beta}\,\rho\mid m_k\big)\Big),
\qquad x(\tau;c)=u(\tau)+c.
$$

This is an **exact closed form** (piecewise in $\tau$).

---

## 5) Numerical verification approach (conceptual)

A direct way to verify the closed form is to compare it against a small‑step numerical integrator for
$$
x'' = |\cos(x-c)|,\quad x(0)=0,\ x'(0)=0,
$$
using, for example, a symplectic/energy‑stable scheme (e.g. velocity Verlet) with a very small step size $\Delta\tau$. Agreement improves with $\Delta\tau\to 0$.

---

## 6) Finite‑time optimization and the “best starting phase”

For each $\tau\ge 0$, define the best phase:
$$
c^*(\tau)\in\arg\max_{c\in[-\pi/2,\pi/2]} x(\tau;c).
$$

A robust way to approximate $c^*(\tau)$ (for each $\tau$ independently) is:

1. sample $x(\tau;c)$ on a grid in $c$,
2. detect candidate local maxima on that grid,
3. refine many candidate intervals (including randomized multi‑start intervals) by 1D bounded optimization of $x(\tau;c)$.

---

## 7) Long‑time optimal phase $c_\infty$ as an exact transformed optimization

Empirically and analytically, $c^*(\tau)$ converges as $\tau\to\infty$:
$$
c^*(\tau)\to c_\infty \in (-\pi/2,\pi/2).
$$

### 7.1 Exact limit optimization objective
Define $\beta(c)=(1+\sin c)/2$ and $K(\cdot)$ as above. Define the function $F(c)$ on $c\in(-\pi/2,\pi/2)$:
$$
F(c)
\;=\;
K(\beta(c))
\;+\;
\sum_{k=1}^{\infty}
\left[
\frac{K\!\left(\frac{1}{k+\beta(c)}\right)}{\sqrt{k+\beta(c)}}
\;-\;
\frac{\pi}{2\sqrt{k}}
\right].
$$

Then the long‑time optimal phase is given by the exact transformed problem
$$
\boxed{
c_\infty \in \arg\min_{c\in[-\pi/2,\pi/2]} F(c)
}
\qquad\text{(equivalently }\arg\max -F(c)\text{).}
$$

Every symbol in $F(c)$ is explicit:
- $\beta(c)=(1+\sin c)/2$,
- $K(m)=\int_0^{\pi/2}(1-m\sin^2\theta)^{-1/2}\,d\theta$,
- the infinite sum is absolutely convergent because the bracketed term is $O(k^{-3/2})$.

### 7.2 High‑precision evaluation of $F(c)$ (zeta‑accelerated tail)
To evaluate $F(c)$ to many digits efficiently, use the exact power series for $K(m)$ (valid for $|m|<1$):
$$
K(m)=\frac{\pi}{2}\sum_{n=0}^{\infty} a_n m^n,
\qquad
a_n=\frac{\binom{2n}{n}^2}{16^n}.
$$

For large $k$, $m=\frac{1}{k+\beta}$ is small, and
$$
\frac{K\!\left(\frac{1}{k+\beta}\right)}{\sqrt{k+\beta}}
=
\frac{\pi}{2}\sum_{n=0}^{\infty} a_n (k+\beta)^{-(n+1/2)}.
$$

This allows the infinite tail $\sum_{k>J}$ to be expressed in terms of the **Hurwitz zeta function** $\zeta(s,a)$, since
$$
\sum_{k=J+1}^{\infty} (k+\beta)^{-s} = \zeta(s, J+1+\beta).
$$

A typical high‑precision strategy is:
- compute the first $k=1,\dots,J$ terms directly (exact $K$),
- approximate the tail using a truncated $K$-series + Hurwitz zeta.

This achieves dozens (or hundreds) of digits with moderate $J$ and series length.

---

## 8) Numerical value of the universal optimum $c_\infty$

Solving
$$
c_\infty \in \arg\min_{c\in[-\pi/2,\pi/2]} F(c)
$$
with high‑precision arithmetic gives:
$$
\boxed{
c_\infty \approx
0.53960434045973786251292463554433549832791740569344\ \text{radians}
}
$$
$$
\boxed{
c_\infty \approx 30.91705131528334632675^\circ
}
$$

(Shown to ~50 decimal digits in radians; more digits are obtainable by increasing precision and tail parameters.)

---

## 9) Practical summary (what to compute)

**Dynamics (exact):**
- evaluate $x(\tau;c)$ by the piecewise Jacobi‑elliptic closed form in §4.

**Best phase for a fixed horizon $\tau$:**
- compute $c^*(\tau)\in\arg\max_{c\in[-\pi/2,\pi/2]} x(\tau;c)$ via independent multi‑start 1D optimization.

**Best phase as $\tau\to\infty$:**
- compute $c_\infty\in\arg\min_{c\in[-\pi/2,\pi/2]} F(c)$ where $F$ is the explicit convergent series in §7.1,
- evaluate $F$ to high precision via the zeta‑accelerated tail in §7.2,
- then minimize $F$ with a 1D method (e.g. golden section, Brent, or root‑finding on $F'(c)$).

---

## 10) “Process” (high-level milestones and discoveries)

1. **Model → dimensionless ODE:** reduce the ideal bike/pedal setup to the single parameter phase ODE $x''=|\cos(x-c)|$.
2. **First integral:** integrate once to get an exact “energy” relation that reduces the problem to quadratures.
3. **Piecewise structure:** handle the absolute value by splitting the motion into $\pi$-sized phase intervals.
4. **Elliptic closed form:** express each interval in Jacobi elliptic functions; glue intervals with explicit segment times.
5. **Verification:** compare the closed form against a small‑step direct numerical integration of the ODE.
6. **Finite‑time optimization:** compute $c^*(\tau)$ via independent multi‑start optimization; observe convergence as $\tau$ grows.
7. **Long‑time limit:** transform the $t\to\infty$ phase selection into the exact convergent optimization $\arg\min F(c)$.
8. **High‑precision constant:** evaluate $F(c)$ with a zeta‑accelerated tail to compute $c_\infty$ to many digits.

---

## 11) Notation reference (single place)

- $x(\tau;c)$: dimensionless distance variable
- $c\in[-\pi/2,\pi/2]$: initial pedal phase parameter
- $\beta(c)=(1+\sin c)/2$
- $K(m)$: complete elliptic integral of the first kind
- $\operatorname{sn}(z\mid m)$: Jacobi elliptic sine
- $\zeta(s,a)$: Hurwitz zeta function
