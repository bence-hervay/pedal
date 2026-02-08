# Optimal “always-stand-on-the-forward-pedal” strategy on an ideal bike

Closed-form dynamics, an exact Jacobi-elliptic solution, and an exact transformed optimization to compute the long‑time optimal starting phase $c_\infty$ to high precision.

---

## Contents

- [1. Physical model](#1-physical-model)
- [2. Parameter removal via nondimensionalization](#2-parameter-removal-via-nondimensionalization)
- [3. Dimensionless mathematical problem](#3-dimensionless-mathematical-problem)
- [4. Exact closed-form solution](#4-exact-closed-form-solution)
- [5. Finite-time optimal starting phase](#5-finite-time-optimal-starting-phase)
- [6. Long-time limit phase $c_\infty$ as an exact transformed optimization](#6-long-time-limit-phase-c_infty-as-an-exact-transformed-optimization)
- [7. High-precision evaluation of the limit objective](#7-high-precision-evaluation-of-the-limit-objective)
- [8. Numerical value of $c_\infty$](#8-numerical-value-of-c_infty)
- [9. Notation reference](#9-notation-reference)

---

## 1. Physical model

An idealized bicycle + rider system moves along a straight line. There are two pedals, 180° apart on the crank. The rider always puts full body weight on the pedal that is *forward* (in the direction of travel). The instantaneous driving torque from that pedal is proportional to the cosine of the crank angle.

Let:

- $M$ be the total mass (bike + rider),
- $g$ be gravitational acceleration,
- $r$ be the crank (pedal) radius,
- $\theta$ be the crank angle measured so that “forward” corresponds to $\theta=0$.

Then the (idealized) torque from standing on a pedal is:

```math
\tau(\theta) = (Mg)\,r\,\cos\theta
```

Because the pedals are separated by $\pi$, the other pedal is at $\theta+\pi$ and contributes $\cos(\theta+\pi)=-\cos\theta$. Choosing the **forward** pedal means taking the positive contribution at all times, which yields an effective driving term proportional to:

```math
|\cos\theta|
```

Assume an ideal drivetrain and rolling without slip, so crank angle advances proportionally to traveled distance.

---

## 2. Parameter removal via nondimensionalization

Let $s(t)$ be the physical distance traveled.

Let:

- $R$ be wheel radius,
- $G$ be the (ideal) gear ratio relating wheel angle $\phi$ and crank angle $\theta$ by $\phi=G\theta$.

Rolling without slip gives $\phi = s/R$, hence:

```math
\theta = \frac{s}{GR} + \theta_0
```

where $\theta_0$ is the initial crank phase.

Ideal torque transfer gives wheel torque $\tau_w=\tau/G$, and ground force $F=\tau_w/R$, so:

```math
s''(t) = \frac{F}{M} = \frac{\tau}{MGR}
```

Substitute $\tau(\theta)=(Mg)r\cos\theta$ and the forward-pedal rule:

```math
s''(t) = \frac{gr}{GR}\left|\cos\!\Big(\frac{s(t)}{GR}+\theta_0\Big)\right|
```

Now define dimensionless distance and time:

```math
x = \frac{s}{GR}, \qquad \tau = t\,\sqrt{\frac{gr}{G^2R^2}}
```

Then the dynamics become **parameter-free**:

```math
\frac{d^2x}{d\tau^2} = \left|\cos\!\big(x(\tau)-c\big)\right|
```

with the phase parameter:

```math
c = -\theta_0
```

**Interpretation:** $M,g,r,R,G$ enter only via the trivial rescalings between $(s,t)$ and $(x,\tau)$. The shape of $x(\tau;c)$, the optimal phase, and the universal long‑time optimum $c_\infty$ are independent of $(M,g,r,R,G)$ in the dimensionless formulation.

---

## 3. Dimensionless mathematical problem

For each fixed $c\in[-\pi/2,\pi/2]$, define $x(\tau;c)$ for $\tau\ge 0$ by:

```math
x(0;c)=0,\qquad x'(0;c)=0,\qquad x''(\tau;c)=\bigl|\cos(x(\tau;c)-c)\bigr|
```

Goal (finite time): for each $\tau\ge 0$,

```math
c^*(\tau)\in\arg\max_{c\in[-\pi/2,\pi/2]} x(\tau;c)
```

Goal (long time): the limiting maximizer

```math
c_\infty = \lim_{\tau\to\infty} c^*(\tau)
```

(typically unique in $(-\pi/2,\pi/2)$).

---

## 4. Exact closed-form solution

### 4.1 Shifted variable

Let:

```math
u(\tau)=x(\tau;c)-c
```

Then:

```math
u''(\tau)=|\cos u(\tau)|,\qquad u(0)=-c,\qquad u'(0)=0
```

and $x(\tau;c)=u(\tau)+c$.

### 4.2 Piecewise Jacobi-elliptic closed form

Define:

```math
\beta=\beta(c)=\frac{1+\sin c}{2}\in(0,1)\quad\text{for }|c|<\pi/2
```

Let $K(m)$ be the complete elliptic integral of the first kind (parameter $m$):

```math
K(m)=\int_0^{\pi/2}\frac{d\theta}{\sqrt{1-m\sin^2\theta}}
```

Let $\text{sn}(z\mid m)$ be the Jacobi elliptic sine.

**Endpoints:** if $c=\pm\pi/2$, then the unique solution is $x(\tau;c)\equiv 0$.

Assume $|c|<\pi/2$.

Define the first segment duration:

```math
T_0 = K(\beta)
```

#### Segment 0 (for $0\le \tau\le T_0$)

```math
u(\tau)=\frac{\pi}{2}-2\arcsin\!\Big(\sqrt{\beta}\;\text{sn}\big(K(\beta)-\tau\mid \beta\big)\Big),
\qquad x(\tau;c)=u(\tau)+c
```

#### Later segments (for $\tau>T_0$)

For each integer $k\ge 1$, define:

```math
m_k = \frac{1}{k+\beta}, \qquad \Delta T_k = \frac{K(m_k)}{\sqrt{k+\beta}}
```

and cumulative times:

```math
T_k = T_0 + \sum_{j=1}^k \Delta T_j
```

Given $\tau>T_0$, find the unique $k\ge 1$ such that $T_{k-1}\le \tau \le T_k$, set $\rho=\tau-T_{k-1}$, and then:

```math
u(\tau)=k\pi+\frac{\pi}{2}-2\arcsin\!\Big(\text{sn}\big(K(m_k)-\sqrt{k+\beta}\,\rho\mid m_k\big)\Big),
\qquad x(\tau;c)=u(\tau)+c
```

This is an **exact** closed form, piecewise in $\tau$.

---

## 5. Finite-time optimal starting phase

For each fixed $\tau$, define:

```math
c^*(\tau)\in\arg\max_{c\in[-\pi/2,\pi/2]} x(\tau;c)
```

A robust computational approach is to optimize in 1D over $c$ (independently per $\tau$) using a mix of grid scans, local peak detection, and randomized multi-start refinement.

---

## 6. Long-time limit phase $c_\infty$ as an exact transformed optimization

Empirically and analytically, $c^*(\tau)$ converges as $\tau\to\infty$ to a constant $c_\infty$.

Define:

```math
\beta(c)=\frac{1+\sin c}{2}
```

Define the exact convergent series objective $F(c)$ for $c\in(-\pi/2,\pi/2)$:

```math
F(c)
=
K(\beta(c))
+
\sum_{k=1}^{\infty}
\left[
\frac{K\!\left(\frac{1}{k+\beta(c)}\right)}{\sqrt{k+\beta(c)}}
-
\frac{\pi}{2\sqrt{k}}
\right]
```

Then the long‑time optimal phase is the solution of the **exact transformed optimization**:

```math
\boxed{
c_\infty \in \arg\min_{c\in[-\pi/2,\pi/2]} F(c)
}
```

Equivalently, $c_\infty\in\arg\max (-F(c))$.

---

## 7. High-precision evaluation of the limit objective

A fast high‑precision route uses the exact series for $K(m)$ (for $|m|<1$):

```math
K(m)=\frac{\pi}{2}\sum_{n=0}^{\infty} a_n m^n,
\qquad
a_n=\frac{\binom{2n}{n}^2}{16^n}
```

For large $k$, the parameter $m=\frac{1}{k+\beta}$ is small, giving:

```math
\frac{K\!\left(\frac{1}{k+\beta}\right)}{\sqrt{k+\beta}}
=
\frac{\pi}{2}\sum_{n=0}^{\infty} a_n (k+\beta)^{-(n+1/2)}
```

This converts the infinite tail $\sum_{k>J}$ into Hurwitz zeta values via:

```math
\sum_{k=J+1}^\infty (k+\beta)^{-s} = \zeta(s, J+1+\beta)
```

So $F(c)$ can be evaluated to many digits by:
- summing the first $k=1,\dots,J$ terms directly,
- replacing the remainder by a truncated $K$-series + Hurwitz zeta tail.

Then minimize $F(c)$ (or solve $F'(c)=0$) using high‑precision arithmetic.

---

## 8. Numerical value of $c_\infty$

High‑precision minimization of $F(c)$ yields:

```math
c_\infty \approx
0.53960434045973786251292463554433549832791740569344\ \text{radians}
```

and:

```math
c_\infty \approx 30.91705131528334632675^\circ
```

---

## 9. Notation reference

- $s(t)$: physical distance traveled
- $x(\tau;c)=s/(GR)$: dimensionless distance
- $\tau=t\sqrt{gr/(G^2R^2)}$: dimensionless time
- $c=-\theta_0$: phase parameter (negative initial crank phase)
- $u(\tau)=x(\tau;c)-c$
- $\beta(c)=(1+\sin c)/2$
- $K(m)$: complete elliptic integral of the first kind
- $\text{sn}(z\mid m)$: Jacobi elliptic sine
- $\zeta(s,a)$: Hurwitz zeta function
