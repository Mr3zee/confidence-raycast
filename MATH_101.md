# Math 101 — What This Plugin Actually Computes

A short reference for every formula in `src/statistics.ts`. Section headers map to function names.

---

## 0. Concepts in one minute

- **H₀ (null hypothesis):** "no real difference between A and B."
- **α (significance level):** the chance we reject H₀ when it's actually true (false positive). Typically 0.05.
- **β:** the chance we *fail* to reject H₀ when there really is an effect (false negative).
- **Power = 1 − β:** chance we detect a real effect. Typically 0.80.
- **p-value:** under H₀, the probability of seeing data at least as extreme as what we observed. Reject H₀ when `p < α`.
- **Confidence interval (CI):** a range built so that — over many repeated experiments — the true effect lies inside it `(1−α)·100%` of the time.
- **MDE (minimum detectable effect):** the smallest effect we want the experiment to be powered to find.

---

## 1. Standard Normal Distribution

The standard normal `Z ~ N(0, 1)` underlies every z-based formula here.

### 1.1 CDF — `normalCdf(x)`

$$
\Phi(x) = \tfrac{1}{2}\!\left[1 + \mathrm{erf}\!\left(\tfrac{x}{\sqrt{2}}\right)\right]
$$

We compute `erf` with the Abramowitz–Stegun 7.1.26 approximation (max abs error ≈ 1.5·10⁻⁷):

$$
\mathrm{erf}(x) \approx 1 - (a_1 t + a_2 t^2 + a_3 t^3 + a_4 t^4 + a_5 t^5)\,e^{-x^2},
\quad t = \tfrac{1}{1 + p\,x}
$$

with the published coefficients `a₁..a₅, p`. We mirror it for negative `x` via `erf(-x) = -erf(x)`.

### 1.2 Inverse CDF — `normalInv(p)`

Given `p ∈ (0,1)`, find `x` with `Φ(x) = p`. We use **Acklam's algorithm**: rational approximations on three regions of `p` (`p < 0.02425`, the central body, and `p > 1 − 0.02425`). Accuracy ≈ 10⁻⁹ across the unit interval — plenty for α/2 quantiles.

This gives us the critical values `z_{α/2}` (CIs) and `z_β` (power).

---

## 2. Student's t-Distribution

Used by Welch's t-test for continuous metrics with unknown, possibly unequal variances.

### 2.1 CDF — `tCdf(t, df)`

The t CDF can be written via the **regularized incomplete beta function** `I_x(a, b)`:

$$
P(T \le t \mid \mathrm{df}=\nu) = 1 - \tfrac{1}{2}\, I_{x}\!\left(\tfrac{\nu}{2},\, \tfrac{1}{2}\right) \quad\text{for } t \ge 0,
\quad x = \tfrac{\nu}{\nu + t^2}
$$

(symmetric for `t < 0`). We compute `I_x(a, b)` via:

$$
I_x(a, b) = \frac{x^a (1-x)^b}{a\,B(a,b)} \cdot \mathrm{cf}(a, b, x)
$$

where `B(a,b)` comes from `lgamma` (Lanczos series), and `cf` is Lentz's continued-fraction expansion (`betacf` in code). We use the standard symmetry `I_x(a,b) = 1 − I_{1−x}(b,a)` near the tails for numerical stability.

### 2.2 Inverse — `tInv(p, df)`

Bisection on `tCdf` over `[−10⁶, 10⁶]`, stopping at width `< 10⁻⁹`. Slow but bullet-proof and called only once per CI.

---

## 3. Two-Proportion z-Test — `twoProportionTest`

Used when each unit is a Bernoulli outcome: converted/not, churned/not.

Inputs: `xA, nA, xB, nB`, confidence level `(1−α)·100%`.

### 3.1 Point estimates

$$
\hat p_A = \tfrac{x_A}{n_A},\quad \hat p_B = \tfrac{x_B}{n_B},\quad \widehat{\Delta} = \hat p_B - \hat p_A
$$

### 3.2 P-value: pooled SE under H₀

Under H₀ both groups share one rate, so we pool:

$$
\hat p = \tfrac{x_A + x_B}{n_A + n_B},\quad
\mathrm{SE}_0 = \sqrt{\hat p\,(1-\hat p)\!\left(\tfrac{1}{n_A} + \tfrac{1}{n_B}\right)}
$$

$$
z = \tfrac{\widehat{\Delta}}{\mathrm{SE}_0},\qquad
p\text{-value} = 2\bigl(1 - \Phi(|z|)\bigr)
$$

### 3.3 Confidence interval: unpooled SE

For the CI we *don't* assume H₀ — we estimate each variance separately:

$$
\mathrm{SE} = \sqrt{\tfrac{\hat p_A(1-\hat p_A)}{n_A} + \tfrac{\hat p_B(1-\hat p_B)}{n_B}}
$$

$$
\mathrm{MoE} = z_{1-\alpha/2}\cdot \mathrm{SE},\qquad
\mathrm{CI} = \widehat{\Delta} \pm \mathrm{MoE}
$$

> **Why two SEs?** The pooled form is the right model under H₀ (the question the p-value answers). The unpooled form is the right estimate of the *actual* difference (the question the CI answers). Standard practice; matches Evan Miller, Optimizely, etc.

---

## 4. Welch's Two-Sample t-Test — `welchTTest`

Used when the metric is continuous (revenue per user, days-to-churn, session length) and variances may differ.

Inputs: means `m_A, m_B`, std devs `s_A, s_B`, sample sizes `n_A, n_B`.

### 4.1 Test statistic

$$
\mathrm{SE} = \sqrt{\tfrac{s_A^2}{n_A} + \tfrac{s_B^2}{n_B}},\qquad
t = \tfrac{m_B - m_A}{\mathrm{SE}}
$$

### 4.2 Welch–Satterthwaite degrees of freedom

$$
\nu = \frac{\bigl(s_A^2/n_A + s_B^2/n_B\bigr)^{2}}
           {\dfrac{(s_A^2/n_A)^2}{n_A - 1} + \dfrac{(s_B^2/n_B)^2}{n_B - 1}}
$$

(generally non-integer, that's expected).

### 4.3 P-value & CI

$$
p\text{-value} = 2\bigl(1 - T_\nu(|t|)\bigr),\qquad
\mathrm{CI} = (m_B - m_A) \pm t_{1-\alpha/2,\,\nu}\cdot \mathrm{SE}
$$

---

## 5. Sample Size for a Two-Proportion Test — `sampleSizeProportion`

Inputs: baseline `p₁`, MDE (relative or absolute) → `p₂`, α, power `1−β`, one/two-sided.

### 5.1 MDE conversion

$$
p_2 = \begin{cases}
p_1\,(1 + \mathrm{MDE}) & \text{if relative} \\
p_1 + \mathrm{MDE} & \text{if absolute}
\end{cases}
$$

with the constraint `0 < p₂ < 1`.

### 5.2 Per-group sample size

We need the test (under H₀, pooled variance) to reject with probability `1 − β` when the truth is `(p₁, p₂)`:

$$
n = \frac{\Bigl(z_{1-\alpha/2}\sqrt{2\bar p\,\bar q} + z_{1-\beta}\sqrt{p_1 q_1 + p_2 q_2}\Bigr)^{2}}{(p_2 - p_1)^2}
$$

where `q = 1 − p`, `p̄ = (p₁ + p₂)/2`, `q̄ = 1 − p̄`. We `Math.ceil` the result and report `total = 2n`. For one-sided tests, replace `z_{1-α/2}` with `z_{1-α}`.

> **Intuition.** Numerator: how strict the gate is (`z_α`) plus how confidently you want to clear it (`z_β`), each weighted by the relevant variance. Denominator: the squared signal you're trying to find. Bigger signal → smaller `n`. Tighter α or higher power → bigger `n`.

---

## 6. Sample Size for a Two-Sample Means Test — `sampleSizeMeans`

Inputs: `s_A, s_B` (treatment defaults to control), absolute MDE `Δ`, α, power.

$$
n = \frac{\bigl(z_{1-\alpha/2} + z_{1-\beta}\bigr)^{2}\,(s_A^2 + s_B^2)}{\Delta^{2}}
$$

This is the **z-approximation** to the t-based formula — fine for `n ≳ 30` per group. For very small expected `n` you'd iterate using `t_{1-α/2,\,ν(n)}`, but planning calculators almost universally stop at the z form.

---

## 7. Worked Example — Churn (used in this README's docs)

- Control: 96 / 1,200 → `p̂_A = 0.080`
- Treatment: 75 / 1,250 → `p̂_B = 0.060`
- Δ̂ = −0.020 absolute, ≈ −25% relative

Pooled: `p̂ = 171/2450 ≈ 0.0698`, `SE₀ ≈ 0.01037`, `z ≈ −1.929`, `p ≈ 0.0537`.

Unpooled: `SE ≈ 0.01036`, `MoE = 1.96·SE ≈ 0.0203`, `95% CI ≈ [−0.0403, +0.0003]`.

So at α = 0.05 the result is right on the line (p just above 0.05; CI just contains 0). You'd want either a larger sample or to accept slightly lower confidence.

Sample size to detect a more conservative 10%-relative MDE (8% → 7.2%, two-sided, α=0.05, power=0.80):

$$
n \approx \frac{(1.96\sqrt{2\cdot 0.076\cdot 0.924} + 0.842\sqrt{0.0736 + 0.0668})^2}{(0.008)^2} \approx 17{,}800
$$

per group — i.e., the experiment above (~1,200 + 1,250) was simply underpowered.

---

## 8. Assumptions & Failure Modes

- **Independence.** All formulas assume one observation = one independent unit. Multi-seat accounts, repeat sessions per user, or network effects break this and inflate apparent significance. Aggregate to the user/account level first.
- **Normal approximation for proportions.** The z-test is reliable when `n·p ≥ 10` *and* `n·(1−p) ≥ 10` for both groups. For very rare events use Fisher's exact instead.
- **Equal allocation.** Sample-size formulas above assume 50/50 splits. Unequal splits inflate total `n`; multiply by `(1 + k)² / (4k)` where `k = n_B/n_A`.
- **Two-sided by default.** One-sided tests double power *only if* you'd genuinely never act on a result in the other direction — rarely true in practice.
- **Peeking.** These formulas assume a single, pre-planned analysis. Repeated checking inflates the false-positive rate; use sequential or Bayesian methods if you must look early.
