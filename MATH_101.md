# Math 101 — What This Plugin Actually Computes

A short reference for every formula in `src/statistics.ts`. Section headers map to function names.

---

## 0. Concepts in one minute

* **`H₀` (null hypothesis):** "No real difference between A and B."
* **`α` (significance level):** The chance we reject `H₀` when it is actually true (false positive). Typically 0.05.
* **`β`:** The chance we *fail* to reject `H₀` when there really is an effect (false negative).
* **Power (`1 − β`):** The chance we detect a real effect. Typically 0.80.
* **`p-value`:** Under `H₀`, the probability of seeing data at least as extreme as what we observed. Reject `H₀` when `p < α`.
* **Confidence interval (CI):** A range built so that—over many repeated experiments—the true effect lies inside it `(1−α)·100%` of the time.
* **MDE (minimum detectable effect):** The smallest effect we want the experiment to be powered to find.

---

## 1. Standard Normal Distribution

The standard normal `Z ~ N(0, 1)` underlies every z-based formula here.

### 1.1 CDF — `normalCdf(x)`

$$
\Phi(x) = \tfrac{1}{2}\!\left[1 + \mathrm{erf}\!\left(\tfrac{x}{\sqrt{2}}\right)\right]
$$

We compute `erf` with the Abramowitz–Stegun 7.1.26 approximation (max absolute error ≈ 1.5·10⁻⁷):

$$
\mathrm{erf}(x) \approx 1 - (a_1 t + a_2 t^2 + a_3 t^3 + a_4 t^4 + a_5 t^5)\,e^{-x^2}, \quad t = \tfrac{1}{1 + p\,x}
$$

This uses the published coefficients `a₁..a₅, p`. We mirror it for negative `x` via `erf(-x) = -erf(x)`.

### 1.2 Inverse CDF — `normalInv(p)`

Given `p ∈ (0,1)`, find `x` with `Φ(x) = p`. We use **Acklam's algorithm**: rational approximations on three regions of `p` (`p < 0.02425`, the central body, and `p > 1 − 0.02425`). Accuracy is ≈ 10⁻⁹ across the unit interval, which is plenty for `α/2` quantiles.

This gives us the critical values `z_{α/2}` (for CIs) and `z_β` (for power).

---

## 2. Student's t-Distribution

This is used by Welch's t-test for continuous metrics with unknown, possibly unequal variances.

### 2.1 CDF — `tCdf(t, df)`

The t-distribution CDF can be written via the **regularized incomplete beta function** `I_x(a, b)`:

$$
P(T \le t \mid \mathrm{df}=\nu) = 1 - \tfrac{1}{2}\, I_{x}\!\left(\tfrac{\nu}{2},\, \tfrac{1}{2}\right) \quad\text{for } t \ge 0, \quad x = \tfrac{\nu}{\nu + t^2}
$$

This is symmetric for `t < 0`. We compute `I_x(a, b)` via:

$$
I_x(a, b) = \frac{x^a (1-x)^b}{a\,B(a,b)} \cdot \mathrm{cf}(a, b, x)
$$

Here, `B(a,b)` comes from `lgamma` (Lanczos series), and `cf` is Lentz's continued-fraction expansion (`betacf` in code). We use the standard symmetry `I_x(a,b) = 1 − I_{1−x}(b,a)` near the tails for numerical stability.

### 2.2 Inverse — `tInv(p, df)`

We use bisection on `tCdf` over `[−10⁶, 10⁶]`, stopping at a width `< 10⁻⁹`. It is slow but bulletproof and is called only once per CI.

---

## 3. One-Proportion z-Test — `oneProportionTest`

Used when you have a single Bernoulli sample and want to compare it against a fixed, known rate `p₀`—such as an SLA you are supposed to hit, last quarter's conversion rate, or textbook prevalence.

Inputs: successes `x`, trials `n`, null proportion `p₀`, and confidence level `(1−α)·100%`.

### 3.1 Point estimate

$$
\hat p = \tfrac{x}{n},\qquad \widehat{\Delta} = \hat p - p_0
$$

### 3.2 P-value: SE under H₀

Under `H₀` the rate is *exactly* `p₀`, so the variance is fixed (no need to estimate it from the sample):

$$
\mathrm{SE}_0 = \sqrt{\tfrac{p_0\,(1-p_0)}{n}},\qquad z = \tfrac{\hat p - p_0}{\mathrm{SE}_0},\qquad p\text{-value} = 2\bigl(1 - \Phi(|z|)\bigr)
$$

### 3.3 Confidence interval (Wald)

The CI here estimates the *true rate `p`*, not the deviation. We use the data—not the null—to estimate variance:

$$
\mathrm{SE} = \sqrt{\tfrac{\hat p\,(1-\hat p)}{n}},\qquad \mathrm{CI} = \hat p \pm z_{1-\alpha/2}\cdot \mathrm{SE}
$$

> **Wald vs. Wilson:** The Wald interval can perform poorly for very small `n` or when `p̂` is near 0/1; in those regimes, the **Wilson** or **Agresti–Coull** intervals are more accurate. We use Wald to match the parallel construction in the two-proportion test; for `n·p̂ ≥ 10` and `n·(1−p̂) ≥ 10`, it works fine.

> **Why two SEs?** Same reason as the two-proportion case: the p-value answers a question about `H₀` (so use `p₀` in the SE), while the CI answers a question about the truth (so use `p̂`).

---

## 4. Two-Proportion z-Test — `twoProportionTest`

Used when each unit is a Bernoulli outcome: converted/not, churned/not.

Inputs: `xA, nA, xB, nB`, and confidence level `(1−α)·100%`.

### 4.1 Point estimates

$$
\hat p_A = \tfrac{x_A}{n_A},\quad \hat p_B = \tfrac{x_B}{n_B},\quad \widehat{\Delta} = \hat p_B - \hat p_A
$$

### 4.2 P-value: pooled SE under H₀

Under `H₀` both groups share one rate, so we pool them:

$$
\hat p = \tfrac{x_A + x_B}{n_A + n_B},\quad \mathrm{SE}_0 = \sqrt{\hat p\,(1-\hat p)\!\left(\tfrac{1}{n_A} + \tfrac{1}{n_B}\right)}
$$

$$
z = \tfrac{\widehat{\Delta}}{\mathrm{SE}_0},\qquad p\text{-value} = 2\bigl(1 - \Phi(|z|)\bigr)
$$

### 4.3 Confidence interval: unpooled SE

For the CI we *do not* assume `H₀`—we estimate each variance separately:

$$
\mathrm{SE} = \sqrt{\tfrac{\hat p_A(1-\hat p_A)}{n_A} + \tfrac{\hat p_B(1-\hat p_B)}{n_B}}
$$

$$
\mathrm{MoE} = z_{1-\alpha/2}\cdot \mathrm{SE},\qquad \mathrm{CI} = \widehat{\Delta} \pm \mathrm{MoE}
$$

> **Why two SEs?** The pooled form is the right model under `H₀` (the question the p-value answers). The unpooled form is the right estimate of the *actual* difference (the question the CI answers). This is standard practice and matches tools like Evan Miller and Optimizely.

---

## 5. Welch's Two-Sample t-Test — `welchTTest`

Used when the metric is continuous (e.g., revenue per user, days-to-churn, session length) and variances may differ.

Inputs: means `m_A, m_B`, standard deviations `s_A, s_B`, and sample sizes `n_A, n_B`.

### 5.1 Test statistic

$$
\mathrm{SE} = \sqrt{\tfrac{s_A^2}{n_A} + \tfrac{s_B^2}{n_B}},\qquad t = \tfrac{m_B - m_A}{\mathrm{SE}}
$$

### 5.2 Welch–Satterthwaite degrees of freedom

$$
\nu = \frac{\bigl(s_A^2/n_A + s_B^2/n_B\bigr)^{2}}{\dfrac{(s_A^2/n_A)^2}{n_A - 1} + \dfrac{(s_B^2/n_B)^2}{n_B - 1}}
$$

(This is generally a non-integer, which is expected).

### 5.3 P-value & CI

$$
p\text{-value} = 2\bigl(1 - T_\nu(|t|)\bigr),\qquad \mathrm{CI} = (m_B - m_A) \pm t_{1-\alpha/2,\,\nu}\cdot \mathrm{SE}
$$

---

## 6. Sample Size for a One-Proportion Test — `sampleSizeOneProportion`

Inputs: null proportion `p₀`, MDE (relative or absolute) → alternative `p₁`, `α`, power `1−β`, and one/two-sided.

### 6.1 MDE conversion

$$
p_1 = \begin{cases} p_0\,(1 + \mathrm{MDE}) & \text{if relative} \\ p_0 + \mathrm{MDE} & \text{if absolute} \end{cases}
$$

With the constraint `0 < p₁ < 1`.

### 6.2 Sample size

We want the one-sample z-test (with SE evaluated under `H₀`) to reject with probability `1 − β` when the truth is `p₁`:

$$
n = \frac{\Bigl(z_{1-\alpha/2}\sqrt{p_0\,(1-p_0)} + z_{1-\beta}\sqrt{p_1\,(1-p_1)}\Bigr)^{2}}{(p_1 - p_0)^2}
$$

We use `Math.ceil`. Unlike the two-sample formulas, there is only one group, so `total = n`. For one-sided tests, replace `z_{1-α/2}` with `z_{1-α}`.

---

## 7. Sample Size for a Two-Proportion Test — `sampleSizeProportion`

Inputs: baseline `p₁`, MDE (relative or absolute) → `p₂`, `α`, power `1−β`, and one/two-sided.

### 7.1 MDE conversion

$$
p_2 = \begin{cases} p_1\,(1 + \mathrm{MDE}) & \text{if relative} \\ p_1 + \mathrm{MDE} & \text{if absolute} \end{cases}
$$

With the constraint `0 < p₂ < 1`.

### 7.2 Per-group sample size

We need the test (under `H₀`, pooled variance) to reject with probability `1 − β` when the truth is `(p₁, p₂)`:

$$
n = \frac{\Bigl(z_{1-\alpha/2}\sqrt{2\bar p\,\bar q} + z_{1-\beta}\sqrt{p_1 q_1 + p_2 q_2}\Bigr)^{2}}{(p_2 - p_1)^2}
$$

Where `q = 1 − p`, `p̄ = (p₁ + p₂)/2`, and `q̄ = 1 − p̄`. We `Math.ceil` the result and report `total = 2n`. For one-sided tests, replace `z_{1-α/2}` with `z_{1-α}`.

> **Intuition:** Numerator: how strict the gate is (`z_α`) plus how confidently you want to clear it (`z_β`), each weighted by the relevant variance. Denominator: the squared signal you are trying to find. Bigger signal → smaller `n`. Tighter `α` or higher power → bigger `n`.

---

## 8. Sample Size for a Two-Sample Means Test — `sampleSizeMeans`

Inputs: `s_A, s_B` (treatment defaults to control), absolute MDE `Δ`, `α`, and power.

$$
n = \frac{\bigl(z_{1-\alpha/2} + z_{1-\beta}\bigr)^{2}\,(s_A^2 + s_B^2)}{\Delta^{2}}
$$

This is the **z-approximation** to the t-based formula—fine for `n ≳ 30` per group. For very small expected `n`, you would iterate using `t_{1-α/2,\,ν(n)}`, but planning calculators almost universally stop at the z-form.

---

## 9. Worked Example — Churn (used in this README's docs)

* Control: 96 / 1,200 → `p̂_A = 0.080`
* Treatment: 75 / 1,250 → `p̂_B = 0.060`
* `Δ̂ = −0.020` absolute, ≈ −25% relative

Pooled: `p̂ = 171/2450 ≈ 0.0698`, `SE₀ ≈ 0.01030`, `z ≈ −1.942`, `p ≈ 0.052`.

Unpooled: `SE ≈ 0.01032`, `MoE = 1.96·SE ≈ 0.0202`, 95% CI ≈ `[−0.0402, +0.0002]`.

So at `α = 0.05`, the result is right on the line (`p` is just above 0.05; CI just contains 0). You would want either a larger sample or to accept slightly lower confidence.

Sample size to detect a more conservative 10%-relative MDE (8% → 7.2%, two-sided, `α=0.05`, `power=0.80`):

$$
n \approx \frac{(1.96\sqrt{2\cdot 0.076\cdot 0.924} + 0.842\sqrt{0.0736 + 0.0668})^2}{(0.008)^2} \approx 17{,}227
$$

Per group—meaning the experiment above (~1,200 + 1,250) was simply underpowered.

---

## 10. Assumptions & Failure Modes

* **Independence:** All formulas assume one observation equals one independent unit. Multi-seat accounts, repeat sessions per user, or network effects break this and inflate apparent significance. Aggregate to the user/account level first.
* **Normal approximation for proportions:** The z-test is reliable when `n·p ≥ 10` *and* `n·(1−p) ≥ 10` for both groups. For very rare events, use Fisher's exact test instead.
* **Equal allocation:** Sample-size formulas above assume 50/50 splits. Unequal splits inflate total `n`; multiply by `(1 + k)² / (4k)` where `k = n_B/n_A`.
* **Two-sided by default:** One-sided tests double power *only if* you would genuinely never act on a result in the other direction—which is rarely true in practice.
* **Peeking:** These formulas assume a single, pre-planned analysis. Repeated checking inflates the false-positive rate; use sequential or Bayesian methods if you must look early.
