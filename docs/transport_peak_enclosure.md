# Scope and derivation of the independent peak verifier

The optimization objective is the maximum heat flux on the closed receiving wall. The verifier operates after the transport solve. It freezes the computed total area A, donor powers w_i, arrival area coordinates c_i, and diffusion widths s_i. It then bounds the maximum of the ideal infinite periodic Gaussian mixture, over the complete receiving coordinate, with arithmetic and series truncation error enclosed.

For normalized receiving coordinate t in [0,1] and m_i=c_i/A:

```
q(t) = sum_i w_i / (sqrt(2*pi)*s_i)
       * sum_{k=-infinity}^{infinity} exp(-0.5 * (A*(t-m_i+k)/s_i)^2).
```

Input binary64 constants are converted exactly to Decimal. When area centers are supplied, c_i/A is divided with outward rounding rather than first rounding the ratio to binary64. Arithmetic uses separate floor/ceiling contexts at 42 significant decimal digits. Correctly rounded exp and sqrt values are enlarged by one adjacent representable value. Pi is enclosed by the existing rational Machin/alternating-series construction. Ambient Decimal precision does not control these operations.

## Omitted periodic images

Both t and m_i lie in [0,1]. If images -K through K are retained, the nearest omitted image is at least K*A away. With r=A_min/s_max, successive omitted exponential terms have ratio at most exp(-(2K+1)*r^2/2). Thus the two omitted tails are bounded above by

```
2 * amplitude_upper * exp(-K^2*r^2/2)
    / (1 - exp(-(2K+1)*r^2/2)).
```

Their lower bound is zero. Far retained images can also use [0, exp(-100)] when their entire exponent interval is below -100; this encloses their contribution rather than deleting or renormalizing it.

For broad kernels a complementary Fourier bound uses a=2*pi^2*(s/A)^2 and

```
sum_{n>=1} exp(-a*n^2) <= exp(-a)/(1-exp(-3*a)).
```

Since |cos|<=1, this bounds the density uniformly in t and avoids poor image-interval conditioning near a uniform distribution. Both bounds enclose the same infinite kernel, so their intersection remains valid.

## Spatial peak enclosure

Every leaf covers a receiving interval. Its density upper enclosure bounds all values in that interval. The largest leaf upper bound therefore bounds the global maximum. Point evaluations provide lower enclosures; the greatest such lower enclosure bounds the maximum from below. Splitting retains both children at the identical Decimal boundary. Upper bounds may be inherited from a parent. No spatial region is removed using sampled heat values.

The stopping criterion is an absolute MW/m² bracket width, checked with upward rounding. A split-budget or coordinate-precision stop returns its actual bracket and does not claim the tolerance was met. Replay reconstructs every leaf and bound from the frozen inputs and rejects modified results.

## Interval-valued transport parameters

The same evaluator accepts explicit ranges of donor powers, arrivals, widths, and period. A lower point enclosure holds for every member of that parameter family, hence it bounds every member's maximum from below. The largest spatial upper enclosure bounds every member's maximum from above. Independent intervals discard correlations, so they can be loose.

An assumed +/-1% width/arrival envelope is a conditional mathematical question. Sampled parameter extrema do not prove that an entire geometry box maps into that envelope. No conditional kernel-family result is promoted to a geometry optimization lower bound without a separate proof of that containment.

## What is and is not certified

- Certified component: complete spatial maximum of the exact infinite kernel defined by the frozen inputs, or a stated interval family, to the recorded bracket/budget.
- Not certified here: exact feasibility of the continuous physical geometry; field interpolation or ODE path error; first-hit branch stability; thermal-solve width error; source quadrature error; production polyline construction error; the production kernel's floating-point/truncation error itself.
- Different tube counts measure numerical sensitivity. Their agreement is not a rigorous error bound on the continuous transport model.
- Missing-hit donors are excluded exactly as in production deposition and remain in the reported power ledger. The verifier concerns the retained mixture; it does not renormalize missing power.
- A fixed-candidate spatial peak enclosure is not a lower bound on the best peak across different geometries. All geometry-global-optimum flags remain false.

Arithmetic reference: [Python Decimal documentation](https://docs.python.org/3.13/library/decimal.html). The documented rounding guarantees are a dependency of this verifier, not proof of the transport model.
