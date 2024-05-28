# Mathematical Background
## What is calibration?
Calibration gives us a way to quantify whether the uncertainty in a series of predictions
is faithful to the distribution of prediction errors.
Given a series of predicted probability distributions $p_{\hat{\xi}_i}$ and corresponding
measurements $\xi_i$, we call the predictions *(marginally) calibrated* if the predicted 
probabilities match the observed frequencies of events.
For example, if each prediction is a univariate normal distribution, we expect about 68% of
the observations to fall within one standard deviation of the prediction.

One way to write this mathematically is

$$\mathit{Pr}\left( \xi_{i} \leq \mathit{cdf}(p_{\hat{\xi}_i}, \rho ) \right) \approx \rho,$$

i.e. defining calibration through the cumulative density function
$\mathit{cdf}(p_{(\cdot)}, \rho)$ of the distribution $p_{(\cdot)}$ evaluated at probability
$\rho$.[^cdf] Notably, however, the cumulative density function is not defined for multivariate
distributions.
For this reason we introduce [The Central Prediction Set](@ref) approach implemented in this package.


## The Central Prediction Set
Let us now construct a simple scheme how we can efficiently construct *centered* prediction
sets $\mathcal{P}$ for multivariate normal distributions $D = \mathcal{N}(\mu, \Sigma)$
with $\mu \in \mathbb{R}^d$ and $\Sigma = \Sigma^\top \in \mathbb{R}^{d \times d}$.
These prediction sets will be *centered* in the sense that they are constructed by growing
them symmetrically outwards from the "center" of the distribution.

Our constructed is based on two key observations:
First, any normal distribution can be "diagonalized" without fundamentally changing the
probability density by considering a new distribution $\tilde{D} = \mathcal{N}(\tilde{\mu}, \tilde{\Sigma})$
where $\tilde{\Sigma} = Q\Sigma Q^\top = diag(\tilde{\sigma}_1^2, \dots, \tilde{\sigma}_d^2)$
with a rotation matrix $Q$ s.t. $\|Q\|_2=0$, and appropriately rotated samples $\tilde{x} = Qx.$
Further, the probability set $\mathcal{P}(\tilde{D}, p)$ for diagonalized normal distributions
can simply be constructed by considering each dimension individually, with

$$\mathcal{P}(\tilde{D}, p) = \left\{ \tilde{x} \mid \tilde{x}_i \in \mathcal{P}\left(\mathcal{N}(\tilde{\mu}_i, \tilde{\sigma}_i^2), {p}^{1/d} \right) \forall i \in 1..d \right\}$$

where we can construct $\mathcal{P}$ for the one dimensional case.
We can do this by again picking a centralized construction, although any other constructions, e.g. based on the cumulative density function, are also valid.
Using the central construction, we can normalize the problem by the "length scale" of each univariate normal and compute

$$\tilde{x}_i \in \mathcal{P}\left(\mathcal{N}(\tilde{\mu}_i, \tilde{\sigma}_i^2), p^{1/d} \right) \Leftrightarrow \left| \frac{\tilde{x}_i - \tilde{\mu}_i}{\tilde{\sigma}_i} \right| < quantile\left(\mathcal{N}(0, 1), \frac{{p}^{1/d}}{2} + \frac{1}{2}\right).$$

Finally, we use the above construction to construct the final prediction set 
 $\mathcal{P}(D, p) = \left\{ x \mid Qx \in P(\tilde{D}, p)\right\}$ that is easy and
efficient to check for any $x$.

Further, we can also define the inverse probability set function $\mathcal{P}^{-1}(D, x) = p$,
i.e. the first probability $p$ such that $\mathcal{P}(D, p)$ includes $x$, by
checking the maximum deviation (normalized by length scale) in each dimension, i.e.

$$\mathcal{P}^{-1}(\tilde{D}, x) = \max_i \left(cdf\left(\mathcal{N}(0, 1), \left| \frac{\tilde{x}_i - \tilde{\mu}_i}{\tilde{\sigma}_i} \right|\right) - \frac{1}{2}\right) \cdot 2.$$

## Calibration and Sharpness
Further, we note that a series of predictions may be perfectly calibrated, but still be "bad"
in the sense that each prediction has large uncertainty.
Therefore, we must also measure sharpness, a measure of the "conciseness" of the predictions.
A good prediction is therefore well calibrated and sharp at the same time.
We present an example in [Measuring Sharpness](@ref).
We refer to [gneitingProbabilisticForecasting2014](@citet) for an overview of the trade off
between calibration and sharpness.

```@bibliography
```

[^cdf]: If we only have samples of $\xi_i$, we may replace $\mathit{cdf}$ with the /empirical/ $\mathit{cdf}$.
