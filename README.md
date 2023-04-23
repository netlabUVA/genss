# genss: Simulate and Fit General State-Space Models using Multiple Iterated Filtering

This package presents a framework for specifying and estimating state-space models, with particular application to psychological and behavioral timeseries. This package is in active development and currently has limited features:

- Estimate the parameters for a discrete time two state model (linear state dynamics) with normally distributed measurements, of the following form:
$$\mathbf{x}_{t+1} = \mathbf{A}\mathbf{x}_t + \boldsymbol{\varepsilon}_t$$
$$\mathbf{y}_t = \mathbf{H}\mathbf{x}_t + \boldsymbol{\varsigma}_t$$
$$\boldsymbol{\varepsilon_t} \sim N(0, \boldsymbol{\Sigma})$$
$$\boldsymbol{\varsigma_t} \sim N(0,\boldsymbol{\Theta})$$

- Estimate the parameters for a discrete time two state model (linear state dynamics) with ordinally distributed measurements (graded response model) of the following form:

$$\mathbf{x}_{t+1} = \mathbf{Ax}_t + \boldsymbol{\varepsilon}_t$$

$$\boldsymbol{\varepsilon}_t \sim N_p(0, \boldsymbol{\Sigma})$$

$$p(y_{it} > j | x_t) = \frac{1}{1+\exp[-(\boldsymbol{H}\mathbf{x}_t - \beta_{ij})]}$$

`genss` currently uses `pomp` to construct and fit these models using Multiple Iterated Filtering. To identify both types of models, the marginal variance of the state process is constrained by solving the following equation for the diagonal elements of $\Sigma$ ($\Gamma$ is the marginal covariance matrix, with diagonal elements constrained to 1, while $\Sigma$ has off diagonal elements constrained to 0):
$$vec(\boldsymbol{\Gamma})=(\mathbf{I}-\mathbf{A}\otimes\mathbf{A} )^{-1} vec(\boldsymbol{\Sigma})$$

## Installation

## Interested in Using `genss`?

Get in touch with Teague Henry (trhenry@virginia.edu) and let us know your use case!
