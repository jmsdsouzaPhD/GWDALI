# **GWDALI (v 1.0)**

A Fisher-Based Software for Parameter Estimation from Gravitational Waves

![Logo](https://github.com/jmsdsouzaPhD/GWDALI/blob/main/docs/source/_static/logo_gwdali.png)

Documentation in [https://gwdali.readthedocs.io/en/latest/](https://gwdali.readthedocs.io/en/latest/)

Main Paper: [arXiv:2307.10154](https://arxiv.org/abs/2307.10154) and [Astronomy and Computing, v. 45, p. 100759, 2023](https://www.sciencedirect.com/science/article/abs/pii/S2213133723000744)

Paper about implementation of automatic differentiation: [arXiv:2510.16955](https://arxiv.org/abs/2510.16955)

References: [arXiv:1401.6892](https://arxiv.org/abs/1401.6892) and [arXiv:2203.02670](https://arxiv.org/abs/2203.02670)

## Installation

To install the software run the command below:

```
pip install gwdali
```

## Requirements

### JAX

The new version of **GWDALI** uses **JAX** to accelerate the computation of waveforms, derivatives, and likelihoods.

We strongly recommend installing **JAX** with *conda* before installing **GWDALI**:

```
conda install -c conda-forge jax
```

Please make sure that **JAX** is installed before running:

```
pip install gwdali
```

Otherwise, pip will attempt to install **JAX** and its dependencies automatically, which may lead to issues with jaxlib on some systems. For more information, see the official [JAX installation guide](https://docs.jax.dev/en/latest/installation.html).

### lalsuite/lalsimulation

To be able to use **LAL** waveforms to compute GW polarizations/strains install the packages [lalsuite, lalsimulation](https://wiki.ligo.org/Computing/LALSuite). It is recommended to use *conda*.

```
conda install lalsuite -c conda-forge
conda install lalsimulation -c conda-forge
```

## Documentation

Available in [https://gwdali.readthedocs.io/en/latest/](https://gwdali.readthedocs.io/en/latest/)
    
## Functionalities

- **get_hphx()**: It returns plus/cross polarizations in the frequency space (SPA);
- **get_strain()**: It retuns detector strains (signals) in the frequency space;
- **get_SNR()**: It retuns detector-network signal-to-noise ratios (individuals and net);
- **draw_detectors()**: It returns a world map showing the chosen detector network configuration;
- **get_derivatives()**: It returns detector signal derivatives;
- **get_tensors()**: It returns DALI tensors including Fisher matrix;
- **Priors()**: Check/Visualize priors to be used in Posterior evaluations;
- **GWDALI()**: Get MCMC/Fisher-Inversion Samples or Posterior-Grid Arrays;

## Authors

- **Josiel Mendonça Soares de Souza** (developer)
- **Riccardo Sturani** (collaborator)
- **Miguel Quartin** (collaborator)

## License

**BSD 3-Clause License**

Copyright (c) 2026, Josiel Mendonça Soares de Souza

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
