# **GWDALI (v 1.0)**

A Fisher-Based Software for Parameter Estimation from Gravitational Waves

![Logo](https://github.com/jmsdsouzaPhD/GWDALI/blob/main/docs/source/_static/logo_gwdali.png)

Documentation in [https://gwdali.readthedocs.io/en/latest/](https://gwdali.readthedocs.io/en/latest/)

Main Paper: [arXiv:2307.10154](https://arxiv.org/abs/2307.10154) and [Astronomy and Computing, v. 45, p. 100759, 2023](https://www.sciencedirect.com/science/article/abs/pii/S2213133723000744)

Paper about implementation of automatic differentiation: [arXiv:2510.16955](https://arxiv.org/abs/2510.16955)

References: [arXiv:1401.6892](https://arxiv.org/abs/1401.6892) and [arXiv:2203.02670](https://arxiv.org/abs/2203.02670)

## Installation

Install GWDALI from PyPI with

```bash
pip install gwdali
```

This command also installs the required Python dependencies, including JAX.

## Additional Requirements

GWDALI depends on **JAX** for automatic differentiation, JIT-accelerated likelihood evaluations, and several internal numerical routines.

For most systems, `pip install gwdali` installs JAX automatically.

If the JAX installation through `pip` fails (for example on older CPUs without AVX support or on specific hardware platforms), install JAX manually following the official JAX installation instructions or using **conda-forge**, and then install GWDALI.

```bash
conda install -c conda-forge jax
```

If you intend to use **LAL waveform models**, install **LALSuite** as well:

```bash
conda install -c conda-forge lalsuite
conda install -c conda-forge lalsimulation
```

The waveform models implemented directly in GWDALI do **not** require LALSuite.

## Documentation

Available in [https://gwdali.readthedocs.io/en/latest/](https://gwdali.readthedocs.io/en/latest/)
    
## Functionalities

- **get_hphx()**: It returns plus/cross polarizations in the frequency space (SPA);
- **get_strain()**: It returns detector strains (signals) in the frequency space;
- **get_SNR()**: It returns detector-network signal-to-noise ratios (individuals and net);
- **draw_detectors()**: It returns a world map showing the chosen detector network configuration;
- **get_derivatives()**: It returns detector signal derivatives;
- **get_tensors()**: It returns DALI tensors including Fisher matrix;
- **Priors()**: Check/Visualize priors to be used in Posterior evaluations;
- **GWDALI()**: Returns MCMC samples, Fisher-inversion samples, or posterior grid arrays.

Check [https://gwdali.readthedocs.io/en/latest/examples.html](https://gwdali.readthedocs.io/en/latest/examples.html) for usage examples.

## References

**[1]** de Souza, J. M. S., & Sturani, R. (2023). GWDALI: A Fisher-matrix based software for gravitational wave parameter-estimation beyond Gaussian approximation. Astronomy and Computing, 45, 100759.

**[2]** de Souza, J. M. S., & Quartin, M. (2026). On the use of the Derivative Approximation for Likelihoods for gravitational wave inference. Journal of Cosmology and Astroparticle Physics, 2026(05), 101.

**[3]** Finn, L. S., & Chernoff, D. F. (1993). Observing binary inspiral in gravitational radiation: One interferometer. Physical Review D, 47(6), 2198.

**[4]** de Souza, J. M. S., & Sturani, R. (2023). Luminosity distance uncertainties from gravitational wave detections of binary neutron stars by third generation observatories. Physical Review D, 108(4), 043027.

**[5]** Sellentin, E., Quartin, M., & Amendola, L. (2014). Breaking the spell of Gaussianity: forecasting with higher order Fisher matrices. Monthly Notices of the Royal Astronomical Society, 441(2), 1831-1840.

**[6]** Wang, Z., Liu, C., Zhao, J., & Shao, L. (2022). Extending the Fisher information matrix in gravitational-wave data analysis. The Astrophysical Journal, 932(2), 102.

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

## Aknowledgements

This work was partialy suported by:
- Coordenação de Aperfeiçoamente de Pessoal de Ensino Superior (CAPES);
- Fundação Carlos Chargas Filho de Amparo à Pesquisa do Estado do Rio de Janeiro (FAPERJ);
- Fundação de Amparo à Pesquisa e Inovação do Espírito Santo (FAPES);

The authors also thank Davi Rodrigues (UFES) for usefull discussions.
