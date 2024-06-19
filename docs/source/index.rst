=================================
GWDALI Software (version 0.1.4)
=================================

Software developed to perform parameter estimations of gravitational waves from compact objects coalescence (CBC) via Gaussian and Beyond-Gaussian approximation of GW likelihood **[1]**. The Gaussian approximation is related to Fisher Matrix, from which it is direct to compute the covariance matrix by inverting the Fisher Matrix **[2]**. GWDALI also deals with the not-so-infrequent cases of Fisher Matrix with zero-determinant, for instance, from Fisher Matrix inversion, the uncertainties of the luminosity distance :math:`\sigma_{d_L}(\iota)` diverges for small values of source inclinations :math:`\iota` (in contrast to what is shown in **[3]**). The Beyond-Gaussian approach uses the `Derivative Approximation for LIkelihoods <https://arxiv.org/abs/1401.6892>`_ (DALI) algorithm proposed in **[4]** and applied to gravitational waves in **[5]**, whose model parameter uncertainties are estimated via Monte Carlo sampling but less costly than using the GW likelihood with no approximation.
Check our paper in `arXiv <https://arxiv.org/abs/2307.10154>`_.

* `Published Paper (Astronomy and Computing) <https://www.sciencedirect.com/science/article/abs/pii/S2213133723000744>`_
* `pypi page <https://pypi.org/project/gwdali/>`_
* `github page <https://github.com/jmsdsouzaPhD/GWDALI/>`_

*********************************
Installation
*********************************

To install the software run the command below:

.. code-block:: console

    $ pip install gwdali

.. toctree::
    :maxdepth: 1
    :caption: Contents:
    
    usage
    api
    dali
    priors
    GW_Sources
    detectors
    license
    author

*********************************  
References
*********************************

    **[1]** de Souza, J. M. S., & Sturani, R. (2023). GWDALI: A Fisher-matrix based software for gravitational wave parameter-estimation beyond Gaussian approximation. Astronomy and Computing, 45, 100759.

    **[2]** Finn, L. S., & Chernoff, D. F. (1993). Observing binary inspiral in gravitational radiation: One interferometer. Physical Review D, 47(6), 2198.

    **[3]** de Souza, J. M. S., & Sturani, R. (2023). Luminosity distance uncertainties from gravitational wave detections of binary neutron stars by third generation observatories. Physical Review D, 108(4), 043027.

    **[4]** Sellentin, E., Quartin, M., & Amendola, L. (2014). Breaking the spell of Gaussianity: forecasting with higher order Fisher matrices. Monthly Notices of the Royal Astronomical Society, 441(2), 1831-1840.

    **[5]** Wang, Z., Liu, C., Zhao, J., & Shao, L. (2022). Extending the Fisher information matrix in gravitational-wave data analysis. The Astrophysical Journal, 932(2), 102.
