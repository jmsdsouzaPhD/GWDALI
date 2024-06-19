=================================
GWDALI Software (version 0.1.4)
=================================

Software developed to perform parameter estimations of gravitational waves from compact objects coalescence (CBC) via Gaussian and Beyond-Gaussian approximation of GW likelihood. The Gaussian approximation is related to Fisher Matrix, from which it is direct to compute the covariance matrix by inverting the Fisher Matrix **[1]**. GWDALI also deals with the not-so-infrequent cases of Fisher Matrix with zero-determinant, for instance, from Fisher Matrix inversion, the uncertainties of the luminosity distance :math:`\sigma_{d_L}(\iota)` diverges for small values of source inclinations :math:`\iota` (in contrast to what is shown in **[2]**). The Beyond-Gaussian approach uses the `Derivative Approximation for LIkelihoods <https://arxiv.org/abs/1401.6892>`_ (DALI) algorithm proposed in **[3]** and applied to gravitational waves in **[4]**, whose model parameter uncertainties are estimated via Monte Carlo sampling but less costly than using the GW likelihood with no approximation.
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

    **[1]** L. S. Finn and D. F. Chernoff, “Observing binary inspiral in gravitational radiation: One interferometer,” Phys. Rev. D, vol. 47, pp. 2198–2219, 1993.

	**[2]** de Souza, Josiel Mendonça Soares, and Riccardo Sturani. "Luminosity distance uncertainties from gravitational wave detections of binary neutron stars by third generation observatories." Physical Review D 108.4 (2023): 043027.

    **[3]** E. Sellentin, M. Quartin, and L. Amendola, “Breaking the spell of gaussianity: forecasting with higher order fisher matrices,” Monthly Notices of the Royal Astronomical Society, vol. 441, no. 2, pp. 1831–1840, 2014.

    **[4]** Z. Wang, C. Liu, J. Zhao, and L. Shao, “Extending the fisher information matrix in gravitational-wave data analysis,” arXiv preprint arXiv:2203.02670, 2022.
