=================================
GWDALI Software (version 0.1.3)
=================================

Software developed to perform parameter estimations of gravitational waves from compact objects coalescence (CBC) via Gaussian and Beyond-Gaussian approximation of GW likelihood. The Gaussian approximation is related to Fisher Matrix, from which it is direct to compute the covariance matrix by inverting the Fisher Matrix **[1]**. GWDALI also deals with the not-so-infrequent cases of Fisher Matrix with zero-determinant. The Beyond-Gaussian approach uses the `Derivative Approximation for LIkelihoods <https://arxiv.org/abs/1401.6892>`_ (DALI) algorithm proposed in **[2]** and applied to gravitational waves in **[3]**, whose model parameter uncertainties are estimated via Monte Carlo sampling but less costly than using the GW likelihood with no approximation.
Check our paper in `arXiv <https://arxiv.org/abs/2307.10154>`_.

.. figure:: ./logo_gwdali.png
   :alt: GWDALI logo
   :align: center

=================================
Installation
=================================

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
    detectors
    license
    author

=================================  
References
=================================

    **[1]** L. S. Finn and D. F. Chernoff, “Observing binary inspiral in gravitational radiation: One interferometer,” Phys. Rev. D, vol. 47, pp. 2198–2219, 1993.

    **[2]** E. Sellentin, M. Quartin, and L. Amendola, “Breaking the spell of gaussianity: forecasting with higher order fisher matrices,” Monthly Notices of the Royal Astronomical Society, vol. 441, no. 2, pp. 1831–1840, 2014.

    **[3]** Z. Wang, C. Liu, J. Zhao, and L. Shao, “Extending the fisher information matrix in gravitational-wave data analysis,” arXiv preprint arXiv:2203.02670, 2022.
