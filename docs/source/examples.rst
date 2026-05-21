Examples
========

The current version of GWDALI deals only whith waveforms in the frequency space.

***************************
Jupyter Notebooks (.ipynb)
***************************

* `1_get_hphx.ipynb <https://github.com/jmsdsouzaPhD/GWDALI/blob/main/Examples_py/1_get_hphx.py>`_
    - Compute polarizations :math:`h_+,\ h_{\times}`
* `2_DrawDetectors.ipynb <https://github.com/jmsdsouzaPhD/GWDALI/blob/main/Examples_py/1_get_hphx.py>`_
    - Vizualization of the GW detectors network choosen by the user
* `3_get_strains.ipynb <https://github.com/jmsdsouzaPhD/GWDALI/blob/main/Examples_py/1_get_hphx.py>`_
    - Compute detectors output :math:`h=F_+h_+ + F_{\times}h_{\times}`
* `4_get_derivatives.ipynb <https://github.com/jmsdsouzaPhD/GWDALI/blob/main/Examples_py/1_get_hphx.py>`_
    - Compute derivative of the GW signals (differention method choosen by the user: numerical/automatic differentiation)
        - Fisrt derivatives :math:`\partial_ih`
        - Second derivatives :math:`\partial_i\partial_jh`
        - Thirds derivatives :math:`\partial_i\partial_j\partial_kh`
* `5_get_tensors.ipynb <https://github.com/jmsdsouzaPhD/GWDALI/blob/main/Examples_py/1_get_hphx.py>`_
    - Compute DALI Tensors:
        - Fisher :math:`F_{ij}\equiv\langle\partial_i h|\partial_j h\rangle`
        - Doublet 
            * :math:`D_{ijk}\equiv\langle\partial_i h|\partial_{jk} h\rangle`
            * :math:`D_{ijkl}\equiv\langle\partial_{ij} h|\partial_{kl} h\rangle`
        - Triplet 
            * :math:`T_{ijkl}\equiv\langle\partial_i h|\partial_{jkl} h\rangle`
            * :math:`T_{ijklm}\equiv\langle\partial_{ij} h|\partial_{klm} h\rangle`
            * :math:`T_{ijklmn}\equiv\langle\partial_{ijk} h|\partial_{lmn} h\rangle`
* `6_GWDALI_Grid.ipynb <https://github.com/jmsdsouzaPhD/GWDALI/blob/main/Examples_py/1_get_hphx.py>`_
    - Compute GW Posteriors from a N-dimensional Grid (choosen by the user)
* `7_GWDALI_MCMC.ipynb <https://github.com/jmsdsouzaPhD/GWDALI/blob/main/Examples_py/1_get_hphx.py>`_
    - Sampling parameters for Posteriors estimation via MCMC or Nested Sampling methods


***************************
Top-Down Python Codes (.py)
***************************

* `1_get_hphx.py <https://github.com/jmsdsouzaPhD/GWDALI/blob/main/Examples_py/1_get_hphx.py>`_
* `2_get_strains.py <https://github.com/jmsdsouzaPhD/GWDALI/blob/main/Examples_py/1_get_hphx.py>`_
* `3_get_derivatives.py <https://github.com/jmsdsouzaPhD/GWDALI/blob/main/Examples_py/1_get_hphx.py>`_
* `4_get_dali-tensors.py <https://github.com/jmsdsouzaPhD/GWDALI/blob/main/Examples_py/1_get_hphx.py>`_
* `5_get_SNR.py <https://github.com/jmsdsouzaPhD/GWDALI/blob/main/Examples_py/1_get_hphx.py>`_
* `6_draw_detectors.py <https://github.com/jmsdsouzaPhD/GWDALI/blob/main/Examples_py/1_get_hphx.py>`_
* `7_get_map.py <https://github.com/jmsdsouzaPhD/GWDALI/blob/main/Examples_py/1_get_hphx.py>`_
* `8_check_priors.py <https://github.com/jmsdsouzaPhD/GWDALI/blob/main/Examples_py/1_get_hphx.py>`_
* `9_gwdali_grid.py <https://github.com/jmsdsouzaPhD/GWDALI/blob/main/Examples_py/1_get_hphx.py>`_
* `10_gwdali_mcmc.py <https://github.com/jmsdsouzaPhD/GWDALI/blob/main/Examples_py/1_get_hphx.py>`_