Waveforms
=========

.. py:function:: GWDALI.get_hphx(detectors, GwPrms, freq, approx, enable_jax_waveforms=True, **kwargs)

   Computes GW polarizations.

   :param list of dicts detectors: Detector network.
   :param dict GwPrms: GW source parameters.
   :param array freq: Frequency array.
   :param str approx: Waveform approximant.

.. py:function:: GWDALI.get_strain(detectors, GwPrms, freq, approx, enable_jax_waveforms=True, **kwargs)

   Computes detector strains.

   :param list of dicts detectors: Detector network.
   :param dict GwPrms: GW source parameters.
   :param array freq: Frequency array.
   :param str approx: Waveform approximant.

.. py:function:: GWDALI.get_SNR(detectors, GwPrms, approx, enable_jax_waveforms=True, **kwargs)

   Computes detector and network signal-to-noise ratios.

.. warning:: 

   The available jax-friendly approximants are:
   - ``"TaylorF2"``
   - ``"TaylorF2_ISCO"``
   - ``"TaylorF2_Spinless"``
   - ``"TaylorF2_Spinless_0PN"``
   - ``"IMRPhenomA"``
   - ``"IMRPhenomB"``
   - ``"IMRPhenomC"``
   - ``"IMRPhenomD"``
   - ``"IMRPhenomHM"``

   Otherwise the software will use the LALSuite approximants via `lalsimulation <https://lscsoft.docs.ligo.org/lalsuite/lalsimulation/>`_.

See the examples bellow also to check how close are the JAX waveforms implementations from the implementations in LAL.

* `Example (Get Polarizations) <https://github.com/jmsdsouzaPhD/GWDALI/blob/main/Examples_jupyter/1_get_hphx.ipynb>`_
* `Example (Get Strains) <https://github.com/jmsdsouzaPhD/GWDALI/blob/main/Examples_jupyter/3_get_strains.ipynb>`_
