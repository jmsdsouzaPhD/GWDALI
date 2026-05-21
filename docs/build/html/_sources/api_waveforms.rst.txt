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