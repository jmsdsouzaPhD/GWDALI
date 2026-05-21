DALI Tensors
============

.. py:function:: GWDALI.get_dali_tensors(GwPrms, detectors, FreeParams, method, approx, enable_jax_waveforms=True, diff_method="autodiff", step_size=[1.e-6,1.e-4,1.e-2], hide_info=False, **kwargs)

   Computes Fisher, Doublet and Triplet tensors.

   :param dict GwPrms: GW source parameters.
   :param list of dicts detectors: Detector network.
   :param list FreeParams: Free parameters.
   :param str method: Tensor expansion method.
   :param str approx: Waveform approximant.
   :param str diff_method: Derivative method.
   :param float or list step_size: Numerical derivative step sizes.

   :returns:
      Dictionary containing ``Fisher``, ``Doublet`` and ``Triplet``.

   :rtype:
      dict
