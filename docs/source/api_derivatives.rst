Derivatives
===========

.. py:function:: GWDALI.get_derivatives(FreeParams, approx, GwPrms, detectors, freq, diff_order="first", diff_method="numdiff", step_size=1.e-6, full_tensor=True, enable_jax_waveforms=True, **kwargs)

   Computes waveform derivatives.

   :param list FreeParams: Parameters used in derivatives.
   :param str approx: Waveform approximant.
   :param dict GwPrms: GW source parameters.
   :param list of dicts detectors: Detector network.
   :param array freq: Frequency array.
   :param str diff_order: Derivative order [``"first"``, ``"second"``, ``"third"``].
   :param str diff_method: Derivative method [``"autodiff"``, ``"numdiff"``].
   :param float or list step_size: Relative numerical derivative step size.
   :param bool full_tensor: Return full symmetric tensors.
   :param bool jitgrad: Enable ``jax.jit()`` on derivatives.
   :param bool EarthRotation: Enable Earth rotation corrections.

   :returns:
      ``Diff_values, time_diff``

   :rtype:
      tuple

* `Example (Get Derivatives) <https://github.com/jmsdsouzaPhD/GWDALI/blob/main/Examples_jupyter/4_get_derivatives.ipynb>`_
