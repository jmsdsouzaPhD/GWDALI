API Reference
=============

Main Interface
===============

.. py:function:: GWDALI.GWDALI(GwPrms, detectors, FreeParams, new_priors=None, approx="TaylorF2", method="Doublet", sampler="nestle", diff_method="autodiff", dali_tensors=None, step_size=[1.e-6], run_sampler=True, hide_info=False, npoints=300, nwalkers=None, ntemps=10, nburn=0., limits=None, pos0=None, npool=None, verbose=True, remove_out=True, output_name=None, save_bilby_path=True, bilby_path="outputs_bilby/", enable_jax_waveforms=True, **kwargs)

   Main interface of GWDALI v1.0.

   Computes:

   - Detector strains
   - Signal-to-noise ratios
   - Fisher Matrix
   - Doublet tensors
   - Triplet tensors
   - Posterior samples

   :param GwPrms: Dictionary containing GW source parameters.
   :param detectors: Detector network configuration.
   :param FreeParams: Free parameters used in the inference.
   :param new_priors: Redefine default priors.
   :param approx: Waveform approximant.
   :param method: DALI method [``"Fisher"``, ``"Doublet"``, ``"Triplet"``].
   :param sampler: Posterior sampler backend from `bilby <https://lscsoft.docs.ligo.org/bilby/>`_.
   :param diff_method: Derivative method [``"autodiff"``, ``"numdiff"``].
   :param dali_tensors: User-defined DALI tensors.
   :param step_size: Numerical derivative step sizes.
   :param run_sampler: Run posterior sampling.
   :param hide_info: Hide terminal outputs.
   :param npoints: Number of posterior samples.
   :param nwalkers: Number of walkers used in MCMC samplers.
   :param ntemps: Number of temperatures in parallel tempering samplers.
   :param nburn: Burn-in fraction.
   :param limits: Parameter limits used in the sampler.
   :param pos0: Initial walker positions.
   :param npool: Number of multiprocessing pools.
   :param verbose: Show sampler outputs.
   :param remove_out: Remove posterior outliers.
   :param output_name: Output filename.
   :param save_bilby_path: Save bilby outputs.
   :param bilby_path: Output directory for bilby files.
   :param enable_jax_waveforms: Enable JAX waveform backend.
   :param disable_jit: Disable ``jax.jit()``.
   :param EarthRotation: Enable Earth rotation corrections.
   :param jitgrad: Enable ``jax.jit()`` on derivatives.

   Supported approximants:

   - ``"TaylorF2"``
   - ``"TaylorF2_ISCO"``
   - ``"TaylorF2_Spinless"``
   - ``"TaylorF2_Spinless_0PN"``
   - ``"IMRPhenomA"``
   - ``"IMRPhenomB"``
   - ``"IMRPhenomC"``
   - ``"IMRPhenomD"``
   - ``"IMRPhenomHM"``

   :returns:
      ``Results, Truths, Tensors, Fisher, runtimes``

   :rtype:
      tuple


Waveforms
==========

.. py:function:: GWDALI.get_hphx(detectors, GwPrms, freq, approx, enable_jax_waveforms=True, **kwargs)

   Computes GW polarizations.

   :param detectors: Detector network.
   :param GwPrms: GW source parameters.
   :param freq: Frequency array.
   :param approx: Waveform approximant.

.. py:function:: GWDALI.get_strain(detectors, GwPrms, freq, approx, enable_jax_waveforms=True, **kwargs)

   Computes detector strains.

   :param detectors: Detector network.
   :param GwPrms: GW source parameters.
   :param freq: Frequency array.
   :param approx: Waveform approximant.

.. py:function:: GWDALI.get_SNR(detectors, GwPrms, approx, enable_jax_waveforms=True, **kwargs)

   Computes detector and network signal-to-noise ratios.


Derivatives
============

.. py:function:: GWDALI.get_derivatives(FreeParams, approx, GwPrms, dets, freq, diff_order="first", diff_method="numdiff", step_size=1.e-6, full_tensor=True, enable_jax_waveforms=True, **kwargs)

   Computes waveform derivatives.

   :param FreeParams: Parameters used in derivatives.
   :param approx: Waveform approximant.
   :param GwPrms: GW source parameters.
   :param dets: Detector network.
   :param freq: Frequency array.
   :param diff_order: Derivative order [``"first"``, ``"second"``, ``"third"``].
   :param diff_method: Derivative method [``"autodiff"``, ``"numdiff"``].
   :param step_size: Relative numerical derivative step size.
   :param full_tensor: Return full symmetric tensors.
   :param jitgrad: Enable ``jax.jit()`` on derivatives.
   :param EarthRotation: Enable Earth rotation corrections.

   :returns:
      ``Diff_values, time_diff``

   :rtype:
      tuple


DALI Tensors
=============

.. py:function:: GWDALI.get_dali_tensors(GwPrms, detectors, FreeParams, method, approx, enable_jax_waveforms=True, diff_method="autodiff", step_size=[1.e-6,1.e-4,1.e-2], hide_info=False, **kwargs)

   Computes Fisher, Doublet and Triplet tensors.

   :param GwPrms: GW source parameters.
   :param detectors: Detector network.
   :param FreeParams: Free parameters.
   :param method: Tensor expansion method.
   :param approx: Waveform approximant.
   :param diff_method: Derivative method.
   :param step_size: Numerical derivative step sizes.

   :returns:
      Dictionary containing ``Fisher``, ``Doublet`` and ``Triplet``.

   :rtype:
      dict


Priors
=======

.. py:function:: GWDALI.Priors(FreeParams, name=None, new_priors=None, plot=False)

   Builds Bayesian priors.

   :param FreeParams: Free parameters.
   :param name: Prior preset name.
   :param new_priors: User-defined priors.
   :param plot: Plot prior distributions.


Detector Utilities
===================

.. py:function:: GWDALI.get_map(detectors, plot_map=False)

   Computes detector antenna response maps.

.. py:function:: GWDALI.draw_detectors(detectors)

   Draws detector locations on Earth maps.