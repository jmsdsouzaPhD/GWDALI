Main Interface
==============

.. py:function:: GWDALI.GWDALI(GwPrms, detectors, FreeParams, new_priors=None, approx="TaylorF2", method="Doublet", sampler="nestle", diff_method="autodiff", dali_tensors=None, step_size=[1.e-6], run_sampler=True, hide_info=False, npoints=300, nwalkers=None, ntemps=10, nburn=0., limits=None, pos0=None, npool=None, verbose=True, remove_out=True, output_name=None, save_bilby_path=True, bilby_path="outputs_bilby/", enable_jax_waveforms=True, **kwargs)

   Main interface of GWDALI v1.0.

   Computes:

   - Detector strains
   - Signal-to-noise ratios
   - Fisher Matrix
   - Doublet tensors
   - Triplet tensors
   - Posterior samples

   :param dictionary GwPrms: Dictionary containing GW source parameters.
   :param list detectors: Detector network configuration.
   :param list FreeParams: Free parameters used in the inference.
   :param dictionary new_priors: Redefine default priors.
   :param str approx: Waveform approximant.
   :param str method: DALI method [``"Fisher"``, ``"Doublet"``, ``"Triplet"``].
   :param str sampler: Posterior sampler backend from `bilby <https://lscsoft.docs.ligo.org/bilby/>`_.
   :param str diff_method: Derivative method [``"autodiff"``, ``"numdiff"``].
   :param list dali_tensors: User-defined DALI tensors (list of arrays).
   :param float or list step_size: Numerical derivative step sizes (float or list).
   :param bool run_sampler: Run posterior sampling.
   :param bool hide_info: Hide terminal outputs.
   :param int npoints: Number of posterior samples.
   :param int nwalkers: Number of walkers used in MCMC samplers.
   :param int ntemps: Number of temperatures in parallel tempering samplers.
   :param float nburn: Burn-in fraction.
   :param dictionary limits: Parameter limits used in the sampler.
   :param array pos0: Initial walker positions.
   :param int npool: Number of multiprocessing pools.
   :param bool verbose: Show sampler outputs.
   :param bool remove_out: Remove posterior outliers.
   :param str output_name: Output filename.
   :param bool save_bilby_path: Save bilby outputs.
   :param str bilby_path: Output directory for bilby files.
   :param bool enable_jax_waveforms: Enable JAX waveform backend.
   :param bool disable_jit: Disable ``jax.jit()``.
   :param bool EarthRotation: Enable Earth rotation corrections.
   :param bool jitgrad: Enable ``jax.jit()`` on derivatives.

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